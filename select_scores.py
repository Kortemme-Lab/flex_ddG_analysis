import pandas as pd
import argparse
import getpass
import MySQLdb
import os
import sys
import time
import numpy as np
import json

from run_descriptions import all_runs

# Comment out structure orders below to make selection faster
structure_orders = [
    'id', # Order structures for averaging based on the random order they were created
    'WildTypeComplex', # Order structures for averaging by the wild type complex score
    # 'MutantComplex', # Order structures for averaging by the mutant type complex score
]

dataframe_cache = '/dbscratch/kyleb/new_query_cache'
desired_summed_and_averaged_columns = [
    'PredictionRunName', 'DataSetID', 'PredictionID', 'ScoreMethodID', 'total', 'ExperimentalDDG', 'StructureOrder',
]

global_score_columns = [
    'total',
    'fa_atr', 'fa_dun', 'fa_elec', 'fa_intra_rep', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb', 'omega', 'p_aa_pp', 'pro_close', 'rama', 'ref', 'yhh_planarity', 'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'fa_intra_atr_xover4', 'fa_intra_elec', 'fa_intra_rep_xover4', 'fa_intra_sol_xover4', 'hxl_tors', 'lk_ball', 'lk_ball_bridge', 'lk_ball_bridge_uncpl', 'lk_ball_iso', 'rama_prepro', 'lk_ball_wtd',
    'cart_bonded',
]

id_columns = [
    'PredictionID', 'DataSetID', 'PDBFileID', 'Subset', 'ExperimentalDDG',
    'PredictionRunName',
    'StructureOrder', # StructureOrder added by this script
    # Last three IDs are special, as they are used for grouping
    'ScoreMethodID', 'StructureID', 'ScoreType',
]

# Coerce the dtype for these
string_cols = ['PredictionRunName', 'Subset', 'PDBFileID', 'ScoreType', 'Mutations', 'AvgRuntime', 'MaxMemGB']
int_cols = ['PredictionID', 'DataSetID', 'StructureID', 'ScoreMethodID']
float_cols = ['ExperimentalDDG', 'total', 'talaris_total_check', 'ref_total_check', 'fa_atr', 'fa_dun', 'fa_elec', 'fa_intra_rep', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb', 'omega', 'p_aa_pp', 'pro_close', 'rama', 'ref', 'yhh_planarity', 'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'fa_intra_atr_xover4', 'fa_intra_elec', 'fa_intra_rep_xover4', 'fa_intra_sol_xover4', 'hxl_tors', 'lk_ball', 'lk_ball_bridge', 'lk_ball_bridge_uncpl', 'lk_ball_iso', 'rama_prepro', 'cart_bonded', 'lk_ball_wtd']

amino_acid_details = {
    'A': {'Polarity': 'H', 'van_der_Waals_volume': 67.0,  'Name': 'Alanine',       'LongCode': 'ALA', 'Aromaticity': 'L', 'Size': 'small'},
    'C': {'Polarity': 'P', 'van_der_Waals_volume': 86.0,  'Name': 'Cysteine',      'LongCode': 'CYS', 'Aromaticity': '-', 'Size': 'small'},
    'E': {'Polarity': 'C', 'van_der_Waals_volume': 109.0, 'Name': 'Glutamic acid', 'LongCode': 'GLU', 'Aromaticity': '-', 'Size': 'large'},
    'D': {'Polarity': 'C', 'van_der_Waals_volume': 91.0,  'Name': 'Aspartic acid', 'LongCode': 'ASP', 'Aromaticity': '-', 'Size': 'small'},
    'G': {'Polarity': 'P', 'van_der_Waals_volume': 48.0,  'Name': 'Glycine',       'LongCode': 'GLY', 'Aromaticity': '-', 'Size': 'small'},
    'F': {'Polarity': 'H', 'van_der_Waals_volume': 135.0, 'Name': 'Phenylalanine', 'LongCode': 'PHE', 'Aromaticity': 'R', 'Size': 'large'},
    'I': {'Polarity': 'H', 'van_der_Waals_volume': 124.0, 'Name': 'Isoleucine',    'LongCode': 'ILE', 'Aromaticity': 'L', 'Size': 'large'},
    'H': {'Polarity': 'C', 'van_der_Waals_volume': 118.0, 'Name': 'Histidine',     'LongCode': 'HIS', 'Aromaticity': '-', 'Size': 'large'},
    'K': {'Polarity': 'C', 'van_der_Waals_volume': 135.0, 'Name': 'Lysine',        'LongCode': 'LYS', 'Aromaticity': '-', 'Size': 'large'},
    'M': {'Polarity': 'H', 'van_der_Waals_volume': 124.0, 'Name': 'Methionine',    'LongCode': 'MET', 'Aromaticity': 'L', 'Size': 'large'},
    'L': {'Polarity': 'H', 'van_der_Waals_volume': 124.0, 'Name': 'Leucine',       'LongCode': 'LEU', 'Aromaticity': 'L', 'Size': 'large'},
    'N': {'Polarity': 'P', 'van_der_Waals_volume': 96.0,  'Name': 'Asparagine',    'LongCode': 'ASN', 'Aromaticity': '-', 'Size': 'small'},
    'Q': {'Polarity': 'P', 'van_der_Waals_volume': 114.0, 'Name': 'Glutamine',     'LongCode': 'GLN', 'Aromaticity': '-', 'Size': 'large'},
    'P': {'Polarity': 'H', 'van_der_Waals_volume': 90.0,  'Name': 'Proline',       'LongCode': 'PRO', 'Aromaticity': '-', 'Size': 'small'},
    'S': {'Polarity': 'P', 'van_der_Waals_volume': 73.0,  'Name': 'Serine',        'LongCode': 'SER', 'Aromaticity': '-', 'Size': 'small'},
    'R': {'Polarity': 'C', 'van_der_Waals_volume': 148.0, 'Name': 'Arginine',      'LongCode': 'ARG', 'Aromaticity': '-', 'Size': 'large'},
    'T': {'Polarity': 'P', 'van_der_Waals_volume': 93.0,  'Name': 'Threonine',     'LongCode': 'THR', 'Aromaticity': '-', 'Size': 'small'},
    'W': {'Polarity': 'H', 'van_der_Waals_volume': 163.0, 'Name': 'Tryptophan',    'LongCode': 'TRP', 'Aromaticity': 'R', 'Size': 'large'},
    'V': {'Polarity': 'H', 'van_der_Waals_volume': 105.0, 'Name': 'Valine',        'LongCode': 'VAL', 'Aromaticity': 'L', 'Size': 'small'},
    'Y': {'Polarity': 'H', 'van_der_Waals_volume': 141.0, 'Name': 'Tyrosine',      'LongCode': 'TYR', 'Aromaticity': 'R', 'Size': 'large'},
}

dtypes = {}
for col in string_cols:
    dtypes[col] = object
for col in int_cols:
    dtypes[col] = np.int64
for col in float_cols:
    dtypes[col] = np.float64

if not os.path.isdir( dataframe_cache ):
    os.makedirs( dataframe_cache )

def fetch_from_db_and_reorder( mysql_con ):
    with open('scores_select.sql', 'r') as f:
        all_scores_query = ' '.join( [ line.strip() for line in f.readlines() if not line.startswith('#') ] )

    # Load data and determine common DataSetIDs
    for benchmark_run in all_runs:
        step_multiplication_factor = benchmark_run.step_multiplication_factor
        prediction_set_name = benchmark_run.prediction_set_name
        print '\n%s' % str(prediction_set_name)

        cached_df_path = os.path.join( dataframe_cache, prediction_set_name + '_all-scores.csv.gz' )
        if os.path.isfile(cached_df_path):
            print 'Cached main df exists:', cached_df_path
            df = pd.read_csv(cached_df_path, compression = 'gzip' )
        else:
            df = pd.read_sql_query( all_scores_query % prediction_set_name, mysql_con)
            # Multiply by step factor
            if step_multiplication_factor:
                df.loc[:,'ScoreMethodID'] *= step_multiplication_factor

            df.to_csv(cached_df_path, compression = 'gzip' )

        print 'PDBs in set: %d' % len( df[['PDBFileID']].drop_duplicates() )
        print 'Datapoints length (without duplicates): %d' % len( df[['DataSetID']].drop_duplicates() )

        # Reorder structure ID's, if specified
        for structure_order in structure_orders:
            # Check for cache
            cached_query_df_path = os.path.join( dataframe_cache, prediction_set_name + '-' + structure_order + '.csv.gz' )
            if os.path.isfile( cached_query_df_path ):
                print 'Cached ordered df exists:', cached_query_df_path
            else:
                if structure_order == 'id':
                    new_query_df = df.assign( StructureOrder = 'id' )
                elif structure_order in df['ScoreType'].drop_duplicates().values:
                    new_df = df.copy()
                    new_df.set_index(['PredictionID', 'ScoreMethodID', 'StructureID'], inplace=True)
                    sorted_df = df.loc[ df['ScoreType'] == structure_order ].sort_values(
                        by = ['PredictionID', 'ScoreMethodID', 'total'],
                    )
                    new_series_values = {
                        'NewStructureID' : [],
                        'PredictionID' : [],
                        'ScoreMethodID' : [],
                        'StructureID' : [],
                    }
                    row_count = 1
                    last_id_tup = None
                    last_prediction_id = None
                    last_score_method_id = None
                    for index, row in sorted_df.iterrows():
                        row_id_tup = (row['PredictionID'], row['ScoreMethodID'])
                        if last_id_tup == None:
                            last_id_tup = row_id_tup
                        if last_id_tup != row_id_tup:
                            row_count = 1
                            last_id_tup = row_id_tup
                        new_series_values['NewStructureID'].append( row_count )
                        new_series_values['PredictionID'].append( row['PredictionID'] )
                        new_series_values['ScoreMethodID'].append( row['ScoreMethodID'] )
                        new_series_values['StructureID'].append( row['StructureID'] )
                        row_count += 1
                    new_column_df = pd.DataFrame( new_series_values )
                    new_column_df.set_index(['PredictionID', 'ScoreMethodID', 'StructureID'], inplace=True)
                    new_df = new_df.merge(new_column_df, how = 'left', left_index = True, right_index = True )
                    new_query_df = new_df.reset_index().assign( StructureOrder = structure_order )
                    new_query_df.drop('StructureID', axis=1, inplace=True)
                    new_query_df.rename( columns = {'NewStructureID' : 'StructureID'}, inplace = True )
                new_query_df.to_csv( cached_query_df_path, compression = 'gzip' )

    mysql_con.close()

def sum_and_average():
    df_paths = [ f for f in os.listdir(dataframe_cache) if '.csv' in f ]
    paths_to_analyze = []
    for benchmark_run in all_runs:
        for structure_order in structure_orders:
            matching_path = None
            for df_path in df_paths:
                if df_path.startswith( benchmark_run.prediction_set_name + '-' + structure_order + '.csv' ):
                    if matching_path != None:
                        print
                        print matching_path
                        print df_path
                        print benchmark_run.prediction_set_name, structure_order
                        print
                    assert( matching_path == None )
                    matching_path = df_path
            if matching_path == None:
                print
                print benchmark_run.prediction_set_name + '-' + structure_order + '.csv'
                print benchmark_run.prediction_set_name, structure_order
                print
            assert( matching_path != None )
            paths_to_analyze.append( (benchmark_run, structure_order, os.path.join( dataframe_cache, matching_path ) ) )

    output_dir = os.path.join( dataframe_cache, 'summed_and_averaged' )
    if not os.path.isdir( output_dir ):
        os.makedirs( output_dir )
    partial_output_dir = os.path.join( dataframe_cache, 'partial_summed_and_averaged' )
    if not os.path.isdir( partial_output_dir ):
        os.makedirs( partial_output_dir )

    for benchmark_run, structure_order, df_path in paths_to_analyze:
        print 'Summing and averaging:', df_path
        df = pd.read_csv(
            df_path, compression = 'gzip',
            dtype = dtypes,
        )

        this_score_columns = [col for col in global_score_columns if col in df.columns]
        score_and_id_columns = this_score_columns + id_columns
        complex_scores_df = df[df.ScoreType.isin( ['MutantLPartner', 'MutantRPartner', 'WildTypeComplex'] )][score_and_id_columns]
        for score_column in this_score_columns:
            complex_scores_df.loc[:,score_column] *= -1
        complex_scores_df = complex_scores_df.append(
            df[df.ScoreType.isin( ['DDG', 'WildTypeLPartner', 'WildTypeRPartner', 'MutantComplex'] )][score_and_id_columns]
        )
        complex_scores_df = complex_scores_df.sort_index() #( columns = ['PredictionID', 'StructureID', ] )
        summed_df_noindex = complex_scores_df.groupby(id_columns[:-1], as_index = False)[this_score_columns].sum()
        # Check that each group of ScoreTypes has no null values
        desired_group_length = len( df['ScoreType'].drop_duplicates() )
        for name, group in complex_scores_df.groupby(id_columns[:-1]):
            if len(group) != desired_group_length:
                print 'len(group) != %d' % desired_group_length
                print name
                print group
                print
        summed_df = complex_scores_df.groupby(id_columns[:-1])[this_score_columns].sum().reset_index()

        df_structure_orders = summed_df['StructureOrder'].drop_duplicates()
        if len(df_structure_orders) != 1:
            print df_structure_orders
            assert( len(df_structure_orders) == 1 )
        structure_order = df_structure_orders[0]
        max_struct_id = summed_df['StructureID'].max()

        for num_structs in xrange( 1, max_struct_id + 1 ):
            new_structure_order = structure_order + '_%02d' % num_structs
            csv_path = os.path.join( output_dir, '%s-%s.csv.gz' % (benchmark_run.prediction_set_name, new_structure_order) )
            partial_csv_path = os.path.join( partial_output_dir, '%s-%s-%s-partial.csv.gz' % (benchmark_run.prediction_set_name, new_structure_order, '%05d' ) )

            subset_summed_df = summed_df.loc[ summed_df['StructureID'] <= num_structs ]
            avg_df = subset_summed_df.groupby(id_columns[:-2])[this_score_columns].mean().round(decimals=4).reset_index()
            avg_df = avg_df.assign( StructureOrder = new_structure_order )

            # Write csv with partial score terms
            partial_columns = list( set( desired_summed_and_averaged_columns ).union( set(global_score_columns) ) )
            for backrub_step in sorted( avg_df['ScoreMethodID'].drop_duplicates() ):
                if backrub_step < 100 or int(backrub_step) % 5000 == 0 and ( (num_structs == 1) or (num_structs % 10 == 0) ):
                    inner_partial_csv_path = partial_csv_path % int(backrub_step)
                    if os.path.isfile( inner_partial_csv_path ):
                        continue
                    partial_avg_df = avg_df.loc[ avg_df['ScoreMethodID'] == backrub_step ][partial_columns].round(4)
                    partial_avg_df.sort_values(
                        ['PredictionRunName', 'DataSetID', 'ScoreMethodID'],
                        inplace = True,
                    )
                    partial_avg_df.to_csv( inner_partial_csv_path, compression = 'gzip' )
                    print 'Saved:', inner_partial_csv_path

            if os.path.isfile( csv_path ):
                continue

            # Write csv with only total score terms
            avg_df = avg_df[desired_summed_and_averaged_columns].round(4)
            avg_df.sort_values(
                ['PredictionRunName', 'DataSetID', 'ScoreMethodID'],
                inplace = True,
            )

            avg_df.to_csv( csv_path, compression = 'gzip' )
            print 'Saved:', csv_path

def fetch_zemu_properties( mysql_con, print_debug = False ):
    with open('dataset_select.sql', 'r') as f:
        dataset_query = ' '.join( [ line.strip() for line in f.readlines() if not line.startswith('#') ] )

    df = pd.read_sql_query( dataset_query, mysql_con)

    antibody_pdbs = [ x.upper() for x in pd.read_csv( 'antibodies.csv', sep = '\t' )['pdb'].values ]

    single = set()
    multiple = set()
    multiple_all_ala = set()
    multiple_none_ala = set()
    all_ala = set()
    sing_ala = set()
    all_s2l = set()
    all_l2s = set()
    some_s2l = set()
    some_l2s = set()

    res_gte25 = set()
    res_gt15_lt25 = set()
    res_lte15 = set()

    antibodies = set()

    by_structure = {}
    by_mutant_type = {}

    for index, row in df.iterrows():
        mutations = [ x.split()[1] for x in row['Mutations'].split(';') ] # Split, and throw away chains
        dataset_id = row['DataSetID']

        if row['PDBFileID'] in antibody_pdbs:
            antibodies.add( dataset_id )

        if row['PDBFileID'] not in by_structure:
            by_structure[ row['PDBFileID'] ] = set()
        by_structure[ row['PDBFileID'] ].add( dataset_id )

        if row['Resolution'] <= 1.5:
            res_lte15.add( dataset_id )
        elif row['Resolution'] < 2.5:
            res_gt15_lt25.add( dataset_id )
        elif row['Resolution'] >= 2.5:
            res_gte25.add( dataset_id )

        if len(mutations) == 1:
            single.add( dataset_id )
        elif len(mutations) > 1:
            multiple.add( dataset_id )

        mutants_all_ala = True
        mutants_any_ala = False

        mutants_some_s2l = False
        mutants_all_s2l = True

        mutants_some_l2s = False
        mutants_all_l2s = True

        mutant_all_aa = None

        for mutation in mutations:
            if mutants_all_ala and mutation[-1] != 'A':
                mutants_all_ala = False

            if not mutants_any_ala and mutation[-1] == 'A':
                mutants_any_ala = True

            if  amino_acid_details[mutation[0]]['van_der_Waals_volume'] < amino_acid_details[mutation[-1]]['van_der_Waals_volume']:
                mutants_some_s2l = True
                mutants_all_l2s = False
            else:
                mutants_all_s2l = False
                mutants_some_l2s = True

            if mutant_all_aa == None:
                mutant_all_aa = mutation[-1]
            elif mutant_all_aa != mutation[-1]:
                mutant_all_aa = False

        if mutant_all_aa:
            if mutant_all_aa not in by_mutant_type:
                by_mutant_type[mutant_all_aa] = set()
            by_mutant_type[mutant_all_aa].add( dataset_id )

        if not mutants_any_ala and len(mutations) > 1:
            multiple_none_ala.add( dataset_id )

        if mutants_all_ala:
            all_ala.add( dataset_id )
            if len(mutations) > 1:
                multiple_all_ala.add( dataset_id )
            elif len(mutations) == 1:
                sing_ala.add( dataset_id )

        if mutants_all_s2l:
            all_s2l.add( dataset_id )

        if mutants_some_s2l:
            some_s2l.add( dataset_id )

        if mutants_all_l2s:
            all_l2s.add( dataset_id )

        if mutants_some_l2s:
            some_l2s.add( dataset_id )

    subsets_dict = {}

    # Load DSSP information
    assert( os.path.isfile( 'dssp_results.csv' ) ) # Generate with dssp.py if doesn't exist
    dssp = pd.read_csv( 'dssp_results.csv' )
    subsets_dict[ 'SS-all-None' ] = sorted( dssp.loc[ dssp['SS'].isnull() ]['DataSetID'] )
    for ss_type in sorted( dssp['SS'].dropna().drop_duplicates() ):
        ss_ids = sorted( dssp.loc[ dssp['SS'] == ss_type ]['DataSetID'] )
        if len(ss_ids) >= 5:
            subsets_dict[ 'SS-' + ss_type] = ss_ids
    for exposure_metric in [x.split('_')[0] for x in dssp.columns if x.endswith('exposure')]:
        for cut_off in [ 0.05, 0.1, 0.25 ]:
            buried_ids = dssp.loc[ dssp[exposure_metric + '_exposure'] < cut_off ]
            if len(buried_ids) > 0:
                subsets_dict[ 'buried_lt' + '%.2f' % cut_off + '_' + exposure_metric ] = sorted( buried_ids['DataSetID'] )

            buried_ids = dssp.loc[ dssp[exposure_metric + '_exposure'] >= cut_off ]
            if len(buried_ids) > 0:
                subsets_dict[ 'buried_gte' + '%.2f' % cut_off + '_' + exposure_metric ] = sorted( buried_ids['DataSetID'] )

    # for key, value in subsets_dict.iteritems():
    #     print key, len(value)

    if print_debug:
        print 'Summary:'
        print 'Single mutations:', len(single)
        print 'Multiple mutations:', len(multiple)
        print 'Multiple (none ala):', len(multiple_none_ala)
        print 'Multiple (all ala):', len(multiple_all_ala)
        print 'All alanines:', len(all_ala)
        print 'Single alanines:', len(sing_ala)
        print 'All small to large:', len(all_s2l)
        print 'All large to small:', len(all_l2s)
        print 'Some small to large:', len(some_s2l)
        print 'Some large to small:', len(some_l2s)
        print 'High res (res_lte15):', len(res_lte15)
        print 'Med res (res_gt15_lt25):', len(res_gt15_lt25)
        print 'Low res(res_gte25):', len(res_gte25)
        print 'Antibodies:', len(antibodies)
        print

    subsets_dict['mult_mut'] = sorted(multiple)
    subsets_dict['mult_all_ala'] = sorted(multiple_all_ala)
    subsets_dict['mult_none_ala'] = sorted(multiple_none_ala)
    subsets_dict['sing_mut'] = sorted(single)
    subsets_dict['s2l'] = sorted(all_s2l)
    subsets_dict['some_s2l'] = sorted(some_s2l)
    subsets_dict['ala'] = sorted(all_ala)
    subsets_dict['sing_ala'] = sorted(sing_ala)
    subsets_dict['l2s'] = sorted(all_l2s)
    subsets_dict['some_l2s'] = sorted(some_l2s)
    subsets_dict['res_lte15'] = sorted(res_lte15)
    subsets_dict['res_gt15_lt25'] = sorted(res_gt15_lt25)
    subsets_dict['res_gte25'] = sorted(res_gte25)
    subsets_dict['antibodies'] = sorted(antibodies)
    for mutant_all_aa in by_mutant_type:
        subsets_dict[ 'all-mut-to-' + mutant_all_aa ] = sorted(by_mutant_type[mutant_all_aa])
    subsets_dict
    for pdb in sorted(by_structure.keys()):
        pdb_dataset_ids = sorted(by_structure[pdb])
        if len(pdb_dataset_ids) >= 5:
            subsets_dict[ 'pdb-' + pdb ] = pdb_dataset_ids

    with open('subsets.json', 'w') as f:
        json.dump(subsets_dict, f, sort_keys = True, indent = 2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--host', default = 'localhost',
                        help = 'MySQL host name')
    parser.add_argument('--port', type = int, default = 3306,
                        help = 'MySQL port')
    parser.add_argument('--user', default = getpass.getuser(),
                        help = 'MySQL user')
    parser.add_argument('--password', default = '',
                        help = 'MySQL password')
    parser.add_argument('--db', default = 'DDG',
                        help = 'MySQL database schema')
    args = parser.parse_args()

    mysql_con = MySQLdb.connect( host = args.host,
                                 port = args.port,
                                 user = args.user, passwd = args.password,
                                 db = args.db
    )

    fetch_zemu_properties( mysql_con )
    fetch_from_db_and_reorder( mysql_con )
    sum_and_average()
