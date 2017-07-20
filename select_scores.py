import pandas as pd
import argparse
import getpass
import MySQLdb
import os
import sys
import time

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
    'fa_atr', 'fa_dun', 'fa_elec', 'fa_intra_rep', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb', 'hbond_sc', 'hbond_sr_bb', 'omega', 'p_aa_pp', 'pro_close', 'rama', 'ref', 'yhh_planarity', 'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'fa_intra_atr_xover4', 'fa_intra_elec', 'fa_intra_rep_xover4', 'fa_intra_sol_xover4', 'hxl_tors', 'lk_ball', 'lk_ball_bridge', 'lk_ball_bridge_uncpl', 'lk_ball_iso', 'rama_prepro',
    'cart_bonded',
]

id_columns = [
    'PredictionID', 'DataSetID', 'PDBFileID', 'Subset', 'ExperimentalDDG',
    'PredictionRunName',
    'StructureOrder', # StructureOrder added by this script
    # Last three IDs are special, as they are used for grouping
    'ScoreMethodID', 'StructureID', 'ScoreType',
]

# Coerce the dtype for these into strings
string_cols = ['PredictionRunName', 'Subset', 'PDBFileID', 'ScoreType', 'Mutations']

if not os.path.isdir( dataframe_cache ):
    os.makedirs( dataframe_cache )

def fetch_from_db_and_reorder( mysql_host, mysql_port, mysql_user, mysql_pass, mysql_db ):
    with open('scores_select.sql', 'r') as f:
        all_scores_query = ' '.join( [ line.strip() for line in f.readlines() if not line.startswith('#') ] )

    mysql_con = MySQLdb.connect( host = mysql_host,
                                 port = mysql_port,
                                 user = mysql_user, passwd = mysql_pass,
                                 db = mysql_db
    )

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
                if benchmark_run.prediction_set_name in df_path and structure_order + '.csv' in df_path:
                    assert( matching_path == None )
                    matching_path = df_path
            assert( matching_path != None )
            paths_to_analyze.append( (benchmark_run, structure_order, os.path.join( dataframe_cache, matching_path ) ) )

    output_dir = os.path.join( dataframe_cache, 'summed_and_averaged' )
    if not os.path.isdir( output_dir ):
        os.makedirs( output_dir )

    for benchmark_run, structure_order, df_path in paths_to_analyze:
        print 'Summing and averaging:', df_path
        df = pd.read_csv( df_path, compression = 'gzip' )

        this_score_columns = [col for col in global_score_columns if col in df.columns]
        score_and_id_columns = this_score_columns + id_columns
        complex_scores_df = df[df.ScoreType.isin( ['MutantLPartner', 'MutantRPartner', 'WildTypeComplex'] )][score_and_id_columns]
        for score_column in this_score_columns:
            complex_scores_df.loc[:,score_column] *= -1
        complex_scores_df = complex_scores_df.append(
            df[df.ScoreType.isin( ['WildTypeLPartner', 'WildTypeRPartner', 'MutantComplex'] )][score_and_id_columns]
        )
        complex_scores_df = complex_scores_df.sort_index() #( columns = ['PredictionID', 'StructureID', ] )
        summed_df_noindex = complex_scores_df.groupby(id_columns[:-1], as_index = False)[this_score_columns].sum()
        for name, group in complex_scores_df.groupby(id_columns[:-1]):
            if len(group) != 6:
                print 'len(group) != 6'
                print name
                print group
                print
        summed_df = complex_scores_df.groupby(id_columns[:-1])[this_score_columns].sum().reset_index()

        df_structure_orders = summed_df['StructureOrder'].drop_duplicates()
        assert( len(df_structure_orders) == 1 )
        structure_order = df_structure_orders[0]
        max_struct_id = summed_df['StructureID'].max()

        for num_structs in xrange( 1, max_struct_id + 1 ):
            new_structure_order = structure_order + '_%02d' % num_structs
            csv_path = os.path.join( output_dir, '%s-%s.csv.gz' % (benchmark_run.prediction_set_name, new_structure_order) )
            if os.path.isfile( csv_path ):
                continue

            subset_summed_df = summed_df.loc[ summed_df['StructureID'] <= num_structs ]
            avg_df = subset_summed_df.groupby(id_columns[:-2])[this_score_columns].mean().round(decimals=4).reset_index()
            avg_df = avg_df.assign( StructureOrder = new_structure_order )
            avg_df = avg_df[desired_summed_and_averaged_columns].round(4)
            avg_df.sort_values(
                ['PredictionRunName', 'DataSetID', 'ScoreMethodID'],
                inplace = True,
            )

            avg_df.to_csv( csv_path, compression = 'gzip' )
            print 'Saved:', csv_path

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
    # fetch_from_db_and_reorder( args.host, args.port, args.user, args.password, args.db)
    sum_and_average()
