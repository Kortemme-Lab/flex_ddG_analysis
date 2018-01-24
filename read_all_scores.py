import cPickle as pickle
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import os
import sys
import shutil
from stats import fraction_correct, mae as calc_mae, bootstrap_xy_stat
import scipy
import json
import subprocess
import collections
import copy
import string
import random
import re
import copy

csv_paths = [
    'data/zemu_1.2-60000_rscript_validated-ref-id_01.csv.gz',
    'data/zemu_1.2-60000_rscript_validated-ref-id_30.csv.gz',
    'data/zemu_1.2-60000_rscript_validated-ref-id_50.csv.gz',
    'data/zemu_1.2-60000_rscript_validated-ref-WildTypeComplex_01.csv.gz',
    'data/zemu_1.2-60000_rscript_validated-ref-WildTypeComplex_30.csv.gz',
    'data/zemu_1.2-60000_rscript_validated-ref-WildTypeComplex_50.csv.gz',
    'data/zemu_control-69aa526-id_01.csv.gz',
    'data/zemu_control-69aa526-id_30.csv.gz',
    'data/zemu_control-69aa526-id_50.csv.gz',
    'data/zemu_control-69aa526-WildTypeComplex_01.csv.gz',
    'data/zemu_control-69aa526-WildTypeComplex_30.csv.gz',
    'data/zemu_control-69aa526-WildTypeComplex_50.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-id_01.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-id_30.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-id_50.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-WildTypeComplex_01.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-WildTypeComplex_30.csv.gz',
    'data/ddg_monomer_16_003-zemu-2-WildTypeComplex_50.csv.gz',
    'data/zemu-brub_1.6-nt10000-id_01.csv.gz',
    'data/zemu-brub_1.6-nt10000-id_30.csv.gz',
    'data/zemu-brub_1.6-nt10000-id_50.csv.gz',
    'data/zemu-brub_1.6-nt10000-WildTypeComplex_01.csv.gz',
    'data/zemu-brub_1.6-nt10000-WildTypeComplex_30.csv.gz',
    'data/zemu-brub_1.6-nt10000-WildTypeComplex_50.csv.gz',
    'data/zemu-values-id_01.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-id_01.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-id_30.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-id_50.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-WildTypeComplex_01.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-WildTypeComplex_30.csv.gz',
    'data/zemu_1.2-60000_rscript_simplified-t14-WildTypeComplex_50.csv.gz'
]

output_dir = 'output'
latex_output_dir = os.path.join( output_dir, 'latex' )
output_fig_path = os.path.join( output_dir, 'figures_and_tables' )
main_template_latex_file = os.path.join( 'latex_templates', 'figs_and_tables.tex' )
if not os.path.isdir( output_fig_path ):
    os.makedirs( output_fig_path )
if not os.path.isdir( latex_output_dir ):
    os.makedirs( latex_output_dir )

# Import subsets
with open('subsets.json') as f:
    subsets = json.load(f)

# Import here as they can be slow, and are unneeded if plots aren't going to be made
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
current_palette = [
    (0.29803921568627451, 0.44705882352941179, 0.69019607843137254),
    (0.33333333333333331, 0.6588235294117647, 0.40784313725490196),
    (0.7686274509803922, 0.30588235294117649, 0.32156862745098042),
    (0.50588235294117645, 0.44705882352941179, 0.69803921568627447),
    (0.80000000000000004, 0.72549019607843135, 0.45490196078431372),
    (0.39215686274509803, 0.70980392156862748, 0.80392156862745101),
]

for csv_path in csv_paths:
    if not os.path.isfile(csv_path):
        print csv_path
        assert( os.path.isfile(csv_path) )
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

display_mut_types = [
    'complete', 's2l', 'sing_ala', 'mult_mut', 'mult_none_ala',
]

mut_types = {
    'complete' : 'Complete dataset',
    'sing_mut' : 'Single mutation',
    'mult_mut' : 'Multiple mutations',
    'mult_all_ala' : 'Multiple mutations, all to alanine',
    'mult_none_ala' : 'Multiple mutations, none to alanine',
    's2l' : 'Small-to-large mutation(s)',
    'l2s' : 'Large-to-small',
    'ala' : 'Mutation(s) to alanine',
    'sing_ala' : 'Single mutation to alanine',
    'res_gte25' : 'Res. $>=$ 2.5 Ang.',
    'res_lte15' : 'Res. $<=$ 1.5 Ang.',
    'res_gt15_lt25' : '1.5 Ang. $<$ Res. $<$ 2.5 Ang.',
    'some_s2l' : 'Some small-to-large',
    'some_l2s' : 'Some large-to-small',
    'antibodies' : 'Antibodies',
    'stabilizing' : 'Stabilizing',
    'neutral' : 'Neutral',
    'positive' : 'Positive',
    's2l_stabilizing' : 'Small-to-large and stabilizing',
    'SS-all-None' : 'Loop',
    'SS-all-H' : 'Alpha Helix',
    'SS-all-E' : 'Strand',
    'SS-all-T' : 'Turn',
    'SS-all-B' : 'Beta Bridge',
    'SS-all-G' : 'Helix-3',
    # 'SS-all-I' : 'Helix-5',
    'SS-all-S' : 'Bend',
}

for cut_off in [0.05, 0.1, 0.25]:
    for stat_type in ['mean']:
        for comparison_type, comparison_type_name in [('lt', '<'), ('gte', '>=')]:
            mut_types[ 'buried_%s%.2f_%s' % (comparison_type, cut_off, stat_type) ] = '%s solvent exposure %s %.2f' % (stat_type.capitalize(), comparison_type_name, cut_off)

run_names = {
    'zemu_1.2-60000_rscript_validated-t14' : 'flex ddG',
    'zemu_1.2-60000_rscript_simplified-t14' : 'flex ddG',
    'zemu_1.2-60000_rscript_validated-ref' : 'flex ddG (REF energy)',
    'zemu-brub_1.6-nt10000' : 'flex ddG (1.6 kT)',
    'ddg_monomer_16_003-zemu-2' : 'ddG monomer',
    'zemu_control' : 'no backrub control',
    'zemu-values' : 'ZEMu paper',
    'zemu_control-69aa526-noglypivot' : 'no backrub control',
    'zemu_control-69aa526' : 'no backrub control',
    'tal_GAM' : 'GAM (Talaris)',
    'ref_GAM' : 'GAM (REF)',
    'control_GAM' : 'GAM (No backrub control, Talaris)',
}

run_colors = {
    run_name : current_palette[i]
    for i, run_name in enumerate([
        'zemu_1.2-60000_rscript_simplified-t14',
        'zemu_control-69aa526',
        'zemu-values',
        'ddg_monomer_16_003-zemu-2',
    ])
}

sorting_type_descriptions = {
    'id' : 'Structures are not sorted, and are randomly added to the ensemble. ',
    'WildTypeComplex' : 'Structures are sorted by their minimized wild-type complex energy. ',
}

sorting_type_descriptions_short = {
    'id' : 'unsorted structures',
    'WildTypeComplex' : 'structures sorted by wild-type complex energy',
}

def make_fig_components_text(fc, run_name, score_method, structure_order):
    matching_uses = set()
    score_method = int(score_method)
    for use_name, uses_list in fc.iteritems():
        for other_run_name, other_structure_order, other_score_method in uses_list:
            if other_run_name == run_name and other_structure_order == structure_order and other_score_method == score_method:
                matching_uses.add( use_name )
    return '; '.join(sorted(list(matching_uses)))

def figure_components():
    # Figure to run map
    def append_all_structs_helper( runs_steps, structure_orders ):
        l = []
        for run_name, step_num in runs_steps:
            for structure_order in structure_orders:
                for i in xrange(1, 51):
                    if i == 50 and structure_order == 'WildTypeComplex':
                        l.append( (run_name, '%s_%02d' % ('id', i), step_num ) )
                    else:
                        l.append( (run_name, '%s_%02d' % (structure_order, i), step_num ) )
        return l

    fig_4 = [ ('zemu_control-69aa526', 'id_50', 8) ]
    for step in np.arange(0,50001,2500):
        fig_4.append( ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', step) )
    total_ddg_score_figure_components = {
        'Figure 2' : [
            ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', 35000),
            ('zemu_control-69aa526', 'id_50', 8),
        ],
        'Figure 3' : append_all_structs_helper([('zemu_1.2-60000_rscript_simplified-t14', 35000), ('zemu_control-69aa526', 8)], ['WildTypeComplex']),
        'Figure 4' : fig_4,
        'Figure 5' : [
            ('tal_GAM', 'id_50', 12),
            ('control_GAM', 'id_50', 12),
        ],
        'Table 2' : [
            ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', 35000),
            ('zemu_control-69aa526', 'id_50', 8),
            ('zemu-values', 'id_01', 11),
            ('ddg_monomer_16_003-zemu-2', 'id_50', 8),
        ],
        'Table S2' : [
            ('zemu-brub_1.6-nt10000', 'id_50', 10000),
            ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', 10000),
        ],
        'Table S3' : [
            ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', 35000),
            ('zemu_control-69aa526', 'id_50', 8),
            ('zemu-values', 'id_01', 11),
            ('ddg_monomer_16_003-zemu-2', 'id_50', 8),
        ],
        'Table S4' : append_all_structs_helper([('zemu_1.2-60000_rscript_simplified-t14', 35000), ('zemu_control-69aa526', 8)], ['WildTypeComplex']),
        'Table S5' : append_all_structs_helper([('ddg_monomer_16_003-zemu-2', 8), ('zemu_control-69aa526', 8)], ['WildTypeComplex']),
        'Table S6' : append_all_structs_helper([('zemu_1.2-60000_rscript_simplified-t14', 35000), ('zemu_control-69aa526', 8)], ['id']),
        'Table S7' : fig_4,
        'Table S8' : [
            ('zemu_1.2-60000_rscript_simplified-t14', 'id_50', 35000),
            ('zemu_1.2-60000_rscript_validated-ref', 'id_50', 35000),
        ],
        'Table S9' : [
            ('tal_GAM', 'id_50', 12),
            ('control_GAM', 'id_50', 12),
            ('ref_GAM', 'id_50', 12),
        ],
        'Figure S2' : append_all_structs_helper([('ddg_monomer_16_003-zemu-2', 8), ('zemu_control-69aa526', 8)], ['WildTypeComplex']),
        'Figure S3' : append_all_structs_helper([('zemu_1.2-60000_rscript_simplified-t14', 35000), ('zemu_control-69aa526', 8)], ['id']),
        'Figure S5' : [
            ('zemu_1.2-60000_rscript_validated-ref', 'id_50', 35000),
            ('ref_GAM', 'id_50', 12),
        ],
        'Figure S6' : [
            ('tal_GAM', 'id_50', 12),
        ],
    }

    partial_ddg_score_figure_components = copy.deepcopy( total_ddg_score_figure_components )
    fig_s1 = []
    for step in np.arange(0,50001,2500):
        fig_s1.extend( append_all_structs_helper( [('zemu_1.2-60000_rscript_simplified-t14', step)], ['WildTypeComplex']) )
    partial_ddg_score_figure_components['Figure S1'] = fig_s1 # Not used now, but could be in the future

    return total_ddg_score_figure_components

cached_loaded_df_initialized = False
cached_loaded_df = None
cached_df_path = '/tmp/read_all_scores.hdf'
def load_df():
    global cached_loaded_df_initialized
    global cached_loaded_df
    if cached_loaded_df_initialized:
        return cached_loaded_df.copy()
    elif os.path.isfile(cached_df_path):
        print 'Loading cached HDF'
        df = pd.read_hdf( cached_df_path )
        cached_loaded_df = df.copy()
        cached_loaded_df_initialized = True
        print 'Done loading HDF\n'
        return df

    print 'Performing initial .csv load'

    # Load CSVs
    df = pd.read_csv( csv_paths[0] )
    for csv_path in csv_paths[1:]:
        df = df.append( pd.read_csv( csv_path ) )

    # Load GAM
    for run_name, csv_path in [
            ('tal_GAM', os.path.join('gam', 'tal_GAM_terms.csv')),
            ('ref_GAM', os.path.join('gam', 'ref_GAM_terms.csv')),
            ('control_GAM', os.path.join('gam', 'control_GAM_terms.csv')),
    ]:
        gam_df = pd.read_csv( csv_path )
        assert( len(gam_df) == 1240 )
        gam_df['DataSetID'] = range(1, 1241)
        gam_df = gam_df.assign( PredictionRunName = run_name )
        gam_df = gam_df.assign( ScoreMethodID = 12 )
        gam_df = gam_df.assign( PredictionID = None )
        gam_df = gam_df.assign( StructureOrder = 'id_50' )
        gam_df['total'] = gam_df[ [x for x in gam_df.columns if x.endswith('_GAM')] ].sum( axis = 1 )
        gam_df = gam_df[ ['PredictionRunName', 'DataSetID', 'ScoreMethodID', 'total', 'exp_data', 'StructureOrder'] ]
        gam_df.columns = [ 'PredictionRunName', 'DataSetID', 'ScoreMethodID', 'total', 'ExperimentalDDG', 'StructureOrder' ]
        df = df.append( gam_df )

    # Finish up
    df = add_score_categories( df )
    df = df.drop_duplicates( ['PredictionRunName', 'DataSetID', 'PredictionID', 'ScoreMethodID', 'MutType', 'total', 'ExperimentalDDG', 'StructureOrder'] )
    cached_loaded_df_initialized = True
    cached_loaded_df = df.copy()
    df.to_hdf( cached_df_path, 'read_all_scores' )
    print 'Done loading csvs\n'
    return df

def add_score_categories(df, mut_type_subsets = None):
    if mut_type_subsets == None or 'complete' in mut_type_subsets or 'MutType' not in df.columns:
        df = df.assign( MutType = 'complete' )

    if mut_type_subsets == None or 'stabilizing' in mut_type_subsets:
        stabilizing = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] <= -1.0) ].copy()
        stabilizing.loc[:,'MutType'] = 'stabilizing'
        df = df.append( stabilizing )

        # Add s2l_stabilizing
        s2l_stabilizing = df.loc[ (df['MutType'] == 's2l') & (df['ExperimentalDDG'] <= -1.0) ].copy()
        s2l_stabilizing.loc[:,'MutType'] = 's2l_stabilizing'
        df = df.append( s2l_stabilizing )

    if mut_type_subsets == None or 'neutral' in mut_type_subsets:
        neutral = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] > -1.0) & (df['ExperimentalDDG'] < 1.0) ].copy()
        neutral.loc[:,'MutType'] = 'neutral'
        df = df.append( neutral )

    if mut_type_subsets == None or 'positive' in mut_type_subsets:
        positive = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] >= 1.0) ].copy()
        positive.loc[:,'MutType'] = 'positive'
        df = df.append( positive )

    for subset_name, subset_keys in subsets.iteritems():
        if mut_type_subsets == None or subset_name in mut_type_subsets:
            subset_df = df.loc[ (df['MutType'] == 'complete') & (df['DataSetID'].isin(subset_keys)) ].copy()
            assert( len(subset_df) > 0 )
            subset_df.loc[:,'MutType'] = subset_name
            df = df.append( subset_df )

    return df

def save_latex( latex_template_file, sub_dict, out_tex_name = None, long_table = False, wrap_command = None ):
    if 'fig-path' in sub_dict:
        shutil.copy( sub_dict['fig-path'], latex_output_dir )
        sub_dict['fig-path'] = os.path.basename( sub_dict['fig-path'] )

    if out_tex_name == None:
        out_tex_name = os.path.basename( latex_template_file )
    if not out_tex_name.endswith( '.tex' ):
        out_tex_name += '.tex'

    with open(latex_template_file, 'r') as f:
        latex_template = f.read()
    for key in sub_dict:
        latex_key = '%%%%%s%%%%' % key
        latex_template = latex_template.replace( latex_key, sub_dict[key] )

    if long_table:
        latex_lines = latex_template.split('\n')
        new_lines = []
        for line in latex_lines:
            if '\\begin{table}' in line or '\\end{tabular}' in line:
                continue
            elif 'begin{tabular}' in line:
                col_sorting = line.strip()[ len('\\begin{tabular}') : ]
                new_lines.append( '\\begin{longtable}' + col_sorting )
            elif 'end{table}' in line:
                new_lines.append( '\\end{longtable}' )
            else:
                new_lines.append(line)
        latex_template = '\n'.join( new_lines )

    if wrap_command != None:
        latex_template = '{' + wrap_command + '\n' + latex_template + '\n}'

    with open( os.path.join( latex_output_dir, out_tex_name ), 'w') as f:
        f.write( latex_template )

def isfloat(x):
    try:
        x = float(x)
    except ValueError:
        return False

    return True

def compile_latex():
    shutil.copy( main_template_latex_file, latex_output_dir )

    for i in xrange( 0, 3 ):
        output = subprocess.check_output( ['xelatex', os.path.basename( main_template_latex_file )[:-4] ], cwd = latex_output_dir )

def make_results_df( generate_plots = False, print_statistics = False, use_cached_if_available = True ):
    results_csv_path = os.path.join(output_dir, 'results.csv')
    if use_cached_if_available and os.path.isfile( results_csv_path ):
        results_df = pd.read_csv( results_csv_path )
        return results_df

    df = load_df()

    # Add 'complete' to front of list so scaling is calculated first
    mut_types = sorted( df['MutType'].drop_duplicates().values )
    mut_types.insert(0, mut_types.pop( mut_types.index('complete') ) )

    structure_orders = sorted( df['StructureOrder'].drop_duplicates().values )
    maes_unscaled = []
    fcs_unscaled = []
    fcs_unscaled_correct = []
    fcs_unscaled_considered = []
    fcs = []
    fcs_correct = []
    fcs_considered = []
    fc_zeros = []
    fc_zeros_correct = []
    fc_zeros_considered = []
    maes = []
    rs = []
    slopes = []
    df_mut_types = []
    df_structure_orders = []
    prediction_runs = []
    steps = []
    ns = []
    scaling_factors = {}
    scaling_slopes = []
    scaling_intercepts = []

    for prediction_run in sorted( df['PredictionRunName'].drop_duplicates().values ):
        for structure_order in structure_orders:
            for mut_type in mut_types:
                outer_sub_df = df.loc[ (df['MutType'] == mut_type) & (df['StructureOrder'] == structure_order) & (df['PredictionRunName'] == prediction_run) ]
                num_steps = float( len(outer_sub_df['ScoreMethodID'].drop_duplicates().values) )
                if generate_plots:
                    outer_fig = plt.figure(figsize=(8.5, 8.5), dpi=600)
                    outer_ax = outer_fig.add_subplot(1, 1, 1)

                for step_i, step in enumerate( sorted( outer_sub_df['ScoreMethodID'].drop_duplicates().values ) ):
                    sub_df = outer_sub_df.loc[ outer_sub_df['ScoreMethodID'] == step ]

                    sub_output_dir = os.path.join(output_dir, prediction_run.replace('.', ''))
                    sub_output_dir = os.path.join(sub_output_dir, mut_type)
                    sub_output_dir = os.path.join(sub_output_dir, structure_order)
                    outer_fig_path = os.path.join(sub_output_dir, 'all_steps.pdf')
                    sub_output_dir = os.path.join(sub_output_dir, str(step))

                    if not os.path.isdir(sub_output_dir) and generate_plots:
                        os.makedirs(sub_output_dir)

                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress( sub_df['total'], sub_df['ExperimentalDDG'] )
                    rs.append(r_value)
                    slopes.append(slope)
                    if mut_type == 'complete':
                        assert( (prediction_run, structure_order, step) not in scaling_factors )
                        scaling_factors[(prediction_run, structure_order, step)] = (slope, intercept)

                    sub_df.loc[:,'scaled_total'] = sub_df.loc[:,'total'] * scaling_factors[(prediction_run, structure_order, step)][0] + scaling_factors[(prediction_run, structure_order, step)][1]
                    scaling_slopes.append( scaling_factors[(prediction_run, structure_order, step)][0] )
                    scaling_intercepts.append( scaling_factors[(prediction_run, structure_order, step)][1] )

                    fc_unscaled, fc_unscaled_correct, fc_unscaled_considered = fraction_correct(
                        sub_df['total'].values,
                        sub_df['ExperimentalDDG'].values,
                    )
                    fcs_unscaled.append(fc_unscaled)
                    fcs_unscaled_correct.append(fc_unscaled_correct)
                    fcs_unscaled_considered.append(fc_unscaled_considered)
                    fc, fc_correct, fc_considered = fraction_correct(
                        sub_df['scaled_total'].values,
                        sub_df['ExperimentalDDG'].values,
                    )
                    fcs.append(fc)
                    fcs_correct.append(fc_correct)
                    fcs_considered.append(fc_considered)
                    fc_zero, fc_zero_correct, fc_zero_considered = fraction_correct(
                        sub_df['total'].values,
                        sub_df['ExperimentalDDG'].values,
                        x_cutoff = 0.00000001,
                        y_cutoff = 0.00000001,
                    )
                    fc_zeros.append(fc_zero)
                    fc_zeros_correct.append(fc_zero_correct)
                    fc_zeros_considered.append(fc_zero_considered)

                    mae_unscaled = calc_mae( sub_df['total'], sub_df['ExperimentalDDG'] )
                    maes_unscaled.append(mae_unscaled)
                    mae = calc_mae( sub_df['scaled_total'], sub_df['ExperimentalDDG'] )
                    maes.append(mae)

                    df_mut_types.append(mut_type)
                    df_structure_orders.append(structure_order)
                    prediction_runs.append(prediction_run)
                    steps.append(step)
                    ns.append(len(sub_df['ExperimentalDDG']))

                    if print_statistics:
                        print mut_type, structure_order, prediction_run
                        print 'R:', '%.3f' % r_value
                        print 'slope:', '%.3f' % slope
                        print 'FC:', '%.3f' % fc
                        print 'MAE:', '%.3f' % mae
                        print
                    if generate_plots:
                        fig = plt.figure(figsize=(8.5, 8.5), dpi=300)
                        ax = fig.add_subplot(1, 1, 1)
                        sns.regplot(y="total", x="ExperimentalDDG", data = sub_df, ax=ax, scatter_kws = { 's' : 4.5 } )

                        sns.regplot(
                            y="total", x="ExperimentalDDG", data = sub_df, ax=outer_ax,
                            color = ( float(num_steps - step_i) / num_steps, 0, ( float(step_i) / num_steps ) ),
                            scatter_kws = { 's' : 2 },
                        )

                        ax.set_xlabel( 'Experimental ddG' )
                        ax.set_ylabel( 'ddG prediction' )

                        fig_path = os.path.join(sub_output_dir, 'total_vs_experimental.pdf')
                        fig.savefig( fig_path )
                        plt.close(fig)

                if generate_plots:
                    outer_fig.savefig(outer_fig_path)
                    plt.close(outer_fig)

    results_df = pd.DataFrame( {
        'PredictionRun' : prediction_runs,
        'MutTypes' : df_mut_types,
        'StructureOrder' : df_structure_orders,
        'FractionCorrect' : fcs,
        'FractionCorrect-count_correct' : fcs_correct,
        'FractionCorrect-count_considered' : fcs_considered,
        'FractionCorrect_unscaled' : fcs_unscaled,
        'FractionCorrect-count_correct-unscaled' : fcs_unscaled_correct,
        'FractionCorrect-count_considered-unscaled' : fcs_unscaled_considered,
        'FractionCorrect-zero_boundary-unscaled' : fc_zeros,
        'FractionCorrect-count_correct-zero_boundary-unscaled' : fc_zeros_correct,
        'FractionCorrect-count_considered-zero_boundary-unscaled' : fc_zeros_considered,
        'MAE' : maes,
        'MAE_unscaled' : maes_unscaled,
        'Slope' : slopes,
        'ScalingSlopes' : scaling_slopes,
        'ScalingIntercepts' : scaling_intercepts,
        'Step' : steps,
        'R' : rs,
        'N' : ns,
    } )

    sort_cols = [
        'FractionCorrect_unscaled',
        'R',
        'MAE_unscaled',
    ]
    ascendings = [False, False, True]
    for sort_col, asc in zip(sort_cols, ascendings):
        results_df.sort_values( sort_col, inplace = True, ascending = asc )
        # print results_df[ ['PredictionRun', 'MutTypes', 'StructureOrder', 'Step', sort_col] ].head(n=20)
    new_column_order = [
        'PredictionRun', 'N', 'MutTypes', 'Step', 'StructureOrder',
        'R', 'MAE', 'MAE_unscaled',
        'Slope', 'ScalingIntercepts', 'ScalingSlopes',
    ]
    for col in sorted(results_df.columns):
        if col not in new_column_order:
            new_column_order.append( col )
    assert( set(new_column_order) == set(results_df.columns) )
    results_df = results_df[new_column_order]
    results_df.sort_values('R', ascending = False).to_csv( results_csv_path )
    if print_statistics:
        print results_csv_path
        print

    return results_df

def figure_mult_non_ala():
    exp_run_name = 'zemu_1.2-60000_rscript_simplified-t14'
    control_run_name = 'zemu_control-69aa526'
    zemu_run_name = 'zemu-values'
    force_backrub_step = 35000
    point_size = 6.0
    alpha = 1.0
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.8 }

    df = load_df()
    exp_colname = 'Experimental ddG'
    pred_colname = 'Rosetta Score'
    display_subset = 'mult_none_ala'
    df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
    df = df.rename( columns = {'total' : pred_colname} )

    sns.set_style("whitegrid")
    fig = plt.figure(
        figsize=(6, 8.5), dpi=600
    )

    df_a = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == display_subset) & (df['ScoreMethodID'] == force_backrub_step) & (df['StructureOrder'] == 'id_50') ]
    ax1 = fig.add_subplot( 2, 1, 1 )

    df_b = df.loc[ (df['PredictionRunName'] == zemu_run_name) & (df['MutType'] == display_subset) ]
    ax2 = fig.add_subplot( 2, 1, 2 )

    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_a, ax = ax1,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[exp_run_name],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_b, ax = ax2,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[exp_run_name],
    )

    # ax1.set_xlabel('')
    ax2.set_ylabel('ZEMu score')
    fig.suptitle(mut_types[display_subset])

    out_path = os.path.join( output_fig_path, 'fig-scatter_mult-non-ala.pdf' )
    fig.savefig( out_path )

def figure_scatter( force_backrub_step = 35000 ):
    exp_run_name = 'zemu_1.2-60000_rscript_simplified-t14'
    control_run_name = 'zemu_control-69aa526'
    point_size = 4.8
    alpha = 0.55
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.8 }

    df = load_df()
    exp_colname = 'Experimental ddG'
    pred_colname = 'Rosetta Score'
    top_subset = 'complete'
    bottom_subset = 's2l'
    df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
    df = df.rename( columns = {'total' : pred_colname} )

    sns.set_style("whitegrid")
    fig = plt.figure(
        figsize=(8.5, 8.5), dpi=600
    )

    complete_corrs = df.loc[ (df['MutType'] == top_subset) & (df['PredictionRunName'] == exp_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_a = force_backrub_step
    df_a = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == top_subset) & (df['ScoreMethodID'] == best_step_a) ]
    ax1 = fig.add_subplot( 2, 2, 1 )

    complete_corrs = df.loc[ (df['MutType'] == top_subset) & (df['PredictionRunName'] == control_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_b = complete_corrs.index[0][0]
    df_b = df.loc[ (df['PredictionRunName'] == control_run_name) & (df['MutType'] == top_subset) & (df['ScoreMethodID'] == best_step_b) ]
    ax2 = fig.add_subplot( 2, 2, 2 )

    complete_corrs = df.loc[ (df['MutType'] == bottom_subset) & (df['PredictionRunName'] == exp_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_c = force_backrub_step
    df_c = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == bottom_subset) & (df['ScoreMethodID'] == best_step_c) ]
    ax3 = fig.add_subplot( 2, 2, 3 )

    complete_corrs = df.loc[ (df['MutType'] == bottom_subset) & (df['PredictionRunName'] == control_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_d = complete_corrs.index[0][0]
    df_d = df.loc[ (df['PredictionRunName'] == control_run_name) & (df['MutType'] == bottom_subset) & (df['ScoreMethodID'] == best_step_d) ]
    ax4 = fig.add_subplot( 2, 2, 4 )

    xmin = min( df_a[pred_colname].min(), df_b[pred_colname].min(), df_c[pred_colname].min(), df_d[pred_colname].min() )
    xmax = max( df_a[pred_colname].max(), df_b[pred_colname].max(), df_c[pred_colname].max(), df_d[pred_colname].max() )
    ymin = min( df_a[exp_colname].min(), df_b[exp_colname].min(), df_c[exp_colname].min(), df_d[exp_colname].min() )
    ymax = max( df_a[exp_colname].max(), df_b[exp_colname].max(), df_c[exp_colname].max(), df_d[exp_colname].max() )
    ax1.set_xlim( xmin, xmax )
    ax1.set_ylim( ymin, ymax )
    ax2.set_xlim( xmin, xmax )
    ax2.set_ylim( ymin, ymax )
    ax3.set_xlim( xmin, xmax )
    ax3.set_ylim( ymin, ymax )
    ax4.set_xlim( xmin, xmax )
    ax4.set_ylim( ymin, ymax )

    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_a, ax = ax1,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[exp_run_name],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_b, ax = ax2,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[control_run_name],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_c, ax = ax3,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[exp_run_name],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_d, ax = ax4,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = run_colors[control_run_name],
    )

    # ax2.set_xticklabels([])
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    # ax2.set_yticklabels([])
    ax4.set_ylabel('')
    ax3.set_xlabel('Experimental $\Delta\Delta G$')
    ax4.set_xlabel('Experimental $\Delta\Delta G$')

    ax1.set_title( '(a) Flex ddG - %s' % (mut_types[top_subset]) )
    ax2.set_title( '(b) Control - %s' % (mut_types[top_subset] ) )
    ax3.set_title( '(c) Flex ddG - %s' % (mut_types[bottom_subset]) )
    ax4.set_title( '(d) Control - %s' % (mut_types[bottom_subset]) )

    # Assert that these lengths are equal since N is displayed for top and bottom together
    assert( len(df_a) == len(df_b) )
    assert( len(df_c) == len(df_d) )

    out_path = os.path.join( output_fig_path, 'fig-scatter.pdf' )
    bottom_subset_label = mut_types[bottom_subset].capitalize()
    if bottom_subset_label.strip() == 'Small-to-large mutation(s)':
        bottom_subset_label = 'Small-to-large mutation(s) subset'
    sub_dict = {
        'exp-method-name' : run_names[exp_run_name],
        'numsteps-a' : str( best_step_a ),
        'numsteps-c' : str( best_step_c ),
        'top-subset' : mut_types[top_subset].capitalize(),
        'top-n' : str( len(df_a) ),
        'control-method-name' : run_names[control_run_name],
        'bottom-subset' : mut_types[bottom_subset].capitalize(),
        'bottom-n' : str( len(df_c) ),
        'fig-path' : out_path,
    }

    fig.savefig( out_path )
    save_latex( 'latex_templates/figure-scatter.tex', sub_dict )
    print out_path

def table_composition():
    # Dataset composition
    df = load_df()
    df = df.drop_duplicates( ['PredictionRunName', 'DataSetID', 'PredictionID', 'ScoreMethodID', 'MutType', 'total', 'ExperimentalDDG', 'StructureOrder'] )
    control_df = df.loc[ (df['PredictionRunName'] == 'zemu_control-69aa526') & (df['ScoreMethodID'] == 8 ) ]

    ns = []
    mut_type_names = []

    for mut_type in display_mut_types:
        ns.append( len( control_df.loc[ control_df['MutType'] == mut_type ] ) )
        mut_type_names.append( mut_types[mut_type] )

    table_df = pd.DataFrame( {
        'Name' : mut_type_names,
        'n' : ns,
    } )

    table_df = table_df[ ['n', 'Name'] ] # Order columns correctly
    table_df.sort_values( 'n', inplace = True, ascending = False )

    # print 'Table 1:'
    # print table_df.head( n = 20 )
    save_latex( 'latex_templates/table-composition.tex', { 'table-comp' : table_df.to_latex( index=False ) } )
    table_df.to_csv( os.path.join(output_fig_path, 'table_comp.csv') )

def table_versions():
    save_latex( 'latex_templates/table-versions.tex', run_names )

def steps_vs_corr( output_figure_name, mut_type_subsets, control_run = 'zemu_control-69aa526', force_control_in_axes = True, max_backrub_steps = 50000, structure_order = 'id_50' ):
    exp_run_name = 'zemu_1.2-60000_rscript_simplified-t14'

    df = load_df()
    exp_colname = 'Experimental ddG'
    pred_colname = 'Rosetta Score'
    df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
    df = df.rename( columns = {'total' : pred_colname} )

    # Filter to desired StructureOrder
    if structure_order != None:
        df = df.loc[ df['StructureOrder'] == structure_order ]

    if max_backrub_steps != None:
        df = df.loc[ df['ScoreMethodID'] <= max_backrub_steps ]

    sns.set_style("white")
    fig = plt.figure(
        figsize=(10.0, 8.5), dpi=600
    )
    fig.subplots_adjust( wspace = 0.6, hspace = 0.3)
    # fig.suptitle('$\Delta\Delta G$ prediction performance vs. number of backrub sampling steps')

    r_axes = []
    r_min = float('inf')
    r_max = float('-inf')

    mae_axes = []
    mae_min = float('inf')
    mae_max = float('-inf')

    ns = []
    legend_lines = []
    legend_labels = []
    data_table = []
    for ax_i, mut_type_subset in enumerate( mut_type_subsets ):
        ax = fig.add_subplot( 2, 2, ax_i + 1 )
        ax.set_title( '(%s) - %s' % (string.ascii_lowercase[ax_i], mut_types[mut_type_subset]) )
        ax.set_ylabel("Pearson's R")
        ax.set_xlabel("Backrub Step")
        ns.append( len(df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) & (df['ScoreMethodID'] == df['ScoreMethodID'].drop_duplicates().values[0]) ]) )
        rs = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
        maes = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )
        assert( len(rs) == len(maes) )
        for name, group in df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID'):
            assert( len(group) == len(subsets[mut_type_subset]) )

        exp_rs_plot, = ax.plot(
            rs['ScoreMethodID'], rs['Experimental ddG'],
            'o',
            linestyle = 'None',
            color = run_colors[exp_run_name],
        )
        r_min = min( r_min, min(rs['Experimental ddG']) )
        r_max = max( r_max, max(rs['Experimental ddG']) )

        ax2 = ax.twinx()
        ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
        ax2.set_ylabel("MAE")
        exp_mae_plot, = ax2.plot(
            maes.index, maes.values,
            'P',
            linestyle = 'None',
            color = ( 1.0 + np.array(run_colors[exp_run_name]) ) / 2.0,
        )
        mae_min = min( mae_min, min(maes.values) )
        mae_max = max( mae_max, max(maes.values) )

        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))

        # Get control performance values
        control_rs = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
        assert( len(control_rs) == 1 )
        control_r = control_rs['Experimental ddG'][0]
        control_rs_plot, = ax.plot(
            [0], [control_r],
            marker = 'o',
            linestyle = 'None',
            color = run_colors[control_run],
        )

        control_maes = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )
        assert( len(control_maes.values) == 1 )
        control_mae = control_maes.values[0]
        control_mae_plot, = ax2.plot(
            [0], [control_mae],
            marker = 'P',
            linestyle = 'None',
            color = ( 1.0 + np.array(run_colors[control_run]) ) / 2.0,
        )

        np.testing.assert_array_equal( rs['ScoreMethodID'], maes.index )
        data_table.append( (run_names[control_run], mut_types[mut_type_subset], 0, control_r, control_mae) )
        for step_num, r, mae in zip(rs['ScoreMethodID'], rs['Experimental ddG'], maes.values):
            if step_num in [2500, 35000] or step_num % 10000 == 0:
                data_table.append( (run_names[exp_run_name], mut_types[mut_type_subset], step_num, r, mae) )


        if force_control_in_axes:
            r_min = min( r_min, control_r )
            r_max = max( r_max, control_r )
            mae_min = min( mae_min, control_mae )
            mae_max = max( mae_max, control_mae )

        if ax_i == 0:
            legend_lines.extend( [exp_rs_plot, control_rs_plot, exp_mae_plot, control_mae_plot] )
            legend_labels.extend( [
                'R - ' + run_names[exp_run_name],
                'R - ' + run_names[control_run],
                'MAE - ' + run_names[exp_run_name],
                'MAE - ' + run_names[control_run],
            ] )

        r_axes.append( ax )
        mae_axes.append( ax2 )

    r_range = r_max - r_min
    r_min -= r_range * 0.1
    r_max += r_range * 0.1
    mae_range = mae_max - mae_min
    mae_min -= mae_range * 0.1
    mae_max += mae_range * 0.1

    for r_ax, mae_ax in zip(r_axes, mae_axes):
        r_ax.set_ylim( [r_min, r_max] )
        mae_ax.set_ylim( [mae_min, mae_max] )
        r_ax.set_zorder( mae_ax.get_zorder() + 10 )
        r_ax.patch.set_visible(False)

    out_path = os.path.join( output_fig_path, '%s.pdf' % output_figure_name )
    underlying_name = '%s-underlying-data' % output_figure_name
    color_description = ''
    if exp_run_name == 'zemu_1.2-60000_rscript_simplified-t14':
        color_description += '\nPredictions generated with the Flex ddG protocol are shown in blue.'
    elif exp_run_name == 'ddg_monomer_16_003-zemu-2':
        color_description += '\nPredictions generated with the ddg\_monomer method are shown in purple.'
    color_description += '\nPredictions generated with the no backrub control protocol are shown in green.'
    sub_dict = {
        'panel-a' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[0] ].capitalize(),  ns[0] ),
        'panel-b' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[1] ].capitalize(),  ns[1] ),
        'panel-c' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[2] ].capitalize(),  ns[2] ),
        'panel-d' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[3] ].capitalize(),  ns[3] ),
        'fig-label' : output_figure_name,
        'fig-path' : out_path,
        'underlying-label' : 'tab:%s' % underlying_name,
        'color-description' : color_description,
    }

    leg = plt.figlegend(
        legend_lines, legend_labels, loc = 'lower center', labelspacing=0.0,
        markerscale = 1.3,
        fontsize = 13,
        ncol = 2,
    )

    fig.savefig( out_path )
    save_latex( 'latex_templates/steps-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )

    fig.savefig( out_path )
    print out_path

    # Save underlying data
    data_table = pd.DataFrame.from_records(
        data_table,
        columns = ['Run', 'Subset', 'Backrub Step', 'R', 'MAE'],
    )
    latex_lines = data_table.to_latex( float_format = '%.2f', index = False ).split('\n')
    latex_lines = [r'\begin{table}'] + latex_lines + ['\caption[]{Selection of key data shown in \cref{fig:%s}}' % output_figure_name, '\label{tab:%s}' % underlying_name, r'\end{table}', '']
    with open( 'output/latex/%s.tex' % underlying_name, 'w' ) as f:
        f.write( '\n'.join(latex_lines) )

def figure_structs_vs_corr(
        exp_run_name = 'zemu_1.2-60000_rscript_simplified-t14', force_backrub_step = 35000, max_backrub_step = 50000,
        mut_type_subsets = ['complete', 's2l', 'mult_none_ala', 'sing_ala'],
        extra_fig_name = None,
):
    sorting_types = ['WildTypeComplex', 'id']
    base_path = 'data/by_struct/%s-%s_%02d.csv.gz'
    control_base_path = 'data/by_struct/%s-%s_%02d.csv.gz'
    number_of_structures = 50
    control_run = 'zemu_control-69aa526'

    for sorting_type in sorting_types:
        structure_orders = []
        df = None
        for struct_id in xrange(1, number_of_structures + 1):
            csv_path = base_path % (exp_run_name, sorting_type, struct_id)
            assert( os.path.isfile( csv_path ) )

            if struct_id == 1:
                df = pd.read_csv( csv_path )
            else:
                df = df.append( pd.read_csv( csv_path ) )

            control_csv_path = control_base_path % (control_run, sorting_type, struct_id)
            df = df.append( pd.read_csv( control_csv_path ) )

            structure_orders.append( '%s_%02d' % (sorting_type, struct_id) )

        df['StructureOrder'] = df['StructureOrder'].apply(
            lambda s: int( s.split('_')[1] )
        )

        df = add_score_categories(df, mut_type_subsets = mut_type_subsets)

        point_size = 4.5
        alpha = 0.6
        scatter_kws = { 's' : point_size, 'alpha' : alpha }
        line_kws = { 'linewidth' : 0.9 }

        exp_colname = 'Experimental ddG'
        pred_colname = 'Rosetta Score'
        df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
        df = df.rename( columns = {'total' : pred_colname} )

        if max_backrub_step != None:
            df = df.loc[ df['ScoreMethodID'] <= max_backrub_step ]

        sns.set_style("white")
        fig = plt.figure(
            figsize = [10.0, 8.5], dpi=600
        )
        fig.subplots_adjust( wspace = 0.6, hspace = 0.3)
        # fig.suptitle('$\Delta\Delta G$ prediction performance vs. number of structural ensemble members', fontsize=20)

        r_axes = []
        r_min = float('inf')
        r_max = float('-inf')

        mae_axes = []
        mae_min = float('inf')
        mae_max = float('-inf')

        # Information for caption/legend and output table
        ns = []
        best_step_ids = []
        legend_lines = []
        legend_labels = []
        data_table = []

        for ax_i, mut_type_subset in enumerate( mut_type_subsets ):
            ax = fig.add_subplot( 2, 2, ax_i + 1 )
            ax.set_ylabel(u"Pearson's R")
            ax.set_xlabel("Number of Models")
            rs = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby(['StructureOrder', 'ScoreMethodID'])[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
            # rs.to_csv( '/tmp/rs_structs-v-corr-%s-%s.csv' % (exp_run_name, mut_type_subset) )
            control_rs = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby(['StructureOrder'])[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()

            # Determine best correlation step
            argmax = rs['Experimental ddG'].argmax()
            best_step_id = rs.ix[argmax]['ScoreMethodID']
            if force_backrub_step != None and best_step_id > 10:
                best_step_id = force_backrub_step
            best_step_ids.append( best_step_id )

            ax.set_title( '(%s) - %s' % (string.ascii_lowercase[ax_i], mut_types[mut_type_subset]) )
            rs = rs.loc[ rs['ScoreMethodID'] == best_step_id ]

            maes = df.loc[ (df['ScoreMethodID'] == best_step_id) & (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('StructureOrder')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )
            control_maes = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby('StructureOrder')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )

            n = df.loc[ (df['ScoreMethodID'] == best_step_id) & (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('StructureOrder')[[pred_colname, exp_colname]].apply( lambda x: len( x[exp_colname] ) )
            assert( len( n.drop_duplicates() ) == 1 )
            ns.append( n.drop_duplicates().values[0] )

            r_min = min( r_min, min(rs['Experimental ddG']), min(control_rs['Experimental ddG']) )
            r_max = max( r_max, max(rs['Experimental ddG']), max(control_rs['Experimental ddG']) )

            ax2 = ax.twinx()
            ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
            ax2.set_ylabel("MAE")

            control_r, = ax.plot(
                control_rs['StructureOrder'], control_rs['Experimental ddG'],
                'o',
                color = run_colors[control_run],
                zorder = 1,
            )
            exp_r, = ax.plot(
                rs['StructureOrder'], rs['Experimental ddG'],
                'o',
                color = run_colors[exp_run_name],
                zorder = 2,
            )
            control_mae, = ax2.plot(
                control_maes.index, control_maes.values,
                'P',
                color = ( 1.0 + np.array(run_colors[control_run]) ) / 2.0,
                markersize = 5,
                zorder = 1,
            )
            exp_mae, = ax2.plot(
                maes.index, maes.values,
                'P',
                color = ( 1.0 + np.array(run_colors[exp_run_name]) ) / 2.0,
                markersize = 5,
                zorder = 2,
            )
            index_arrays = [np.array(x) for x in [control_rs['StructureOrder'].values, rs['StructureOrder'].values, control_maes.index.values, maes.index.values]]
            for i in index_arrays:
                for j in index_arrays:
                    np.testing.assert_array_equal( i, j )

            for num_structs_to_save, r_to_save, mae_to_save in zip(rs['StructureOrder'], rs['Experimental ddG'], maes.values):
                if num_structs_to_save in [1, 20, 30, 40, 50]:
                    data_table.append( (run_names[exp_run_name], mut_types[mut_type_subset], num_structs_to_save, r_to_save, mae_to_save) )
            for num_structs_to_save, r_to_save, mae_to_save in zip(control_rs['StructureOrder'], control_rs['Experimental ddG'], control_maes.values):
                if num_structs_to_save in [1, 20, 30, 40, 50]:
                    data_table.append( (run_names[control_run], mut_types[mut_type_subset], num_structs_to_save, r_to_save, mae_to_save) )

            if ax_i == 0:
                legend_lines.extend( [exp_r, control_r, exp_mae, control_mae] )
                legend_labels.extend( [
                    'R - ' + run_names[exp_run_name],
                    'R - ' + run_names[control_run],
                    'MAE - ' + run_names[exp_run_name],
                    'MAE - ' + run_names[control_run],
                ] )

            mae_min = min( mae_min, min(maes.values), min(control_maes.values) )
            mae_max = max( mae_max, max(maes.values), max(control_maes.values) )

            ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))

            r_axes.append( ax )
            mae_axes.append( ax2 )

        r_range = r_max - r_min
        r_min -= r_range * 0.1
        r_max += r_range * 0.1
        mae_range = mae_max - mae_min
        mae_min -= mae_range * 0.1
        mae_max += mae_range * 0.1

        for r_ax, mae_ax in zip(r_axes, mae_axes):
            r_ax.set_ylim( [r_min, r_max] )
            mae_ax.set_ylim( [mae_min, mae_max] )
            r_ax.set_zorder( mae_ax.get_zorder() + 10 )
            r_ax.patch.set_visible(False)

        if extra_fig_name != None:
            output_figure_name = 'structs-v-corr-%s-%s-%s' % (sorting_type, exp_run_name, extra_fig_name)
        else:
            output_figure_name = 'structs-v-corr-%s-%s' % (sorting_type, exp_run_name)
        output_figure_name = output_figure_name.replace('_', '-').replace('.', '')
        underlying_name = '%s-underlying-data' % output_figure_name
        out_path = os.path.join( output_fig_path, output_figure_name + '.pdf' )
        color_description = ''
        if exp_run_name == 'zemu_1.2-60000_rscript_simplified-t14':
            color_description += '\nPredictions generated with the Flex ddG protocol are shown in blue.'
        elif exp_run_name == 'ddg_monomer_16_003-zemu-2':
            color_description += '\nPredictions generated with the ddg\_monomer protocol are shown in purple.'
        color_description += '\nPredictions generated with the no backrub control protocol are shown in green.'
        if sorting_type == 'id' and exp_run_name == 'zemu_1.2-60000_rscript_simplified-t14':
            supp_note = ' This plot is analogous to \cref{fig:structs-v-corr-WildTypeComplex-zemu-12-60000-rscript-simplified-t14} in the main manuscript, except that the models used for averaging are not sorted by score.'
        else:
            supp_note = ''
        if force_backrub_step != None and best_step_id > 10:
            backrub_steps_note = ' Flex ddG is run with %d backrub steps.' % force_backrub_step
        else:
            backrub_steps_note = ''
        sub_dict = {
            'backrub-steps-note' : backrub_steps_note,
            'supp-note' : supp_note,
            'fig-label' : output_figure_name,
            'fig-path' : out_path,
            'exp-run-name' : run_names[exp_run_name],
            'control-name' : run_names[control_run],
            'sorting-type' : sorting_type_descriptions[sorting_type],
            'underlying-label' : 'tab:%s' % underlying_name,
            'color-description' : color_description,
        }
        for alpha_i, alpha in enumerate( string.ascii_lowercase[:4] ):
            if best_step_ids[alpha_i] >= 10 and force_backrub_step == None:
                sub_dict[ 'panel-' + alpha ] = '%s (n = %d, backrub steps = %d)' % ( mut_types[ mut_type_subsets[alpha_i] ].capitalize(),  ns[alpha_i], best_step_ids[alpha_i] )
            else:
                sub_dict[ 'panel-' + alpha ] = '%s (n = %d)' % ( mut_types[ mut_type_subsets[alpha_i] ].capitalize(),  ns[alpha_i] )

        leg = plt.figlegend(
            legend_lines, legend_labels, loc = 'lower center', labelspacing=0.0,
            markerscale = 1.3,
            fontsize = 13,
            ncol = 2,
        )
        fig.savefig( out_path )
        save_latex( 'latex_templates/structs-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )
        print out_path

        # Save underlying data
        data_table = pd.DataFrame.from_records(
            data_table,
            columns = ['Run', 'Subset', 'Models', 'R', 'MAE'],
        )
        latex_lines = data_table.to_latex( float_format = '%.2f', index = False ).split('\n')
        latex_lines = [r'\begin{table}'] + latex_lines + ['\caption[]{Selection of key data shown in \cref{fig:%s}}' % output_figure_name, '\label{tab:%s}' % underlying_name, r'\end{table}', '']
        with open( 'output/latex/%s.tex' % underlying_name, 'w' ) as f:
            f.write( '\n'.join(latex_lines) )

def table_main( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_03'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('zemu-values', 11, 'id_01'),
    ]

    short_caption = ""
    caption_text = short_caption + "Summary of prediction performance. Flex ddG predictions used 50 models and %d backrub steps. ddG monomer predictions were taken by averaging the \ddg\ scores of the three lowest scoring output models. N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps # TODo and dataset. ddg monomer from text

    subset_table(
        'table-main', results_df, display_runs, caption_text, short_caption, table_mut_types = [
            'complete', 's2l', 'ala', 'sing_ala', 'mult_mut', 'mult_all_ala', 'mult_none_ala', 'antibodies',
        ],
    )

def table_ref( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu_1.2-60000_rscript_validated-ref', backrub_steps, 'id_50'),
    ]

    short_caption = "REF results"
    caption_text = 'Performance comparison of the standard flex ddG protocol (using Rosetta\'s Talaris energy function) with flex ddG run with the REF score function, and %d backrub steps. Data for the flex ddG method with the Talaris energy function are as in \cref{tab:table-main} in the main text. res $<=$ 1.5 Ang." indicates data points for which the resolution of the input wild-type crystal structure is less than or equal to 1.5 \AA. N = number of cases in the dataset or subset. R = Pearson\'s R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold.' % backrub_steps

    table_mut_types = [
        'complete', 's2l', 'sing_ala', 'mult_mut',
        'res_lte15', 'res_gte25',
    ]

    subset_table( 'table-ref', results_df, display_runs, caption_text, short_caption, table_mut_types = table_mut_types )

def backrub_temp_table( results_df ):
    backrub_steps = 10000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu-brub_1.6-nt10000', backrub_steps, 'id_50'),
    ]

    short_caption = "Comparison of backrub temperature results"
    caption_text = "Flex ddG performance comparison, when backrub is run with a sampling temperature (kT) of 1.2 or 1.6 and %d backrub steps (this differs from \cref{table-main}, which is shown at 35,000 backrub steps). N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps

    subset_table( 'table-temperature', results_df, display_runs, caption_text, short_caption )

def ddg_monomer_table( results_df ):
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_03'),
        ('ddg_monomer_16_003-zemu-2', 7, 'WildTypeComplex_03'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('ddg_monomer_16_003-zemu-2', 7, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]


    short_caption = 'ddG monomer results'
    caption_text = short_caption + ". N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold."

    subset_table( 'table-ddG-monomer', results_df, display_runs, caption_text, short_caption )

def multiple_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Multiple mutations results'
    caption_text = short_caption + " (backrub steps = %d). N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps

    subset_table( 'table-mult', results_df, display_runs, caption_text, short_caption, table_mut_types = ['mult_mut', 'mult_all_ala', 'mult_none_ala', 'ala'] )

def antibodies_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on antibodies'
    caption_text = "Performance of the Rosetta flex ddG method on the subset of complexes containing an antibody binding partner (flex ddG run with %d backrub steps). N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps

    subset_table( 'table-antibodies', results_df, display_runs, caption_text, short_caption, table_mut_types = ['complete', 'antibodies'] )

def stabilizing_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on stabilizing mutations'
    caption_text = "Performance of the Rosetta flex ddG method on the subset of mutations experimentally determined to be stabilizing ($\Delta\Delta G <= -1$), neutral ($-1 < \Delta\Delta G < 1$), or destabilizing ($\Delta\Delta G >= 1$). Flex ddG was run with %d backrub steps. N = number of cases in the dataset or subset. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps

    subset_table( 'table-stabilizing', results_df, display_runs, caption_text, short_caption, table_mut_types = ['stabilizing', 'neutral', 'positive'] )

def by_pdb_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_simplified-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on PDB '
    caption_text = ". Backrub steps = %d. N = number of cases (variants) for each complex. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct. Best performance for each metric and dataset is shown in bold." % backrub_steps

    display_columns = collections.OrderedDict( [
        ('MutTypes', 'PDB'),
        ('PredictionRun', 'Prediction Method'),
        ('N', 'N'),
        # ('StructureOrder', 'Num/Sorting of Structs'),
        ('R', 'R'),
        ('MAE_unscaled', 'MAE'),
        ('FractionCorrect_unscaled', 'FC'),
    ] )

    all_middle_lines = []
    for pdb_mut_type in sorted( [ m for m in results_df['MutTypes'].drop_duplicates().values if m.startswith('pdb-') ] ):
        filled_short_caption = short_caption + pdb_mut_type[4:]
        header_lines, group_rows_lines, footer_lines = subset_table( 'table-' + pdb_mut_type, results_df, display_runs, filled_short_caption + caption_text, filled_short_caption, table_mut_types = [pdb_mut_type] )
        all_middle_lines.append( group_rows_lines )

    new_lines = []
    new_lines.extend( header_lines )
    for i, rows in enumerate(all_middle_lines):
        new_lines.extend( rows )
        if i < len(all_middle_lines) - 1:
            new_lines.append( '\hline' )
    new_lines.extend( footer_lines )

    short_caption = 'Flex ddG performance by PDB structure for all complexes with 5 or more cases'
    table_name = 'table-by-structure'
    save_latex(
        'latex_templates/subset-table.tex',
        {
            'table-text' : '\n'.join(new_lines),
            'caption' : short_caption + caption_text,
            'short-caption' : short_caption,
            'table-label' : table_name,
        },
        out_tex_name = table_name,
        long_table = True,
        wrap_command = '\\small',
    )

def subset_table( table_name, results_df, display_runs, caption_text, short_caption, table_mut_types = None, display_columns = None ):
    if table_mut_types == None:
        table_mut_types = display_mut_types

    # Determine if StructureOrder needs to be included as a column
    unique_structure_orders = {}
    for run_name, step, structure_order in display_runs:
        key_tup = (run_name, step)
        if key_tup not in unique_structure_orders:
            unique_structure_orders[key_tup] = set()
        unique_structure_orders[key_tup].add( structure_order )
    include_structure_order_column = False
    for s in unique_structure_orders.values():
        if len( s ) > 1:
            include_structure_order_column = True
            break

    if display_columns == None:
        if include_structure_order_column:
            display_columns = collections.OrderedDict( [
                ('MutTypes', 'Mutation Category'),
                ('PredictionRun', 'Prediction Method'),
                ('N', 'N'),
                ('StructureOrder', 'Num/Sorting of Structs'),
                ('R', 'R'),
                ('MAE_unscaled', 'MAE'),
                ('FractionCorrect_unscaled', 'FC'),
            ] )
        else:
            display_columns = collections.OrderedDict( [
                ('MutTypes', 'Mutation Category'),
                ('PredictionRun', 'Prediction Method'),
                ('N', 'N'),
                ('R', 'R'),
                ('MAE_unscaled', 'MAE'),
                ('FractionCorrect_unscaled', 'FC'),
            ] )


    # These run names will have step explanations added in, if appropriate
    annotated_run_names = copy.copy( run_names )

    results_subset = pd.DataFrame()
    for mut_type_i, mut_type in enumerate(table_mut_types):
        for run_name, step, structure_order in display_runs:
            # print run_name, results_df.loc[results_df['PredictionRun'] == run_name].head()
            # print step, results_df.loc[results_df['Step'] == step].head()
            # print structure_order, results_df.loc[ results_df['StructureOrder'] == structure_order ].head()
            # print mut_type, results_df.loc[results_df['MutTypes'] == mut_type].head()
            new_row = results_df.loc[ (results_df['PredictionRun'] == run_name) & (results_df['Step'] == step) & (results_df['StructureOrder'] == structure_order ) & (results_df['MutTypes'] == mut_type) ]

            # Annotate run name further if necessary
            if 'ddg_monomer' in run_name and step == 7:
                new_run_name = run_name + '-7'
                annotated_run_names[new_run_name] = annotated_run_names[run_name] + ' (soft-rep)'
                new_row = new_row.replace( { run_name : new_run_name } )
            elif 'ddg_monomer' in run_name and step == 8:
                new_run_name = run_name + '-8'
                annotated_run_names[new_run_name] = annotated_run_names[run_name]
                new_row = new_row.replace( { run_name : new_run_name } )

            if len(new_row) != 1:
                print run_name, step, structure_order, mut_type
                print new_row.head()
                assert( len(new_row) == 1 )
            if len(results_subset) == 0:
                results_subset = new_row
            else:
                results_subset = results_subset.append( new_row )

            # Extra caption text (only for first loop through mut_types)
            # if mut_type_i == 0:
            #     if step >= 50:
            #         caption_text += ' %s steps = %d.' % ( run_names[run_name], step )


    out_path = os.path.join( output_fig_path, '%s.csv' % table_name.replace('-', '_') )
    beautified_results = results_subset[ display_columns.keys() ].rename( columns = display_columns ).replace( annotated_run_names ).replace( mut_types )
    beautified_results.to_csv(out_path, float_format = '%.2f' )

    # Generate latex and add hlines
    latex_lines = beautified_results.to_latex( index=False, float_format = '%.2f' ).split('\n')
    header_lines = []
    group_rows = [ [] ]
    footer_lines = []
    in_group_lines = 0
    for line in latex_lines:
        line = line.strip()
        if len(line) == 0:
            continue

        if in_group_lines == 0:
            header_lines.append( line )
            if 'midrule' in line:
                in_group_lines = 1
        elif in_group_lines == 1:
            if 'bottomrule' in line or 'end{tabular}' in line:
                footer_lines.append( line )
            else:
                group_rows[-1].append( line )
                if len( group_rows[-1] ) >= len(display_runs):
                    group_rows.append( [] )
        elif in_group_lines == 2:
            footer_lines.append( line )

    if len(group_rows[-1]) == 0:
        group_rows.pop()

    # Bold best numeric in each column of each set of group rows

    # better_columns defines which metric is "best", since a higher R is better, whereas a lower MAE is better
    # These and first_numeric_columns will need to be changed if display_columns is changed
    lower_better_columns = [ display_columns.keys().index('MAE_unscaled') ]
    higher_better_columns = [ display_columns.keys().index('R'), display_columns.keys().index('FractionCorrect_unscaled') ]
    first_numeric_column = display_columns.keys().index('R')
    new_group_rows = []
    for rows in group_rows:
        new_rows = []

        best_row_results = {}
        col_vals = {}
        # First pass to find best
        for row in rows:
            for i, result in enumerate( [x.strip().replace('\\', '').strip() for x in row.split('&')]):
                if i >= first_numeric_column and (
                        ( i in lower_better_columns and (i not in best_row_results or float(result) <= float(best_row_results[i])) )
                        or
                        ( i in higher_better_columns and (i not in best_row_results or float(result) >= float(best_row_results[i])) )
                ):
                    best_row_results[i] = result
                if i not in col_vals:
                    col_vals[i] = set()
                col_vals[i].add( result )

        # Second pass to add bolding
        for row_i, row in enumerate(rows):
            row_data = [x.strip().replace('\\', '').strip() for x in row.split('&')]
            new_row_data = []
            for i, cell in enumerate( row_data ):
                if i < first_numeric_column and len( col_vals[i] ) == 1:
                    if row_i == 0:
                        new_row_data.append( ' \multirow{ %d}{*}{%s} ' % (len(rows), cell) )
                    else:
                        new_row_data.append( ' ' )
                else:
                    if i in best_row_results and best_row_results[i] in cell:
                        new_row_data.append(' %s ' % cell.replace( best_row_results[i], '\\textbf{' + best_row_results[i] + '}' ) )
                    else:
                        new_row_data.append( ' %s ' %cell )
            new_rows.append( '&'.join( new_row_data ).replace('_', '\\_') + ' \\\\' )
        new_group_rows.append( new_rows )

    group_rows = new_group_rows

    new_lines = []
    group_rows_lines = []
    new_lines.extend( header_lines )
    for i, rows in enumerate(group_rows):
        new_lines.extend( rows )
        group_rows_lines.extend( rows )
        if i < len(group_rows) - 1:
            new_lines.append( '\hline' )
            group_rows_lines.append( '\hline' )
    new_lines.extend( footer_lines )

    save_latex(
        'latex_templates/subset-table.tex',
        {
            'table-text' : '\n'.join(new_lines),
            'caption' : caption_text,
            'short-caption' : short_caption,
            'table-label' : table_name
        },
        out_tex_name = table_name,
    )

    # print table_name
    # print beautified_results.head( n = 30 )
    # print
    return (header_lines, group_rows_lines, footer_lines)

def prediction_error( score_method_id = 35000, prediction_run = 'zemu_1.2-60000_rscript_simplified-t14' ):
    # Outputs individual errors for other figure purposes
    df = load_df()
    df = df.loc[ df['MutType'] == 'complete' ]
    df['Error'] = df['total'] - df['ExperimentalDDG']
    df['AbsError'] = np.abs( df['Error'] )
    df.sort_values( ['AbsError'], inplace = True, ascending = False )
    df.to_csv( os.path.join(output_fig_path, 'prediction_error.csv') )

    error_df = df.loc[ ( df['StructureOrder'] == 'id_50') & (df['ScoreMethodID'] == score_method_id)  & (df['PredictionRunName'] == prediction_run) ]
    # print error_df.head()
    print prediction_run, len(error_df)
    error_by_subset = []
    for subset, dataset_ids in subsets.iteritems():
        if not subset.startswith('pdb-'):
            mean_error = df.loc[ df['DataSetID'].isin(dataset_ids) ]['Error'].mean()
            error_by_subset.append( ( mean_error, len( dataset_ids ), subset ) )
    error_by_subset.sort()
    for mean_error, dataset_ids_len, subset in error_by_subset:
        print subset, dataset_ids_len, '%.2f' % mean_error

def make_supp_csv( only_fig_columns = True ):
    df = load_df()

    fc = figure_components()

    display_run_names = [
        'zemu_1.2-60000_rscript_simplified-t14',
        'zemu_1.2-60000_rscript_validated-ref',
        'zemu-brub_1.6-nt10000',
        'zemu_control-69aa526',
        'ddg_monomer_16_003-zemu-2',
        'zemu-values',
        'tal_GAM',
        'ref_GAM',
        'control_GAM',
    ]
    display_mut_types = ['s2l', 'ala', 'sing_ala', 'mult_mut', 'mult_all_ala', 'mult_none_ala', 'antibodies']
    structure_orders = {
        'id_01' : '1 randomly selected model',
        'id_30' : '3 randomly selected models',
        'id_30' : '30 randomly selected models',
        'id_50' : 'All 50 models',
        'WildTypeComplex_01' : 'lowest scoring model by wild-type complex energy',
        'WildTypeComplex_03' : '3 lowest scoring models by wild-type complex energy',
        'WildTypeComplex_30' : '30 lowest scoring models by wild-type complex energy',
        'WildTypeComplex_50' : 'All 50 models',
    }
    backrub_steps = range(8,21) + [2500, 10000, 35000, 50000]
    column_name_map = collections.OrderedDict( [
        ('PredictionRunName', 'Run'),
        ('DataSetID', 'ID'),
        ('ScoreMethodID', 'Backrub Steps'),
        ('total', 'Total Score'),
        ('ExperimentalDDG', 'Experimental ddG'),
        ('StructureOrder', 'Number of Models'),
        ('MutType', 'Subset'),
        ('PredictionID', 'PredictionID'),
    ] )
    supp_df = df.loc[
        ( df['PredictionRunName'].isin(display_run_names) ) &
        ( df['ScoreMethodID'].isin(backrub_steps) ) &
        ( df['StructureOrder'].isin(structure_orders.keys()) ) &
        ( df['MutType'] == 'complete' )
    ]

    supp_df_renamed = supp_df[ column_name_map.keys() ]
    supp_df_renamed['ScoreMethodID'] = supp_df_renamed['ScoreMethodID'].map(lambda x: int(x) if x > 100 else None)
    supp_df_renamed['PredictionRunName'] = supp_df_renamed['PredictionRunName'].map(lambda x: run_names[x])
    supp_df_renamed['StructureOrder'] = supp_df_renamed['StructureOrder'].map(lambda x: structure_orders[x])
    supp_df_renamed.columns = column_name_map.values()
    supp_df_renamed.to_csv( os.path.join(output_dir, 'supp-df-rows.csv') )

    rev_df = pd.read_csv( 'table-zemu-filtered.csv' )[ ['DataSetID', 'PDBFileID', 'Mutations', 'ExperimentalDDG'] ]
    for mut_type in display_mut_types:
        rev_df[ 'Subset: ' + mut_types[mut_type] + '?' ] = rev_df['DataSetID'].isin( subsets[mut_type] )

    rev_df = rev_df.append(
        pd.DataFrame( {'DataSetID':[-1]} ),
        ignore_index = True,
    )
    rev_df.sort_values( ['DataSetID'], inplace = True )
    supp_df.sort_values( ['PredictionRunName', 'StructureOrder', 'ScoreMethodID'], inplace = True )
    for name, group in supp_df.groupby(['PredictionRunName', 'ScoreMethodID', 'StructureOrder']):
        run_name, score_method, structure_order = name
        fig_components_text = make_fig_components_text(fc, run_name, score_method, structure_order)
        if only_fig_columns and len(fig_components_text.strip()) == 0:
            continue
        run_name = run_names[run_name]
        if score_method > 100:
            score_method = '%05d' % int(score_method)
        elif 'GAM' in run_name:
            score_method = '35000' # HARD CODE
        else:
            score_method = 'N/A'
        structure_order = structure_orders[structure_order]
        name = (run_name, score_method, structure_order)

        score_df = group[ ['DataSetID', 'total'] ]
        new_col = 'Run: %s - Backrub steps: %s - Number of models: %s' % name
        if new_col in rev_df.columns:
            print name, new_col
            print rev_df.columns
            assert( new_col not in rev_df.columns )
        score_df.columns = ['DataSetID', new_col ]
        score_df = score_df.append( pd.DataFrame( {
            'DataSetID' : [-1],
            new_col : [fig_components_text],
        } ) )
        rev_df = rev_df.merge( score_df, left_on = 'DataSetID', right_on = 'DataSetID', how = 'left', suffixes = ('','_y') )#.drop( 'DataSetID_y', axis = 1 )

    #rev_df.sort_values( ['PredictionRunName', 'StructureOrder', 'ScoreMethodID'], inplace = True )
    new_cols = []
    for col in rev_df.columns:
        if col in column_name_map:
            new_cols.append( column_name_map[col] )
        else:
            new_cols.append( col )
    rev_df.columns = new_cols
    rev_df.to_csv( os.path.join(output_dir, 'flex-ddG-data.csv') )

def make_subsets_report():
    df = pd.DataFrame( [(mut_types[key] if key in mut_types else key, len(value)) for key, value in subsets.iteritems()], columns = ['subset', 'count'] )
    df.sort_values( 'count', ascending = False, inplace = True )
    df.to_csv( 'subset_report.csv' )

if __name__ == '__main__':
    steps_vs_corr( 'steps-v-corr_burial', ['buried_lt0.10_mean', 'buried_gte0.10_mean', 'buried_lt0.25_mean', 'buried_gte0.25_mean'] )
    steps_vs_corr( 'steps-v-corr_burial2', ['buried_lt0.10_mean', 'buried_gte0.10_mean', 'buried_lt0.05_mean', 'buried_gte0.05_mean'] )
    steps_vs_corr( 'steps-v-corr_ss', ['SS-all-None', 'SS-all-H', 'SS-all-E', 'SS-all-T'] )
    steps_vs_corr( 'steps-v-corr_ss2', ['SS-all-B', 'SS-all-G', 'SS-all-None', 'SS-all-T'] )
    figure_structs_vs_corr(
        mut_type_subsets = ['buried_lt0.10_mean', 'buried_gte0.10_mean', 'buried_lt0.25_mean', 'buried_gte0.25_mean'],
        extra_fig_name = 'burial',
    )
    figure_structs_vs_corr(
        mut_type_subsets = ['SS-all-None', 'SS-all-H', 'SS-all-E', 'SS-all-T'],
        extra_fig_name = 'SS',
    )
    sys.exit()
    figure_mult_non_ala()
    make_supp_csv()

    prediction_error()

    table_composition()
    table_versions()
    figure_scatter()
    make_subsets_report()
    steps_vs_corr( 'steps-v-corr', ['complete', 's2l', 'mult_none_ala', 'sing_ala'] )
    steps_vs_corr( 'steps-v-corr_mult', ['mult_mut', 'ala', 'mult_all_ala', 'mult_none_ala'] )
    steps_vs_corr( 'steps-v-corr_resolution', ['complete', 'res_gte25', 'res_lte15', 'res_gt15_lt25'] )
    steps_vs_corr( 'steps-v-corr_some_sizes', ['some_s2l', 's2l', 'some_l2s', 'l2s'] )
    figure_structs_vs_corr()
    figure_structs_vs_corr( 'ddg_monomer_16_003-zemu-2' )

    results_df = make_results_df()

    by_pdb_table( results_df )
    table_ref( results_df )
    table_main( results_df )
    backrub_temp_table( results_df )
    ddg_monomer_table( results_df )
    multiple_table( results_df )
    antibodies_table( results_df )
    stabilizing_table( results_df )

    # compile_latex()
