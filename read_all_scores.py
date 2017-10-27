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

# TODO: move this into repository
csv_paths = [
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-t14-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-ref-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-ref-WildTypeComplex_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_control-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_control-69aa526-noglypivot-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_control-69aa526-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/ddg_monomer_16_003-zemu-2-WildTypeComplex_03.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/ddg_monomer_16_003-zemu-2-WildTypeComplex_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu-brub_1.6-nt10000-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu-values-id_01.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-ref-cart-id_50.csv.gz' ),
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
current_palette = sns.color_palette()

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
    'mult_all_ala' : 'Multiple, all alanine',
    'mult_none_ala' : 'Multiple, none alanine',
    's2l' : 'Small-To-Large Mutation(s)',
    'l2s' : 'Large-To-Small',
    'ala' : 'Mutation(s) to alanine',
    'sing_ala' : 'Single mutation to alanine',
    'res_gte25' : 'Res. $>=$ 2.5 Ang.',
    'res_lte15' : 'Res. $<=$ 1.5 Ang.',
    'res_gt15_lt25' : '1.5 Ang. $<$ Res. $<$ 2.5 Ang.',
    'some_s2l' : 'Some Small-To-Large',
    'some_l2s' : 'Some Large-To-Small',
    'antibodies' : 'Antibodies',
    'stabilizing' : 'Stabilizing',
    'neutral' : 'Neutral',
    'positive' : 'Positive',
}

run_names = {
    'zemu_1.2-60000_rscript_validated-t14' : 'flex ddG',
    'zemu_1.2-60000_rscript_validated-ref' : 'flex ddG (REF energy)',
    'zemu-brub_1.6-nt10000' : 'flex ddG (1.6 kT)',
    'ddg_monomer_16_003-zemu-2' : 'ddG monomer',
    'zemu_control' : 'no backrub control',
    'zemu-values' : 'ZEMu paper',
    'zemu_control-69aa526-noglypivot' : 'no backrub control',
    'zemu_control-69aa526' : 'no backrub control',
}

sorting_type_descriptions = {
    'id' : 'no sorting of structures',
    'WildTypeComplex' : 'structures sorted by minimized wild-type complex energy',
}


cached_loaded_df_initialized = False
cached_loaded_df = None
def load_df():
    global cached_loaded_df_initialized
    global cached_loaded_df
    if cached_loaded_df_initialized:
        return cached_loaded_df.copy()

    print 'Performing initial .csv load'

    df = pd.read_csv( csv_paths[0] )
    for csv_path in csv_paths[1:]:
        df = df.append( pd.read_csv( csv_path ) )
    df = add_score_categories( df )
    df = df.drop_duplicates( ['PredictionRunName', 'DataSetID', 'PredictionID', 'ScoreMethodID', 'MutType', 'total', 'ExperimentalDDG', 'StructureOrder'] )
    cached_loaded_df_initialized = True
    cached_loaded_df = df.copy()
    print 'Done loading csvs\n'
    return df

def add_score_categories(df, mut_type_subsets = None):
    if mut_type_subsets == None or 'complete' in mut_type_subsets:
        df = df.assign( MutType = 'complete' )

    if mut_type_subsets == None or 'stabilizing' in mut_type_subsets:
        stabilizing = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] <= -1.0) ].copy()
        stabilizing.loc[:,'MutType'] = 'stabilizing'
        df = df.append( stabilizing )

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
            subset_df.loc[:,'MutType'] = subset_name
            df = df.append( subset_df )

    return df

def save_latex( latex_template_file, sub_dict, out_tex_name = None ):
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

    sort_cols = ['FractionCorrect', 'R', 'MAE']
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

def figure_scatter():
    exp_run_name = 'zemu_1.2-60000_rscript_validated-t14'
    control_run_name = 'zemu_control-69aa526'
    point_size = 4.5
    alpha = 0.6
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.9 }

    df = load_df()
    exp_colname = 'Experimental ddG'
    pred_colname = 'Rosetta Score'
    top_subset = 'complete'
    bottom_subset = 's2l'
    df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
    df = df.rename( columns = {'total' : pred_colname} )

    fig = plt.figure(
        figsize=(8.5, 8.5), dpi=600
    )

    complete_corrs = df.loc[ (df['MutType'] == top_subset) & (df['PredictionRunName'] == exp_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_a = complete_corrs.index[0][0]
    df_a = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == top_subset) & (df['ScoreMethodID'] == best_step_a) ]
    ax1 = fig.add_subplot( 2, 2, 1 )

    complete_corrs = df.loc[ (df['MutType'] == top_subset) & (df['PredictionRunName'] == control_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_b = complete_corrs.index[0][0]
    df_b = df.loc[ (df['PredictionRunName'] == control_run_name) & (df['MutType'] == top_subset) & (df['ScoreMethodID'] == best_step_b) ]
    ax2 = fig.add_subplot( 2, 2, 2 )

    complete_corrs = df.loc[ (df['MutType'] == bottom_subset) & (df['PredictionRunName'] == exp_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_c = complete_corrs.index[0][0]
    df_c = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == bottom_subset) & (df['ScoreMethodID'] == best_step_c) ]
    ax3 = fig.add_subplot( 2, 2, 3 )

    complete_corrs = df.loc[ (df['MutType'] == bottom_subset) & (df['PredictionRunName'] == control_run_name) ].groupby( 'ScoreMethodID' )[[pred_colname,exp_colname]].corr().ix[0::2,exp_colname].sort_values( ascending = False )
    best_step_d = complete_corrs.index[0][0]
    df_d = df.loc[ (df['PredictionRunName'] == control_run_name) & (df['MutType'] == bottom_subset) & (df['ScoreMethodID'] == best_step_d) ]
    ax4 = fig.add_subplot( 2, 2, 4 )

    # Override b, c, d best steps to just be the same as a
    best_step_b = best_step_a
    best_step_c = best_step_a
    best_step_d = best_step_a

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
        color = sns.color_palette()[0],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_b, ax = ax2,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_c, ax = ax3,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )
    sns.regplot(
        y = pred_colname, x = exp_colname,
        data = df_d, ax = ax4,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )

    # ax2.set_xticklabels([])
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    # ax2.set_yticklabels([])
    ax4.set_ylabel('')

    ax1.set_title( '(a) Backrub - %s' % (mut_types[top_subset]) )
    ax2.set_title( '(b) Control - %s' % (mut_types[top_subset] ) )
    ax3.set_title( '(c) Backrub - %s' % (mut_types[bottom_subset]) )
    ax4.set_title( '(d) Control - %s' % (mut_types[bottom_subset]) )

    # Assert that these lengths are equal since N is displayed for top and bottom together
    assert( len(df_a) == len(df_b) )
    assert( len(df_c) == len(df_d) )

    out_path = os.path.join( output_fig_path, 'fig-scatter.pdf' )
    sub_dict = {
        'exp-method-name' : run_names[exp_run_name].capitalize(),
        'numsteps-a' : str( best_step_a ),
        'numsteps-c' : str( best_step_c ),
        'top-subset' : mut_types[top_subset].capitalize(),
        'top-n' : str( len(df_a) ),
        'control-method-name' : run_names[control_run_name].capitalize(),
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

def steps_vs_corr( output_figure_name, mut_type_subsets, control_run = 'zemu_control-69aa526', force_control_in_axes = True ):
    exp_run_name = 'zemu_1.2-60000_rscript_validated-t14'
    point_size = 4.5
    alpha = 0.6
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.9 }

    df = load_df()
    exp_colname = 'Experimental ddG'
    pred_colname = 'Rosetta Score'
    df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
    df = df.rename( columns = {'total' : pred_colname} )

    sns.set_style("white")
    fig = plt.figure(
        figsize=(10.0, 8.5), dpi=600
    )
    fig.subplots_adjust( wspace = 0.6, hspace = 0.3)


    r_axes = []
    r_min = float('inf')
    r_max = float('-inf')

    mae_axes = []
    mae_min = float('inf')
    mae_max = float('-inf')

    ns = []

    for ax_i, mut_type_subset in enumerate( mut_type_subsets ):
        ax = fig.add_subplot( 2, 2, ax_i + 1 )
        ax.set_title( mut_types[mut_type_subset] )
        ax.set_ylabel("Pearson's R")
        ax.yaxis.label.set_color(current_palette[0] )
        ax.set_xlabel("Backrub Step")
        ns.append( len(df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) & (df['ScoreMethodID'] == df['ScoreMethodID'].drop_duplicates().values[0]) ]) )
        rs = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
        maes = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )

        ax.plot( rs['ScoreMethodID'], rs['Experimental ddG'], 'o', color = current_palette[0] )
        r_min = min( r_min, min(rs['Experimental ddG']) )
        r_max = max( r_max, max(rs['Experimental ddG']) )

        ax2 = ax.twinx()
        ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
        ax2.set_ylabel("MAE")
        ax2.plot( maes.index, maes.values, 's', color = current_palette[1] )
        mae_min = min( mae_min, min(maes.values) )
        mae_max = max( mae_max, max(maes.values) )

        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
        ax2.yaxis.label.set_color(current_palette[1] )

        # Get control performance values
        control_rs = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
        assert( len(control_rs) == 1 )
        control_r = control_rs['Experimental ddG'][0]
        ax.plot( [0], [control_r], marker = 'X', color = current_palette[0] )

        control_maes = df.loc[ (df['PredictionRunName'] == control_run) & (df['MutType'] == mut_type_subset ) ].groupby('ScoreMethodID')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )
        assert( len(control_maes.values) == 1 )
        control_mae = control_maes.values[0]
        ax2.plot( [0], [control_mae], marker = 'X', color = current_palette[1] )

        if force_control_in_axes:
            r_min = min( r_min, control_r )
            r_max = max( r_max, control_r )
            mae_min = min( mae_min, control_mae )
            mae_max = max( mae_max, control_mae )

        r_axes.append( ax )
        mae_axes.append( ax2 )

    r_range = r_max - r_min
    r_min -= r_range * 0.1
    r_max += r_range * 0.1
    mae_range = mae_max - mae_min
    mae_min -= mae_range * 0.1
    mae_max += mae_range * 0.1

    for ax in r_axes:
        ax.set_ylim( [r_min, r_max] )

    for ax in mae_axes:
        ax.set_ylim( [mae_min, mae_max] )

    out_path = os.path.join( output_fig_path, '%s.pdf' % output_figure_name )
    sub_dict = {
        'panel-a' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[0] ].capitalize(),  ns[0] ),
        'panel-b' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[1] ].capitalize(),  ns[1] ),
        'panel-c' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[2] ].capitalize(),  ns[2] ),
        'panel-d' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[3] ].capitalize(),  ns[3] ),
        'fig-label' : output_figure_name,
        'fig-path' : out_path,
    }

    fig.savefig( out_path )
    save_latex( 'latex_templates/steps-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )

    fig.savefig( out_path )
    print out_path

def figure_structs_vs_corr( exp_run_name = 'zemu_1.2-60000_rscript_validated-t14', force_backrub_step = 35000 ):
    sorting_types = ['WildTypeComplex', 'id']
    base_path = '/dbscratch/kyleb/new_query_cache/summed_and_averaged/%s-%s_%02d.csv.gz'
    control_base_path = '/dbscratch/kyleb/new_query_cache/summed_and_averaged/%s-%s_%02d.csv.gz'
    number_of_structures = 50
    mut_type_subsets = ['complete', 's2l', 'mult_none_ala', 'sing_ala']
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

        sns.set_style("white")
        fig = plt.figure(
            figsize=(10.0, 8.5), dpi=600
        )
        fig.subplots_adjust( wspace = 0.6, hspace = 0.3)
        fig.suptitle('Number of Structures Performance\n(%s)' % sorting_type_descriptions[sorting_type], fontsize=20)

        r_axes = []
        r_min = float('inf')
        r_max = float('-inf')

        mae_axes = []
        mae_min = float('inf')
        mae_max = float('-inf')

        # Information for caption
        ns = []
        best_step_ids = []

        for ax_i, mut_type_subset in enumerate( mut_type_subsets ):
            ax = fig.add_subplot( 2, 2, ax_i + 1 )
            ax.set_ylabel("Pearson's R")
            ax.yaxis.label.set_color(current_palette[0] )
            ax.set_xlabel("Number of Structures")
            rs = df.loc[ (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby(['StructureOrder', 'ScoreMethodID'])[[pred_colname, exp_colname]].corr().ix[0::2, exp_colname].reset_index()
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

            ax2.plot( control_maes.index, control_maes.values, 's', color = (1.0 + np.array( current_palette[1] )) / 2.0 )
            ax2.plot( maes.index, maes.values, 's', color = current_palette[1] )
            ax.plot( control_rs['StructureOrder'], control_rs['Experimental ddG'], 'o', color = (1.0 + np.array( current_palette[0] )) / 2.0 )
            ax.plot( rs['StructureOrder'], rs['Experimental ddG'], 'o', color = current_palette[0] )

            mae_min = min( mae_min, min(maes.values), min(control_maes.values) )
            mae_max = max( mae_max, max(maes.values), max(control_maes.values) )

            ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
            ax2.yaxis.label.set_color(current_palette[1] )

            r_axes.append( ax )
            mae_axes.append( ax2 )

        r_range = r_max - r_min
        r_min -= r_range * 0.1
        r_max += r_range * 0.1
        mae_range = mae_max - mae_min
        mae_min -= mae_range * 0.1
        mae_max += mae_range * 0.1

        for ax in r_axes:
            ax.set_ylim( [r_min, r_max] )

        for ax in mae_axes:
            ax.set_ylim( [mae_min, mae_max] )

        output_figure_name = 'structs-v-corr-%s-%s' % (sorting_type, exp_run_name)
        output_figure_name = output_figure_name.replace('_', '-').replace('.', '')
        out_path = os.path.join( output_fig_path, output_figure_name + '.pdf' )
        sub_dict = {
            'fig-label' : output_figure_name,
            'fig-path' : out_path,
            'exp-run-name' : run_names[exp_run_name],
            'control-name' : run_names[control_run],
        }
        for alpha_i, alpha in enumerate( string.ascii_lowercase[:4] ):
            if best_step_ids[alpha_i] >= 10:
                sub_dict[ 'panel-' + alpha ] = '%s (n = %d, backrub step = %d)' % ( mut_types[ mut_type_subsets[alpha_i] ].capitalize(),  ns[alpha_i], best_step_ids[alpha_i] )
            else:
                sub_dict[ 'panel-' + alpha ] = '%s (n = %d)' % ( mut_types[ mut_type_subsets[alpha_i] ].capitalize(),  ns[alpha_i] )

        fig.savefig( out_path )
        save_latex( 'latex_templates/structs-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )
        print out_path

def table_main( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_03'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('zemu-values', 11, 'id_01'),
    ]

    short_caption = ""
    caption_text = short_caption + "Main results table. Backrub steps = %d. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    subset_table(
        'table-main', results_df, display_runs, caption_text, short_caption, table_mut_types = [
            'complete', 's2l', 'ala', 'sing_ala', 'mult_mut', 'mult_all_ala', 'mult_none_ala', 'antibodies',
        ],
    )

def table_ref( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu_1.2-60000_rscript_validated-ref', backrub_steps, 'id_50'),
    ]

    short_caption = "REF results"
    caption_text = short_caption + ". Backrub steps = %d. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    table_mut_types = [
        'complete', 's2l', 'sing_ala', 'mult_mut',
        'res_lte15', 'res_gte25',
    ]

    subset_table( 'table-ref', results_df, display_runs, caption_text, short_caption, table_mut_types = table_mut_types )

def backrub_temp_table( results_df ):
    backrub_steps = 10000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu-brub_1.6-nt10000', backrub_steps, 'id_50'),
    ]

    short_caption = "Comparison of backrub temperature results"
    caption_text = short_caption + ". Backrub steps = %d. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

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
    caption_text = short_caption + ". R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct."

    subset_table( 'table-temperature', results_df, display_runs, caption_text, short_caption )

def multiple_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Multiple mutations results'
    caption_text = short_caption + " (backrub steps = %d). R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    subset_table( 'table-mult', results_df, display_runs, caption_text, short_caption, table_mut_types = ['mult_mut', 'mult_all_ala', 'mult_none_ala', 'ala'] )

def antibodies_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on antibodies'
    caption_text = "Performance of the Rosetta flex ddG method on the subset of complexes containing an antibody binding partner (backrub steps = %d). R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    subset_table( 'table-antibodies', results_df, display_runs, caption_text, short_caption, table_mut_types = ['complete', 'antibodies'] )

def stabilizing_table( results_df ):
    backrub_steps = 35000
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on stabilizing mutations'
    caption_text = "Performance of the Rosetta flex ddG method on the subset of mutations experimentally determined to be stabilizing ($\Delta\Delta$G $< 0$). Backrub steps = %d. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    subset_table( 'table-stabilizing', results_df, display_runs, caption_text, short_caption, table_mut_types = ['stabilizing', 'neutral', 'positive'] )

def by_pdb_table( results_df ):
    backrub_steps = 32500
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', backrub_steps, 'id_50'),
        ('zemu_control-69aa526', 8, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    short_caption = 'Flex ddG performance on PDB '
    caption_text = ". Backrub steps = %d. R = Pearson's R. MAE = Mean Absolute Error. FC = Fraction Correct." % backrub_steps

    display_columns = collections.OrderedDict( [
        ('MutTypes', 'PDB'),
        ('PredictionRun', 'Prediction Method'),
        ('N', 'N'),
        # ('StructureOrder', 'Num/Sorting of Structs'),
        ('R', 'R'),
        ('MAE', 'MAE'),
        ('FractionCorrect', 'FC'),
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

    short_caption = 'Flex ddG performance by structure'
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
                ('MAE', 'MAE'),
                ('FractionCorrect', 'FC'),
            ] )
        else:
            display_columns = collections.OrderedDict( [
                ('MutTypes', 'Mutation Category'),
                ('PredictionRun', 'Prediction Method'),
                ('N', 'N'),
                ('R', 'R'),
                ('MAE', 'MAE'),
                ('FractionCorrect', 'FC'),
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
                annotated_run_names[new_run_name] = annotated_run_names[run_name] + ' (hard-rep)'
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
    lower_better_columns = [ display_columns.keys().index('MAE') ]
    higher_better_columns = [ display_columns.keys().index('R'), display_columns.keys().index('FractionCorrect') ]
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

def prediction_error( score_method_id = 35000, prediction_run = 'zemu_1.2-60000_rscript_validated-t14' ):
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

if __name__ == '__main__':
    results_df = make_results_df()

    prediction_error()

    table_composition()
    table_versions()
    figure_scatter()
    steps_vs_corr( 'steps-v-corr', ['complete', 's2l', 'mult_none_ala', 'sing_ala'] )
    steps_vs_corr( 'steps-v-corr_mult', ['mult_mut', 'ala', 'mult_all_ala', 'mult_none_ala'] )
    steps_vs_corr( 'steps-v-corr_resolution', ['complete', 'res_gte25', 'res_lte15', 'res_gt15_lt25'] )
    steps_vs_corr( 'steps-v-corr_some_sizes', ['some_s2l', 's2l', 'some_l2s', 'l2s'] )
    figure_structs_vs_corr()
    figure_structs_vs_corr( 'ddg_monomer_16_003-zemu-2' )

    by_pdb_table( results_df )
    table_ref( results_df )
    table_main( results_df )
    backrub_temp_table( results_df )
    ddg_monomer_table( results_df )
    multiple_table( results_df )
    antibodies_table( results_df )
    stabilizing_table( results_df )

    compile_latex()
