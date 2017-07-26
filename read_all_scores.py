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

csv_paths = [
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-t14-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_control-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/ddg_monomer_16_003-zemu-2-WildTypeComplex_03.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu-brub_1.6-nt10000-id_50.csv.gz' ),
    os.path.expanduser( '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu-values-id_01.csv.gz' ),
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
    assert( os.path.isfile(csv_path) )
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

display_mut_types = [
    'complete', 'ala', 'sing_mut', 'mult_mut', 's2l',
]

mut_types = {
    'complete' : 'Complete dataset',
    'sing_mut' : 'Single mutation',
    'mult_mut' : 'Multiple mutations',
    's2l' : 'Small-To-Large Mutation(s)',
    'l2s' : 'Large-To-Small',
    'ala' : 'Mutation(s) to Alanine',
    'res_gte25' : 'Res. >= 2.5$\AA$',
    'res_lte15' : 'Res. <= 1.5$\AA$',
    'res_gt15_lt25' : '1.5$\AA$ < Res. < 2.5$\AA$',
    'some_s2l' : 'Some Small-To-Large',
    'some_l2s' : 'Some Large-To-Small',
}

run_names = {
    'zemu_1.2-60000_rscript_validated-t14' : 'ddG backrub',
    'ddg_monomer_16_003-zemu-2' : 'ddG monomer',
    'zemu_control' : 'minimization only control',
    'zemu-values' : 'ZEMu paper',
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

def add_score_categories(df):
    df = df.assign( MutType = 'complete' )

    stabilizing = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] <= -1.0) ].copy()
    stabilizing.loc[:,'MutType'] = 'stabilizing'

    neutral = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] > -1.0) & (df['ExperimentalDDG'] < 1.0) ].copy()
    neutral.loc[:,'MutType'] = 'neutral'

    positive = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] >= 1.0) ].copy()
    positive.loc[:,'MutType'] = 'positive'

    df = df.append( stabilizing )
    df = df.append( positive )
    df = df.append( neutral )

    for subset_name, subset_keys in subsets.iteritems():
        subset_df = df.loc[ (df['MutType'] == 'complete') & (df['DataSetID'].isin(subset_keys)) ].copy()
        subset_df.loc[:,'MutType'] = subset_name
        df = df.append( subset_df )

    return df

def save_latex( latex_template_file, sub_dict, out_tex_name = None ):
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

def compile_latex():
    shutil.copy( main_template_latex_file, latex_output_dir )

    for i in xrange( 0, 3 ):
        output = subprocess.check_output( ['xelatex', os.path.basename( main_template_latex_file )[:-4] ], cwd = latex_output_dir )

def make_results_df( generate_plots = False, print_statistics = False ):
    df = load_df()

    # Add 'complete' to front of list so scaling is calculated first
    mut_types = sorted( df['MutType'].drop_duplicates().values )
    mut_types.insert(0, mut_types.pop( mut_types.index('complete') ) )

    structure_orders = sorted( df['StructureOrder'].drop_duplicates().values )
    fcs_unscaled = []
    maes_unscaled = []
    fcs = []
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

                    fc_unscaled = fraction_correct(
                        sub_df['total'].values,
                        sub_df['ExperimentalDDG'].values,
                    )
                    fcs_unscaled.append(fc_unscaled)
                    fc = fraction_correct(
                        sub_df['scaled_total'].values,
                        sub_df['ExperimentalDDG'].values,
                    )
                    fcs.append(fc)

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
                        fig = plt.figure(figsize=(8.5, 8.5), dpi=600)
                        ax = fig.add_subplot(1, 1, 1)
                        sns.regplot(x="total", y="ExperimentalDDG", data = sub_df, ax=ax, scatter_kws = { 's' : 4.5 } )

                        sns.regplot(
                            x="total", y="ExperimentalDDG", data = sub_df, ax=outer_ax,
                            color = ( float(num_steps - step_i) / num_steps, 0, ( float(step_i) / num_steps ) ),
                            scatter_kws = { 's' : 2 },
                        )

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
        'FractionCorrect_unscaled' : fcs_unscaled,
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
    results_csv_path = os.path.join(output_dir, 'results.csv')
    new_column_order = [
        'PredictionRun', 'N', 'MutTypes', 'Step', 'StructureOrder',
        'R', 'MAE', 'MAE_unscaled', 'FractionCorrect', 'FractionCorrect_unscaled',
        'Slope', 'ScalingIntercepts', 'ScalingSlopes',
    ]
    assert( set(new_column_order) == set(results_df.columns) )
    results_df = results_df[new_column_order]
    results_df.sort_values('R', ascending = False).to_csv( results_csv_path )
    if print_statistics:
        print results_csv_path
        print

    return results_df

def figure_2():
    exp_run_name = 'zemu_1.2-60000_rscript_validated-t14'
    control_run_name = 'zemu_control'
    point_size = 4.5
    alpha = 0.6
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.9 }

    df = load_df()
    exp_colname = 'Experimental DDG'
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
        x = pred_colname, y = exp_colname,
        data = df_a, ax = ax1,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )
    sns.regplot(
        x = pred_colname, y = exp_colname,
        data = df_b, ax = ax2,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )
    sns.regplot(
        x = pred_colname, y = exp_colname,
        data = df_c, ax = ax3,
        scatter_kws = scatter_kws,
        line_kws = line_kws,
        ci = None,
        color = sns.color_palette()[0],
    )
    sns.regplot(
        x = pred_colname, y = exp_colname,
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

    ax1.set_title( '(a) Backrub (%d steps) - %s' % (best_step_a, mut_types[top_subset]) )
    ax2.set_title( '(b) Control - %s' % (mut_types[top_subset] ) )
    ax3.set_title( '(c) Backrub (%d steps) - %s' % (best_step_c, mut_types[bottom_subset]) )
    ax4.set_title( '(d) Control - %s' % (mut_types[bottom_subset]) )

    # Assert that these lengths are equal since N is displayed for top and bottom together
    assert( len(df_a) == len(df_b) )
    assert( len(df_c) == len(df_d) )

    out_path = os.path.join( output_fig_path, 'fig2.pdf' )
    sub_dict = {
        'exp-method-name' : run_names[exp_run_name].capitalize(),
        'numsteps-a' : str( best_step_a ),
        'numsteps-c' : str( best_step_c ),
        'top-subset' : mut_types[top_subset].capitalize(),
        'top-n' : str( len(df_a) ),
        'control-method-name' : run_names[control_run_name].capitalize(),
        'bottom-subset' : mut_types[bottom_subset].capitalize(),
        'bottom-n' : str( len(df_c) ),
        'fig-path' : os.path.relpath(out_path, latex_output_dir),
    }

    fig.savefig( out_path )
    save_latex( 'latex_templates/figure-2.tex', sub_dict )
    print out_path

def table_1():
    # Dataset composition
    df = load_df()
    df = df.drop_duplicates( ['PredictionRunName', 'DataSetID', 'PredictionID', 'ScoreMethodID', 'MutType', 'total', 'ExperimentalDDG', 'StructureOrder'] )
    control_df = df.loc[ (df['PredictionRunName'] == 'zemu_control') & (df['ScoreMethodID'] == 8 ) ]

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

    print 'Table 1:'
    print table_df.head( n = 20 )
    save_latex( 'latex_templates/table-1.tex', { 'table-1' : table_df.to_latex( index=False ) } )
    table_df.to_csv( os.path.join(output_fig_path, 'table_1.csv') )

def steps_vs_corr( output_figure_name, mut_type_subsets ):
    exp_run_name = 'zemu_1.2-60000_rscript_validated-t14'
    point_size = 4.5
    alpha = 0.6
    scatter_kws = { 's' : point_size, 'alpha' : alpha }
    line_kws = { 'linewidth' : 0.9 }

    df = load_df()
    exp_colname = 'Experimental DDG'
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

        ax.plot( rs['ScoreMethodID'], rs['Experimental DDG'], 'o', color = current_palette[0] )
        r_min = min( r_min, min(rs['Experimental DDG']) )
        r_max = max( r_max, max(rs['Experimental DDG']) )

        ax2 = ax.twinx()
        ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
        ax2.set_ylabel("MAE")
        ax2.plot( maes.index, maes.values, 's', color = current_palette[1] )
        mae_min = min( mae_min, min(maes.values) )
        mae_max = max( mae_max, max(maes.values) )

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

    out_path = os.path.join( output_fig_path, '%s.pdf' % output_figure_name )
    sub_dict = {
        'panel-a' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[0] ].capitalize(),  ns[0] ),
        'panel-b' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[1] ].capitalize(),  ns[1] ),
        'panel-c' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[2] ].capitalize(),  ns[2] ),
        'panel-d' : '%s (n=%d)' % ( mut_types[ mut_type_subsets[3] ].capitalize(),  ns[3] ),
        'fig-label' : output_figure_name,
        'fig-path' : os.path.relpath(out_path, latex_output_dir),
    }

    fig.savefig( out_path )
    save_latex( 'latex_templates/steps-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )

    fig.savefig( out_path )
    print out_path

def figure_4():
    exp_run_name = 'zemu_1.2-60000_rscript_validated-t14'
    sorting_types = ['WildTypeComplex', 'id']
    base_path = '/dbscratch/kyleb/new_query_cache/summed_and_averaged/%s-%s_%02d.csv.gz'
    number_of_structures = 50

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

            structure_orders.append( '%s_%02d' % (sorting_type, struct_id) )

        df['StructureOrder'] = df['StructureOrder'].apply(
            lambda s: int( s.split('_')[1] )
        )

        df = add_score_categories(df)

        point_size = 4.5
        alpha = 0.6
        scatter_kws = { 's' : point_size, 'alpha' : alpha }
        line_kws = { 'linewidth' : 0.9 }

        exp_colname = 'Experimental DDG'
        pred_colname = 'Rosetta Score'
        df = df.rename( columns = {'ExperimentalDDG' : exp_colname} )
        df = df.rename( columns = {'total' : pred_colname} )

        sns.set_style("white")
        fig = plt.figure(
            figsize=(10.0, 8.5), dpi=600
        )
        fig.subplots_adjust( wspace = 0.6, hspace = 0.3)
        fig.suptitle('Number of Structures Performance (%s)' % sorting_type_descriptions[sorting_type], fontsize=20)

        mut_type_subsets = ['complete', 's2l', 'sing_mut', 'ala']

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

            # Determine best correlation step
            argmax = rs['Experimental DDG'].argmax()
            best_step_id = rs.ix[argmax]['ScoreMethodID']
            best_step_ids.append( best_step_id )
            ax.set_title( '%s' % (mut_types[mut_type_subset]) )
            rs = rs.loc[ rs['ScoreMethodID'] == best_step_id ]

            maes = df.loc[ (df['ScoreMethodID'] == best_step_id) & (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('StructureOrder')[[pred_colname, exp_colname]].apply( lambda x: calc_mae( x[exp_colname], x[pred_colname] ) )

            n = df.loc[ (df['ScoreMethodID'] == best_step_id) & (df['PredictionRunName'] == exp_run_name) & (df['MutType'] == mut_type_subset ) ].groupby('StructureOrder')[[pred_colname, exp_colname]].apply( lambda x: len( x[exp_colname] ) )
            assert( len( n.drop_duplicates() ) == 1 )
            ns.append( n.drop_duplicates().values[0] )

            ax.plot( rs['StructureOrder'], rs['Experimental DDG'], 'o', color = current_palette[0] )
            r_min = min( r_min, min(rs['Experimental DDG']) )
            r_max = max( r_max, max(rs['Experimental DDG']) )

            ax2 = ax.twinx()
            ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
            ax2.set_ylabel("MAE")
            ax2.plot( maes.index, maes.values, 's', color = current_palette[1] )
            mae_min = min( mae_min, min(maes.values) )
            mae_max = max( mae_max, max(maes.values) )

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

        output_figure_name = 'fig4-%s' % sorting_type
        out_path = os.path.join( output_fig_path, output_figure_name + '.pdf' )
        sub_dict = {
            'panel-a' : '%s (n = %d, backrub step = %d)' % ( mut_types[ mut_type_subsets[0] ].capitalize(),  ns[0], best_step_ids[0] ),
            'panel-b' : '%s (n = %d, backrub step = %d)' % ( mut_types[ mut_type_subsets[1] ].capitalize(),  ns[1], best_step_ids[1] ),
            'panel-c' : '%s (n = %d, backrub step = %d)' % ( mut_types[ mut_type_subsets[2] ].capitalize(),  ns[2], best_step_ids[2] ),
            'panel-d' : '%s (n = %d, backrub step = %d)' % ( mut_types[ mut_type_subsets[3] ].capitalize(),  ns[3], best_step_ids[3] ),
            'fig-label' : output_figure_name,
            'fig-path' : os.path.relpath(out_path, latex_output_dir),
        }

        fig.savefig( out_path )
        save_latex( 'latex_templates/structs-vs-corr.tex', sub_dict, out_tex_name = output_figure_name )
        print out_path

def table_2( results_df ):
    # PredictionRun, Step, StructureOrder
    display_runs = [
        ('zemu_1.2-60000_rscript_validated-t14', 35000, 'id_50'),
        ('ddg_monomer_16_003-zemu-2', 8, 'WildTypeComplex_03'),
        ('zemu_control', 8, 'id_50'),
        ('zemu-values', 11, 'id_01'),
    ]
    display_columns = {
        'PredictionRun' : 'Prediction Method',
        'N' : 'N',
        'MutTypes' : 'Mutation Category',
        'R' : 'R',
        'MAE' : 'MAE',
        'FractionCorrect' : 'FC',
    }
    display_column_order = ['PredictionRun', 'MutTypes', 'R', 'MAE', 'FractionCorrect']

    results_subset = pd.DataFrame()
    for mut_type in display_mut_types:
        for run_name, step, structure_order in display_runs:
            new_row = results_df.loc[ (results_df['PredictionRun'] == run_name) & (results_df['Step'] == step) & (results_df['StructureOrder'] == structure_order ) & (results_df['MutTypes'] == mut_type) ]
            if len(new_row) != 1:
                print run_name, step, structure_order, mut_type
                print new_row.head()
                assert( len(new_row) == 1 )
            if len(results_subset) == 0:
                results_subset = new_row
            else:
                results_subset = results_subset.append( new_row )

    out_path = os.path.join( output_fig_path, 'table_2.csv' )
    beautified_results = results_subset[ display_column_order ].rename( columns = display_columns ).replace( run_names ).replace( mut_types )
    beautified_results.to_csv(out_path, float_format = '%.2f' )

    # Generate latex and add hlines
    latex_lines = beautified_results.to_latex( index=False, float_format = '%.2f' ).split('\n')
    header_lines = []
    group_rows = [ [] ]
    footer_lines = []
    in_group_lines = 0
    for line in latex_lines:
        if in_group_lines == 0:
            header_lines.append( line )
            if line.startswith( '\midrule' ):
                in_group_lines = 1
        elif in_group_lines == 1:
            if line.startswith( '\bottomrule' ):
                footer_lines.append( line )
            else:
                group_rows[-1].append( line )
                if len( group_rows[-1] ) >= len(display_runs):
                    group_rows.append( [] )
        elif in_group_lines == 2:
            footer_lines.append( line )

    if len(group_rows[-1]) == 0:
        group_rows.pop()

    new_lines = []
    new_lines.extend( header_lines )
    for i, rows in enumerate(group_rows):
        new_lines.extend( rows )
        if i < len(group_rows) - 1:
            new_lines.append( '\hline' )
    new_lines.extend( footer_lines )

    save_latex( 'latex_templates/table-2.tex', {
        'table-2' : '\n'.join(new_lines),
    } )

    print 'Table 2:'
    print beautified_results.head( n = 30 )
    print


if __name__ == '__main__':
    # results_df = make_results_df()
    # table_1()
    # table_2( results_df )
    # figure_2()
    # steps_vs_corr( 'fig3', ['complete', 's2l', 'mult_mut', 'ala'] )
    # steps_vs_corr( 'fig3_resolution', ['complete', 'res_gte25', 'res_lte15', 'res_gt15_lt25'] )
    # steps_vs_corr( 'fig3_some_sizes', ['some_s2l', 's2l', 'some_l2s', 'l2s'] )
    figure_4()
    compile_latex()
