import cPickle as pickle
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import os
import sys
import shutil
from stats import fraction_correct, mae as calc_mae, bootstrap_xy_stat
import scipy

csv_path = os.path.expanduser( '~/data/ddg/interface_ddg_paper/publication_analyze_output/60k_with_control/split/id.csv.gz' )
output_dir = 'output'
generate_plots = False
print_statistics = False

if generate_plots:
    # Import here as they can be slow, and are unneeded if plots aren't going to be made
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns

assert( os.path.isfile(csv_path) )
if os.path.isdir(output_dir):
    shutil.rmtree(output_dir)

def add_score_categories(df):
    stabilizing = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] <= -1.0) ].copy()
    stabilizing.loc[:,'MutType'] = 'stabilizing'

    neutral = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] > -1.0) & (df['ExperimentalDDG'] < 1.0) ].copy()
    neutral.loc[:,'MutType'] = 'neutral'

    positive = df.loc[ (df['MutType'] == 'complete') & (df['ExperimentalDDG'] >= 1.0) ].copy()
    positive.loc[:,'MutType'] = 'positive'

    df = df.append( stabilizing )
    df = df.append( positive )
    df = df.append( neutral )
    return df

def main():
    df = add_score_categories( pd.read_csv(csv_path) )

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

                    if not os.path.isdir(sub_output_dir):
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
        print results_df[ ['PredictionRun', 'MutTypes', 'StructureOrder', 'Step', sort_col] ].head(n=20)
    results_csv_path = os.path.join(output_dir, 'results.csv')
    results_df.sort_values('R', ascending = False).to_csv( results_csv_path )
    print results_csv_path

if __name__ == '__main__':
    main()
