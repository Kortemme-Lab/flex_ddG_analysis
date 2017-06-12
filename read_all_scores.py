import cPickle as pickle
import pandas as pd
import numpy as np
import os
import sys
from stats import fraction_correct, mae as calc_mae, bootstrap_xy_stat
import scipy

csv_path = os.path.expanduser( '~/data/ddg/interface_ddg_paper/publication_analyze_output/60k_with_control/split/id.csv.gz' )
assert( os.path.isfile(csv_path) )

def main():
    df = pd.read_csv(csv_path)

    mut_types = sorted( df['MutType'].drop_duplicates().values )
    structure_orders = sorted( df['StructureOrder'].drop_duplicates().values )
    fcs = []
    maes = []
    rs = []
    slopes = []
    df_mut_types = []
    df_structure_orders = []
    prediction_runs = []

    for prediction_run in sorted( df['PredictionRunName'].drop_duplicates().values ):
        for mut_type in mut_types:
            for structure_order in structure_orders:
                sub_df = df.loc[ (df['MutType'] == mut_type) & (df['StructureOrder'] == structure_order) & (df['PredictionRunName'] == prediction_run) ]
                print mut_type, structure_order, prediction_run

                fc = fraction_correct(
                    sub_df['total'].values,
                    sub_df['ExperimentalDDG'].values,
                )
                print 'FC:', '%.3f' % fc
                fcs.append(fc)

                mae = calc_mae( sub_df['total'], sub_df['ExperimentalDDG'] )
                print 'MAE:', '%.3f' % mae
                maes.append(mae)

                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress( sub_df['total'], sub_df['ExperimentalDDG'] )
                print 'R:', '%.3f' % r_value
                rs.append(r_value)
                print 'slope:', '%.3f' % slope
                slopes.append(slope)
                df_mut_types.append(mut_type)
                df_structure_orders.append(structure_order)
                prediction_runs.append(prediction_run)
                print

    results_df = pd.DataFrame( {
        'PredictionRun' : prediction_runs,
        'MutTypes' : df_mut_types,
        'StructureOrder' : df_structure_orders,
        'FractionCorrect' : fcs,
        'MAE' : maes,
        'Slope' : rs,
        'R' : rs,
    } )
    sort_cols = ['FractionCorrect', 'R', 'MAE']
    ascendings = [False, False, True]
    for sort_col, asc in zip(sort_cols, ascendings):
        results_df.sort_values( sort_col, inplace = True, ascending = asc )
        print results_df[ ['PredictionRun', 'MutTypes', 'StructureOrder', sort_col] ].head(n=10)

if __name__ == '__main__':
    main()
