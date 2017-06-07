import cPickle as pickle
import pandas
import numpy as np
import os
import sys
from stats import fraction_correct, mae, fraction_correct_fuzzy_linear, bootstrap_xy_stat
import scipy

csv_path = os.path.expanduser( '~/data/ddg/interface_ddg_paper/publication_analyze_output/zemu_12-60000_rscript_validated-t14/joint/all_scores.csv' )
assert( os.path.isfile(csv_path) )

def main():
    df = pandas.read_csv(csv_path)

    mut_types = df['MutType'].drop_duplicates().values
    structure_orders = [ df['StructureOrder'].drop_duplicates().values[0] ]

    for mut_type in mut_types:
        for structure_order in structure_orders:
            sub_df = df.loc[ (df['MutType'] == mut_type) & (df['StructureOrder'] == structure_order) ]
            print mut_type, structure_order
            print 'FC:', fraction_correct(
                sub_df['total'].values,
                sub_df['ExperimentalDDG'].values,
            )

            print bootstrap_xy_stat(
                sub_df['total'].values,
                sub_df['ExperimentalDDG'].values,
                fraction_correct,
            )

            print 'FC fuzzy:', fraction_correct_fuzzy_linear(
                sub_df['total'].values,
                sub_df['ExperimentalDDG'].values,
            )

            print 'MAE:', mae( sub_df['total'], sub_df['ExperimentalDDG'] )
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress( sub_df['total'], sub_df['ExperimentalDDG'] )
            print 'R:', r_value
            print 'slope:', slope
            print get_xy_dataset_statistics_pandas(
                sub_df,
                x_series = 'total',
                y_series = 'ExperimentalDDG',
                bootstrap_data = True
            )
            print

if __name__ == '__main__':
    main()
