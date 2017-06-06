import cPickle as pickle
import pandas
import numpy as np
import os
import sys
from klab.stats.misc import fraction_correct_pandas

pickle_path = os.path.expanduser( '~/data/ddg/interface_ddg_paper/publication_analyze_output/zemu_12-60000_rscript_validated-t14/joint/all_scores.pickle' )
assert( os.path.isfile(pickle_path) )

def main():
    with open(pickle_path, 'rb') as f:
        df = pickle.load(f)
    print df.head()
    print
    print df['MutType'].drop_duplicates()
    print df['StructureOrder'].drop_duplicates()
    print df.loc[ df['DataSetID'] == 1 ].head()
    for mut_type, structure_order in [('complete', 'id')]:
        print mut_type, structure_order
        print fraction_correct_pandas(
            df.loc[ (df['MutType'] == mut_type) & (df['StructureOrder'] == structure_order) ],
            x_series = 'total',
            y_series = 'ExperimentalDDG',
        )
        print

if __name__ == '__main__':
    main()
