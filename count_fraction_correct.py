import pandas as pd
import numpy as np

from read_all_scores import add_score_categories, mut_types

input_data= [
    ('data/zemu_control-69aa526-noglypivot-id_50.csv.gz', 8),
    ('data/zemu_1.2-60000_rscript_simplified-t14-id_50.csv.gz', 35000),
    # ('data/zemu_1.2-60000_rscript_validated-t14-id_50.csv.gz', 35000),
    # ('data/zemu_control-69aa526-id_50.csv.gz', 8),
]

for csv_path, backrub_steps in input_data:
    print csv_path

    df = pd.read_csv( csv_path )
    df = add_score_categories( pd.read_csv( csv_path ) )
    df = df.loc[ (df['ScoreMethodID'] == backrub_steps) ]
    df = df.drop_duplicates()
    for mut_type in ['complete', 's2l']:
        sub_df = df.loc[ df['MutType'] == mut_type ]
        if mut_type == 'complete':
            assert( len(sub_df) == 1240 )
        print '%s subset (n=%d, backrub steps=%d)' % (mut_types[mut_type], len(sub_df), backrub_steps)

        sub_df = sub_df.loc[ sub_df['ExperimentalDDG'] > 0.0 ]
        print 'Number of mutations with exp. ddG > 0: %d' % len( sub_df )
        print 'Number of the above predicted to be <= -1REU: %d' % len( sub_df.loc[ sub_df['total'] <= -1.0 ] )
        print 'Number of the above predicted to be < 0 REU: %d' % len( sub_df.loc[ sub_df['total'] < 0.0 ] )
        # print 'Number of the above predicted to be > 0 REU: %d' % len( sub_df.loc[ sub_df['total'] > 0.0 ] )
        # print 'Number of the above predicted to be 0 <= DDG < 1: %d' % len( sub_df.loc[ (sub_df['total'] >= 0.0) & (sub_df['total'] < 1.0) ] )
        # print 'Number of the above predicted to be >= 1 REU: %d' % len( sub_df.loc[ (sub_df['total'] >= 1.0) ] )
        print

        sub_df = df.loc[ df['MutType'] == mut_type ]
        sub_df = sub_df.loc[ sub_df['ExperimentalDDG'] >= 1.0 ]
        print 'Number of mutations with exp. ddG >= 1.0: %d' % len( sub_df )
        print 'Number of the above predicted to be <= -1REU: %d' % len( sub_df.loc[ sub_df['total'] <= -1.0 ] )
        print 'Number of the above predicted to be < 0 REU: %d' % len( sub_df.loc[ sub_df['total'] < 0.0 ] )
        # print 'Number of the above predicted to be > 0 REU: %d' % len( sub_df.loc[ sub_df['total'] > 0.0 ] )
        # print 'Number of the above predicted to be 0 <= DDG < 1: %d' % len( sub_df.loc[ (sub_df['total'] >= 0.0) & (sub_df['total'] < 1.0) ] )
        # print 'Number of the above predicted to be >= 1 REU: %d' % len( sub_df.loc[ (sub_df['total'] >= 1.0) ] )
        print
