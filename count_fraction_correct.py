import pandas as pd
import numpy as np

from read_all_scores import add_score_categories, mut_types

csv_path = '/dbscratch/kyleb/new_query_cache/summed_and_averaged/zemu_1.2-60000_rscript_validated-t14-id_50.csv.gz'
backrub_steps = 35000

df = pd.read_csv( csv_path )
df = add_score_categories( pd.read_csv( csv_path ) )
df = df.loc[ (df['ScoreMethodID'] == backrub_steps) ]

for mut_type in ['complete', 's2l', 'l2s']:
    sub_df = df.loc[ df['MutType'] == mut_type ]
    print '%s subset (n=%d, backrub steps=%d)' % (mut_types[mut_type], len(sub_df), backrub_steps)
    sub_df = sub_df.loc[ sub_df['ExperimentalDDG'] <= -1.0 ]
    print 'Number of stabilizing mutations in this subset: %d' % len( sub_df )
    print 'Number of the above predicted to be <= -1REU: %d' % len( sub_df.loc[ sub_df['total'] <= -1.0 ] )
    print 'Number of the above predicted to be < 0 REU: %d' % len( sub_df.loc[ sub_df['total'] < 0.0 ] )
    print 'Number of the above predicted to be 0 <= DDG < 1: %d' % len( sub_df.loc[ (sub_df['total'] >= 0.0) & (sub_df['total'] < 1.0) ] )
    print 'Number of the above predicted to be >= 1 REU: %d' % len( sub_df.loc[ (sub_df['total'] >= 1.0) ] )
    print
