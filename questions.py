import pandas as pd
import numpy as np

# hdf_path = '/dbscratch/kyleb/new_query_cache/zemu_1.2-60000_rscript_simplified-t14-id.hdf'

# #df = pd.read_csv( '/dbscratch/kyleb/new_query_cache/zemu_1.2-60000_rscript_simplified-t14-id.csv.gz' )
# # df.to_hdf( hdf_path, 'df' )
# df = pd.read_hdf( hdf_path )
# print df.head()
# print df['ScoreType'].drop_duplicates()

df = pd.read_csv( '/dbscratch/kyleb/new_query_cache/zemu_control-69aa526-id.csv.gz' )
df = df.loc[ (df['ScoreType'] == 'MutantComplex') | (df['ScoreType'] == 'WildTypeComplex') ]
results = {}
count_50 = 0
count_45 = 0
for name, group in df.groupby( ['PredictionID', 'ScoreType', 'total'] ):
    i, struct_type, score = name
    if struct_type not in results:
        results[struct_type] = {}
    if i not in results[struct_type]:
        results[struct_type][i] = []
    results[struct_type][i].append( len(group) )
for struct_type in results:
    for i in results[struct_type]:
        results[struct_type][i] = len( results[struct_type][i] )
    results[struct_type] = np.mean( results[struct_type].values() )

for struct_type, avg_num in results.iteritems():
    print 'There are an average of %.1f unique structures per %s model in the control' % (avg_num, struct_type)

# df = pd.read_csv( '/dbscratch/kyleb/new_query_cache/zemu_control-69aa526_all-scores.csv.gz' )
# print df.head()
