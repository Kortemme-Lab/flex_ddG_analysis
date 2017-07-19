import pandas as pd
import argparse
import getpass
import MySQLdb
import os
import time

from run_descriptions import all_runs

# Comment out structure orders below to make selection faster
structure_orders = [
    'id', # Order structures for averaging based on the random order they were created
    'WildTypeComplex', # Order structures for averaging by the wild type complex score
    # 'MutantComplex', # Order structures for averaging by the mutant type complex score
]

dataframe_cache = '/dbscratch/kyleb/new_query_cache'
categorized_query_cache = os.path.join( dataframe_cache, 'categorized_queries' )
if not os.path.isdir( categorized_query_cache ):
    os.makedirs( categorized_query_cache )

def fetch_from_db_and_reorder( mysql_host, mysql_port, mysql_user, mysql_pass, mysql_db ):
    with open('scores_select.sql', 'r') as f:
        all_scores_query = ' '.join( [ line.strip() for line in f.readlines() if not line.startswith('#') ] )

    mysql_con = MySQLdb.connect( host = mysql_host,
                                 port = mysql_port,
                                 user = mysql_user, passwd = mysql_pass,
                                 db = mysql_db
    )

    # Load data and determine common DataSetIDs
    for benchmark_run in all_runs:
        step_multiplication_factor = benchmark_run.step_multiplication_factor
        prediction_set_name = benchmark_run.prediction_set_name
        print '\n%s' % str(prediction_set_name)

        cached_df_path = os.path.join( dataframe_cache, prediction_set_name + '_all-scores.csv.gz' )
        if os.path.isfile(cached_df_path):
            print 'Cached main df exists:', cached_df_path
            df = pd.read_csv(cached_df_path, compression = 'gzip' )
        else:
            df = pd.read_sql_query( all_scores_query % prediction_set_name, mysql_con)
            # Multiply by step factor
            if step_multiplication_factor:
                df.loc[:,'ScoreMethodID'] *= step_multiplication_factor

            df.to_csv(cached_df_path, compression = 'gzip' )

        print 'PDBs in set: %d' % len( df[['PDBFileID']].drop_duplicates() )
        print 'Datapoints length (without duplicates): %d' % len( df[['DataSetID']].drop_duplicates() )

        # Reorder structure ID's, if specified
        for structure_order in structure_orders:
            # Check for cache
            cached_query_df_path = os.path.join( dataframe_cache, prediction_set_name + '-' + structure_order + '.csv.gz' )
            if os.path.isfile( cached_query_df_path ):
                print 'Cached ordered df exists:', cached_query_df_path
            else:
                if structure_order == 'id':
                    new_query_df = df.assign( StructureOrder = 'id' )
                elif structure_order in df['ScoreType'].drop_duplicates().values:
                    new_df = df.copy()
                    new_df.set_index(['PredictionID', 'ScoreMethodID', 'StructureID'], inplace=True)
                    sorted_df = df.loc[ df['ScoreType'] == structure_order ].sort_values(
                        by = ['PredictionID', 'ScoreMethodID', 'total'],
                    )
                    new_series_values = {
                        'NewStructureID' : [],
                        'PredictionID' : [],
                        'ScoreMethodID' : [],
                        'StructureID' : [],
                    }
                    row_count = 1
                    last_id_tup = None
                    last_prediction_id = None
                    last_score_method_id = None
                    for index, row in sorted_df.iterrows():
                        row_id_tup = (row['PredictionID'], row['ScoreMethodID'])
                        if last_id_tup == None:
                            last_id_tup = row_id_tup
                        if last_id_tup != row_id_tup:
                            row_count = 1
                            last_id_tup = row_id_tup
                        new_series_values['NewStructureID'].append( row_count )
                        new_series_values['PredictionID'].append( row['PredictionID'] )
                        new_series_values['ScoreMethodID'].append( row['ScoreMethodID'] )
                        new_series_values['StructureID'].append( row['StructureID'] )
                        row_count += 1
                    new_column_df = pd.DataFrame( new_series_values )
                    new_column_df.set_index(['PredictionID', 'ScoreMethodID', 'StructureID'], inplace=True)
                    new_df = new_df.merge(new_column_df, how = 'left', left_index = True, right_index = True )
                    new_query_df = new_df.reset_index().assign( StructureOrder = structure_order )
                    new_query_df.drop('StructureID', axis=1, inplace=True)
                    new_query_df.rename( columns = {'NewStructureID' : 'StructureID'}, inplace = True )
                new_query_df.to_csv( cached_query_df_path, compression = 'gzip' )

            assert( 'MutType' not in new_query_df.columns )

    mysql_con.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--host', default = 'localhost',
                        help = 'MySQL host name')
    parser.add_argument('--port', type = int, default = 3306,
                        help = 'MySQL port')
    parser.add_argument('--user', default = getpass.getuser(),
                        help = 'MySQL user')
    parser.add_argument('--password', default = '',
                        help = 'MySQL password')
    parser.add_argument('--db', default = 'DDG',
                        help = 'MySQL database schema')

    args = parser.parse_args()
    fetch_from_db_and_reorder( args.host, args.port, args.user, args.password, args.db)
