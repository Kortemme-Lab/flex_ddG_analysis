import pandas as pd
import argparse
import getpass
import MySQLdb
import os
import sys
import time
import numpy as np
import json
import cPickle as pickle

def get_df( mysql_con ):
    query = "SELECT ID, UserPPDataSetExperimentID  FROM PredictionPPI WHERE PredictionSet='zemu_1.2-60000_rscript_simplified-t14'"
    full_df = pd.read_sql_query( query, mysql_con)

    query = "SELECT ID, UserPPDataSetExperimentID  FROM PredictionPPI WHERE PredictionSet='zemu_1.2-60000_struct-t14'"
    new_df = pd.read_sql_query( query, mysql_con )

    data_df = pd.read_csv( 'data/zemu_1.2-60000_rscript_simplified-t14-id_30.csv.gz' )
    data_df = data_df.loc[ data_df['ScoreMethodID'] == 35000 ]

    data_df['abs_error'] = np.abs( data_df['total'] - data_df['ExperimentalDDG'] )
    data_df.sort_values( 'abs_error', inplace = True )
    data_ids = data_df[['DataSetID', 'abs_error']]
    data_ids = data_ids.merge( new_df, left_on = 'DataSetID', right_on = 'UserPPDataSetExperimentID' )

    return data_ids

def main( mysql_con ):
    df = get_df( mysql_con )
    df.to_csv( 'interesting.csv' )

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

    mysql_con = MySQLdb.connect( host = args.host,
                                 port = args.port,
                                 user = args.user, passwd = args.password,
                                 db = args.db
    )

    main( mysql_con )
