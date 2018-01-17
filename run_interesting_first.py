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

def main( mysql_con ):
    query = "SELECT ID, UserPPDataSetExperimentID  FROM PredictionPPI WHERE PredictionSet='zemu_1.2-60000_rscript_simplified-t14'"
    full_df = pd.read_sql_query( query, mysql_con)

    query = "SELECT ID, UserPPDataSetExperimentID  FROM PredictionPPI WHERE PredictionSet='zemu_1.2-60000_struct-t14'"
    new_df = pd.read_sql_query( query, mysql_con )

    data_df = pd.read_csv( 'data/zemu_1.2-60000_rscript_simplified-t14-id_30.csv.gz' )
    data_df = data_df.loc[ data_df['ScoreMethodID'] == 35000 ]

    data_df['abs_error'] = np.abs( data_df['total'] - data_df['ExperimentalDDG'] )
    data_df.sort_values( 'abs_error', inplace = True )
    data_ids = data_df[['DataSetID']]
    data_ids = data_ids.merge( new_df, left_on = 'DataSetID', right_on = 'UserPPDataSetExperimentID' )
    ids_in_order = list( data_ids['ID'] )

    first_to_run = []
    last_to_run = []
    grab_first = True
    while( len(ids_in_order) > 0 ):
        if grab_first:
            ppi_id = ids_in_order.pop(0)
        else:
            ppi_id = ids_in_order.pop()
        grab_first = not grab_first

        for struct_id in range(1, 31):
            first_to_run.append( '%d/%04d' % (ppi_id, struct_id) )

        for struct_id in range(31, 51):
            last_to_run.append( '%d/%04d' % (ppi_id, struct_id) )

    ids_to_run = first_to_run + last_to_run

    jd_path = '/home/kyleb/gits/ddg/job_output/180115-kyleb_zemu_1.2-60000_struct-t14/data-0/job_dict-0.pickle'
    with open(jd_path, 'r') as f:
        jobs = pickle.load(f)

    assert( len(jobs) == len(ids_to_run) )
    assert( sorted(jobs.keys()) == sorted(ids_to_run) )

    with open( os.path.join( os.path.dirname(jd_path), 'sorted_ids.pickle' ), 'wb' ) as f:
        pickle.dump( ids_to_run, f )

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
