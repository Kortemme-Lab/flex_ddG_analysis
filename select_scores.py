import pandas as pd
import argparse
import getpass
import MySQLdb
import os
import time

cache_dir = 'query_cache'
if not os.path.isdir( cache_dir ):
    os.makedirs( cache_dir )

def main( mysql_host, mysql_port, mysql_user, mysql_pass, mysql_db ):
    with open('scores_select.sql', 'r') as f:
        all_scores_query = ' '.join( [ line.strip() for line in f.readlines() if not line.startswith('#') ] )

    mysql_con = MySQLdb.connect( host = mysql_host,
                                 port = mysql_port,
                                 user = mysql_user, passwd = mysql_pass,
                                 db = mysql_db
    )

    print 'Loading all scores dataframe'
    all_scores_cache = os.path.join( cache_dir, 'all_scores.csv.gz' )
    if os.path.isfile( all_scores_cache ):
        all_scores = pd.read_csv( all_scores_cache )
    else:
        all_scores = pd.read_sql( all_scores_query, con = mysql_con )
        all_scores.to_csv( all_scores_cache, compression = 'gzip' )

    print 'Loaded all scores dataframe. Number rows:', len(all_scores)
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
    main( args.host, args.port, args.user, args.password, args.db)
