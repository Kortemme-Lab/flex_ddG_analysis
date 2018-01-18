#!/usr/bin/python3

import sys
import os
import sqlite3
import shutil
import tempfile
from pprint import pprint
import pandas as pd
import numpy as np
import re
import datetime
import sys
import collections
import threading

from zip_util import unzip_file

rosetta_output_file_name = 'rosetta.out'
output_database_name = 'ddG.db3'
struct_output_database_name = 'struct.db3'
trajectory_stride = 15000
script_output_folder = 'output'

zemu_gam_params = {
    'fa_sol' :      (6.940, -6.722),
    'hbond_sc' :    (1.902, -1.999),
    'hbond_bb_sc' : (0.063,  0.452),
    'fa_rep' :      (1.659, -0.836),
    'fa_elec' :     (0.697, -0.122),
    'hbond_lr_bb' : (2.738, -1.179),
    'fa_atr' :      (2.313, -1.649),
}

resfile_mutations_cache = {}

# The Reporter class is useful for printing output for tasks which will take a long time
# It can even predict a finish time!


class Reporter:
    # Time in seconds function
    # Converts datetime timedelta object to number of seconds
    @staticmethod
    def ts(td):
        return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6

    @staticmethod
    def mean(l):
        # Not using numpy mean to avoid dependency
        return float( sum(l) ) / float( len(l) )

    def __init__( self, task, entries = 'files', print_output = True, eol_char = '\r', total_count = None ):
        self._lock = threading.Lock()
        self.print_output = print_output
        self.start = datetime.datetime.now()
        self.entries = entries
        self.lastreport = self.start
        self.task = task
        self.report_interval = datetime.timedelta( seconds = 1 ) # Interval to print progress
        self.n = 0
        self.completion_time = None
        if self.print_output:
            print('\nStarting ' + task)
        self.total_count = None # Total tasks to be processed
        self.maximum_output_string_length = 0
        self.rolling_est_total_time = collections.deque( maxlen = 50 )
        self.kv_callback_results = {}
        self.list_results = []
        self.eol_char = eol_char

        if total_count != None:
            self.set_total_count( total_count )

    def set_total_count(self, x):
        self.total_count = x
        self.rolling_est_total_time = collections.deque( maxlen = max(1, int( .05 * x )) )

    def decrement_total_count(self):
        if self.total_count:
            self.total_count -= 1

    def report(self, n):
        with self._lock:
            self.n = n
            time_now = datetime.datetime.now()
            if self.print_output and self.lastreport < (time_now - self.report_interval):
                self.lastreport = time_now
                if self.total_count:
                    percent_done = float(self.n) / float(self.total_count)
                    est_total_time_seconds = self.ts(time_now - self.start) * (1.0 / percent_done)
                    self.rolling_est_total_time.append( est_total_time_seconds )
                    est_total_time = datetime.timedelta( seconds = self.mean(self.rolling_est_total_time) )
                    time_remaining = est_total_time - (time_now - self.start)
                    eta = time_now + time_remaining
                    time_remaining_str = 'ETA: %s Est. time remaining: ' % eta.strftime("%Y-%m-%d %H:%M:%S")

                    time_remaining_str += str( datetime.timedelta( seconds = int(self.ts(time_remaining)) ) )

                    output_string = "  Processed: %d %s (%.1f%%) %s" % (n, self.entries, percent_done*100.0, time_remaining_str)
                else:
                    output_string = "  Processed: %d %s" % (n, self.entries)

                output_string += self.eol_char

                if len(output_string) > self.maximum_output_string_length:
                    self.maximum_output_string_length = len(output_string)
                elif len(output_string) < self.maximum_output_string_length:
                    output_string = output_string.ljust(self.maximum_output_string_length)
                sys.stdout.write( output_string )
                sys.stdout.flush()

    def increment_report(self):
        self.report(self.n + 1)

    def decrement_report(self):
        self.report(self.n - 1)

    def add_to_report(self, x):
        self.report(self.n + x)

    def done(self):
        self.completion_time = datetime.datetime.now()
        if self.print_output:
            print('Done %s, processed %d %s, took %s\n' % (self.task, self.n, self.entries, self.completion_time-self.start))

    def elapsed_time(self):
        if self.completion_time:
            return self.completion_time - self.start
        else:
            return time.time() - self.start

def gam_function(x, score_term = None ):
    return -1.0 * np.exp( zemu_gam_params[score_term][0] ) + 2.0 * np.exp( zemu_gam_params[score_term][0] ) / ( 1.0 + np.exp( -1.0 * x * np.exp( zemu_gam_params[score_term][1] ) ) )

def apply_zemu_gam(scores):
    new_columns = list(scores.columns)
    new_columns.remove('total_score')
    scores = scores.copy()[ new_columns ]
    for score_term in zemu_gam_params:
        assert( score_term in scores.columns )
        scores[score_term] = scores[score_term].apply( gam_function, score_term = score_term )
    scores[ 'total_score' ] = scores[ list(zemu_gam_params.keys()) ].sum( axis = 1 )
    scores[ 'score_function_name' ] = scores[ 'score_function_name' ] + '-gam'
    return scores

def rosetta_output_succeeded( potential_struct_dir ):
    path_to_rosetta_output = os.path.join( potential_struct_dir, rosetta_output_file_name )
    if not os.path.isfile(path_to_rosetta_output) and os.path.isfile(path_to_rosetta_output + '.gz'):
        unzip_file( path_to_rosetta_output + '.gz' )
    if not os.path.isfile(path_to_rosetta_output):
        return False

    db3_file = os.path.join( potential_struct_dir, output_database_name )
    if not os.path.isfile( db3_file ) and os.path.isfile( db3_file + '.gz' ):
        unzip_file( db3_file + '.gz' )
    if not os.path.isfile( db3_file ):
        return False

    struct_db3_file = os.path.join( potential_struct_dir, struct_output_database_name )
    if not os.path.isfile( struct_db3_file ) and os.path.isfile( struct_db3_file + '.gz' ):
        unzip_file( struct_db3_file + '.gz' )
    if not os.path.isfile( struct_db3_file ):
        return False

    success_line_found = False
    no_more_batches_line_found = False
    with open( path_to_rosetta_output, 'r' ) as f:
        for line in f:
            if line.startswith( 'protocols.jd2.JobDistributor' ) and 'reported success in' in line:
                success_line_found = True
            if line.startswith( 'protocols.jd2.JobDistributor' ) and 'no more batches to process' in line:
                no_more_batches_line_found = True

    return no_more_batches_line_found and success_line_found

def find_finished_jobs( output_folder ):
    return_dict = {}
    job_dirs = [ os.path.abspath(os.path.join(output_folder, d)) for d in os.listdir(output_folder) if os.path.isdir( os.path.join(output_folder, d) )]
    for job_dir in job_dirs:
        completed_struct_dirs = []
        for potential_struct_dir in sorted([ os.path.abspath(os.path.join(job_dir, d)) for d in os.listdir(job_dir) if os.path.isdir( os.path.join(job_dir, d) )]):
            if rosetta_output_succeeded( potential_struct_dir ):
                completed_struct_dirs.append( potential_struct_dir )
        if len(completed_struct_dirs) > 0:
            return_dict[job_dir] = completed_struct_dirs

    return return_dict

def get_scores_from_db3_file(db3_file, struct_number, case_name):
    if db3_file.endswith('.gz'):
        tmp_dir = tempfile.mkdtemp(prefix='unzip_db3_')
        new_db3_path = os.path.join(tmp_dir, os.path.basename(db3_file))
        shutil.copy(db3_file, new_db3_path)
        db3_file = unzip_file(new_db3_path)
    else:
        tmp_dir = None

    conn = sqlite3.connect(db3_file)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    num_batches = c.execute('SELECT max(batch_id) from batches').fetchone()[0]

    scores = pd.read_sql_query('''
    SELECT batches.name, structure_scores.struct_id, score_types.score_type_name, structure_scores.score_value, score_function_method_options.score_function_name from structure_scores
    INNER JOIN batches ON batches.batch_id=structure_scores.batch_id
    INNER JOIN score_function_method_options ON score_function_method_options.batch_id=batches.batch_id
    INNER JOIN score_types ON score_types.batch_id=structure_scores.batch_id AND score_types.score_type_id=structure_scores.score_type_id
    ''', conn)

    def renumber_struct_id( struct_id ):
        return trajectory_stride * ( int(struct_id-1) // num_batches )

    scores['struct_id'] = scores['struct_id'].apply( renumber_struct_id )
    scores['name'] = scores['name'].apply( lambda x: x[:-9] if x.endswith('_dbreport') else x )
    scores = scores.pivot_table( index = ['name', 'struct_id', 'score_function_name'], columns = 'score_type_name', values = 'score_value' ).reset_index()
    scores.rename( columns = {
        'name' : 'state',
        'struct_id' : 'backrub_steps',
    }, inplace=True)
    scores['struct_num'] = struct_number
    scores['case_name'] = case_name

    conn.close()

    if tmp_dir:
        shutil.rmtree(tmp_dir)

    return scores

def process_finished_struct( output_path, case_name ):
    db3_file = os.path.join( output_path, output_database_name )
    assert( os.path.isfile( db3_file ) )
    struct_number = int( os.path.basename(output_path) )
    scores_df = get_scores_from_db3_file( db3_file, struct_number, case_name )

    struct_db3_file = os.path.join( output_path, struct_output_database_name )
    assert( os.path.isfile( struct_db3_file ) )
    struct_number = int( os.path.basename(output_path) )
    struct_df = get_struct_info_from_db3_file( struct_db3_file, struct_number, case_name )

    return (scores_df, struct_df)

def calc_ddg( scores ):
    total_structs = np.max( scores['struct_num'] )

    nstructs_to_analyze = set([1, total_structs])
    for x in range(1, total_structs):
        if x % 10 == 0:
            nstructs_to_analyze.add(x)
    nstructs_to_analyze = sorted(nstructs_to_analyze)

    all_ddg_scores = []
    for nstructs in nstructs_to_analyze:
        ddg_scores = scores.loc[ ((scores['state'] == 'unbound_mut') | (scores['state'] == 'bound_wt')) & (scores['struct_num'] <= nstructs) ].copy()
        for column in ddg_scores.columns:
            if column not in ['state', 'case_name', 'backrub_steps', 'struct_num', 'score_function_name']:
                ddg_scores.loc[:,column] *= -1.0
        ddg_scores = ddg_scores.append( scores.loc[ ((scores['state'] == 'unbound_wt') | (scores['state'] == 'bound_mut')) & (scores['struct_num'] <= nstructs) ].copy() )
        ddg_scores = ddg_scores.groupby( ['case_name', 'backrub_steps', 'struct_num', 'score_function_name'] ).sum().reset_index()

        if nstructs == total_structs:
            struct_scores = ddg_scores.copy()

        ddg_scores = ddg_scores.groupby( ['case_name', 'backrub_steps', 'score_function_name'] ).mean().round(decimals=5).reset_index()
        new_columns = list(ddg_scores.columns.values)
        new_columns.remove( 'struct_num' )
        ddg_scores = ddg_scores[new_columns]
        ddg_scores[ 'scored_state' ] = 'ddG'
        ddg_scores[ 'nstruct' ] = nstructs
        all_ddg_scores.append(ddg_scores)

    return (pd.concat(all_ddg_scores), struct_scores)

def calc_dgs( scores ):
    l = []

    total_structs = np.max( scores['struct_num'] )

    nstructs_to_analyze = set([1, total_structs])
    for x in range(1, total_structs):
        if x % 10 == 0:
            nstructs_to_analyze.add(x)
    nstructs_to_analyze = sorted(nstructs_to_analyze)

    for state in ['mut', 'wt']:
        for nstructs in nstructs_to_analyze:
            dg_scores = scores.loc[ (scores['state'].str.endswith(state)) & (scores['state'].str.startswith('unbound')) & (scores['struct_num'] <= nstructs) ].copy()
            for column in dg_scores.columns:
                if column not in ['state', 'case_name', 'backrub_steps', 'struct_num', 'score_function_name']:
                    dg_scores.loc[:,column] *= -1.0
            dg_scores = dg_scores.append( scores.loc[ (scores['state'].str.endswith(state)) & (scores['state'].str.startswith('bound')) & (scores['struct_num'] <= nstructs) ].copy() )
            dg_scores = dg_scores.groupby( ['case_name', 'backrub_steps', 'struct_num', 'score_function_name'] ).sum().reset_index()
            dg_scores = dg_scores.groupby( ['case_name', 'backrub_steps', 'score_function_name'] ).mean().round(decimals=5).reset_index()
            new_columns = list(dg_scores.columns.values)
            new_columns.remove( 'struct_num' )
            dg_scores = dg_scores[new_columns]
            dg_scores[ 'scored_state' ] = state + '_dG'
            dg_scores[ 'nstruct' ] = nstructs
            l.append( dg_scores )
    return l

def mutations_from_resfile( resfile_path ):
    abs_resfile_path = os.path.abspath( resfile_path )
    if abs_resfile_path in resfile_mutations_cache:
        return resfile_mutations_cache[abs_resfile_path]

    mutations = set()
    with open(abs_resfile_path, 'r') as f:
        for line in f:
            m = re.match( r'(\d+)([A-Z]?)(?: )([A-Z])(?: .*?)', line )
            if m:
                resnum = int( m.group(1) )
                insertion_code = m.group(2)
                if len(insertion_code) == 0:
                    insertion_code = ' '
                chain = m.group(3)

                mutations.add( (resnum, insertion_code, chain) )

    resfile_mutations_cache[abs_resfile_path] = sorted(mutations)
    return( resfile_mutations_cache[abs_resfile_path] )

def get_struct_info_from_db3_file( struct_db3_file, struct_number, case_name ):
    structures_per_stride = [ 'backrub', 'wt_bound', 'mut_bound' ]

    task_id = os.path.basename( os.path.dirname( os.path.dirname(stuct_db3_file) ) )
    data_dir = os.path.join( os.path.join( os.path.dirname( os.path.dirname( os.path.dirname(stuct_db3_file) ) ), 'data' ), task_id )
    resfile_path = os.path.join( data_dir, 'mutations.resfile' )
    mutations = mutations_from_resfile( os.path.join( os.path.dirname(struct_db3_file), resfile_path ) )

    conn = sqlite3.connect(struct_db3_file)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    select_query = '''
    SELECT
    residue_pdb_identification.struct_id, residue_pdb_identification.pdb_residue_number,
    residue_pdb_identification.chain_id, residue_pdb_identification.insertion_code,
    residue_pdb_identification.residue_number,
	residues.name3,
	residue_rotamers.rotamer_bin_probability,
	residue_burial.sasa_r140,
    phi, psi, omega, chi1, chi2, chi3, chi4
    FROM protein_residue_conformation

    LEFT JOIN residue_pdb_identification ON residue_pdb_identification.struct_id=protein_residue_conformation.struct_id
	AND residue_pdb_identification.residue_number=protein_residue_conformation.seqpos

	LEFT JOIN residue_rotamers ON residue_rotamers.struct_id=protein_residue_conformation.struct_id
	AND residue_rotamers.residue_number=protein_residue_conformation.seqpos

	LEFT JOIN residues on residues.resNum=protein_residue_conformation.seqpos
	AND residues.struct_id=protein_residue_conformation.struct_id

	LEFT JOIN residue_burial ON residue_burial.struct_id=protein_residue_conformation.struct_id
	AND residue_burial.resNum=protein_residue_conformation.seqpos

    WHERE
    '''

    where_query = '\nOR\n'.join( ["(pdb_residue_number = %d AND chain_id = '%s' AND insertion_code = '%s')" % (mutation[0], mutation[2], mutation[1]) for mutation in mutations] )

    select_query += where_query

    df = pd.read_sql_query(select_query, conn)
    df['chi1_chi2'] = df['chi1'] + df['chi2']

    def state_apply( struct_id ):
        if struct_id == 1:
            return 'input_wildtype'
        elif struct_id == 2:
            return 'input_minimized'
        else:
            return structures_per_stride[(struct_id-3) % len(structures_per_stride)]

    df['state'] = df['struct_id'].apply( state_apply )

    def renumber_struct_id( struct_id ):
        if struct_id == 1:
            return -2
        elif struct_id == 2:
            return -1
        else:
            return trajectory_stride * (( struct_id - 2 - 1 ) // len(structures_per_stride))

    df['struct_id'] = df['struct_id'].apply( renumber_struct_id )

    df.rename( columns = {
        'struct_id' : 'backrub_steps',
    }, inplace=True)
    df['struct_num'] = struct_number
    df['case_name'] = case_name
    sort_columns = ['case_name', 'struct_num', 'backrub_steps', 'pdb_residue_number', 'insertion_code', 'state']
    df.sort_values( by = sort_columns, inplace = True )
    other_columns = [ x for x in df.columns.values if x not in sort_columns ]
    df = df[ sort_columns + other_columns ]

    conn.close()

    return df

def analyze_output_folder( output_folder ):
    # Pass in an outer output folder. Subdirectories are considered different mutation cases, with subdirectories of different structures.
    finished_jobs = find_finished_jobs( output_folder )
    if len(finished_jobs) == 0:
        print( 'No finished jobs found' )
        return

    ddg_scores_dfs = []
    struct_scores_dfs = []
    structs_dfs = []
    r = Reporter( 'analyzing finished jobs', entries = 'jobs', total_count = len(finished_jobs) )
    for finished_job, finished_structs in finished_jobs.items():
        inner_scores_list = []
        inner_structs_list = []
        for finished_struct in finished_structs:
            inner_scores, inner_structs = process_finished_struct( finished_struct, os.path.basename(finished_job) )
            inner_scores_list.append( inner_scores )
            inner_structs_list.append( inner_structs )
        scores = pd.concat( inner_scores_list )
        structs = pd.concat( inner_structs_list )
        structs_dfs.append( structs )
        ddg_scores, struct_scores = calc_ddg( scores )
        struct_scores_dfs.append( struct_scores )
        ddg_scores_dfs.append( ddg_scores )
        ddg_scores_dfs.append( apply_zemu_gam(ddg_scores) )
        ddg_scores_dfs.extend( calc_dgs( scores ) )
        r.increment_report()
    r.done()

    if not os.path.isdir(script_output_folder):
        os.makedirs(script_output_folder)
    basename = os.path.basename(output_folder)

    pd.concat( struct_scores_dfs ).to_csv( os.path.join(script_output_folder, basename + '-struct_scores_results.csv' ) )

    pd.concat( structs_dfs ).to_csv( os.path.join(script_output_folder, basename + '-structs.csv' ) )

    df = pd.concat( ddg_scores_dfs )
    df.to_csv( os.path.join(script_output_folder, basename + '-results.csv') )
    pivot_df = df.copy( ['backrub_steps', 'case_name', 'scored_state', 'score_function_name', 'nstruct', 'total_score' ] )
    pivot_df = pivot_df.pivot_table( columns = ['case_name', 'backrub_steps', 'score_function_name', 'scored_state', 'nstruct'], values = 'total_score' ).reset_index()
    pivot_df.rename( columns = { 0 : 'total_score' }, inplace = True )
    pivot_df.to_csv( os.path.join(script_output_folder, basename + '-pivot_results.csv') )
    print( pivot_df.head( n = 20 ) )

if __name__ == '__main__':
    for folder_to_analyze in sys.argv[1:]:
        if os.path.isdir( folder_to_analyze ):
            analyze_output_folder( folder_to_analyze )
