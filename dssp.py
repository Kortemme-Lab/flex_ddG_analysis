import os
import sys
import pandas as pd
import re
import numpy as np

from klab.bio.dssp import ComplexDSSP

zemu_dataset_path = 'table-zemu-filtered.csv'

def main():
    zemu = pd.read_csv(zemu_dataset_path, dtype = {'PDBFileID' : str} )
    dssp_results = []
    for pdb_id in zemu['PDBFileID'].drop_duplicates():
        d = ComplexDSSP.from_RCSB( pdb_id )
        # raw_mutations = [ x.split(';') for x in zemu.loc[ zemu['PDBFileID'] == pdb_id ]['Mutations'] ]
        for index, row in zemu.loc[ zemu['PDBFileID'] == pdb_id ].iterrows():
            data_set_id = row['DataSetID']
            mutations = []
            for mut in row['Mutations'].split(';'):
                m = re.match( '([A-Z])\s+[A-Z](\d+)([a-zA-Z]??)[A-Z]', mut )
                if not m:
                    print mut
                    assert( m )
                chain = m.group(1)
                resi = m.group(2)
                insertion_code = m.group(3)
                if len(insertion_code) == 0:
                    resi = resi + ' '
                else:
                    resi = resi + insertion_code.upper()
                mutations.append( (chain, resi) )

            ss = set()
            exposures = []
            for chain, resi in mutations:
                # print d.dssp[chain].keys()
                # print d.dssp[chain][ ('%s' % resi).rjust(5) ]
                r = d.dssp[chain][ ('%s' % resi).rjust(5) ]
                exposures.append( r['exposure'] )
                if len(r['ss'].strip()) > 0:
                    ss.add( r['ss'].strip() )
            if len(ss) == 1:
                ss = 'all-' + ss.pop()
            else:
                ss = None

            exposures = np.array(exposures)
            dssp_results.append( (
                data_set_id, ss,
                np.min(exposures), np.median(exposures), np.mean(exposures), np.max(exposures),
            ) )
    dssp_results = pd.DataFrame( dssp_results, columns = ['DataSetID', 'SS', 'min_exposure', 'median_exposure', 'mean_exposure', 'max_exposure'] )
    dssp_results.to_csv( 'dssp_results.csv' )

    print dssp_results.head()

if __name__ == '__main__':
    main()
