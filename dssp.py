import os
import sys
import pandas as pd
import re

from klab.bio.dssp import ComplexDSSP

zemu_dataset_path = 'table-zemu-filtered.csv'

def main():
    zemu = pd.read_csv(zemu_dataset_path)
    for pdb_id in zemu['PDBFileID'].drop_duplicates():
        dssp = ComplexDSSP.from_RCSB( pdb_id )
        raw_mutations = [ x.split(';') for x in zemu.loc[ zemu['PDBFileID'] == pdb_id ]['Mutations'] ]
        mutations = []

        for index, row in zemu.loc[ zemu['PDBFileID'] == pdb_id ].iterrows():
            mutations = []
            for mut in row['Mutations']:
                m = re.match( '([A-Z])\s+[A-Z](\d+)[A-Z]', mut )
                assert( m )
                chain = m.group(1)
                resi = int( m.group(2) )
                mutations.append( (chain, resi) )
            print mutations
        print dssp.keys()
        print zemu.head()
        sys.exit()

if __name__ == '__main__':
    main()
