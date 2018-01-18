#!/usr/bin/python3

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


# Configure Seaborn
sns.set_style("whitegrid")

top_x = 21

def torsions( run_name, ids_to_use, sub_name ):
    output_dir = os.path.join( os.path.join( '/dbscratch/kyleb/local_data/output', run_name ), sub_name )
    if not os.path.isdir( output_dir ):
        os.makedirs( output_dir )

    input_df_path = os.path.join('/dbscratch/kyleb/local_data/output', run_name + '-structs.csv')

    assert( len(ids_to_use) == top_x )

    input_df = pd.read_csv( input_df_path )

    input_df = input_df.loc[ input_df['case_name'].isin(ids_to_use) ]
    if not len(input_df['case_name'].drop_duplicates()) == top_x:
        print( 'WARNING', len(input_df['case_name'].drop_duplicates()), top_x )
    # assert( len(input_df['case_name'].drop_duplicates()) == top_x )

    ##########################################################################################

    df = input_df.copy()[[
        'case_name', 'state', 'backrub_steps',
        'rotamer_bin_probability',
    ]].copy()
    df = df.loc[ df['backrub_steps'] >= 0 ]

    g = sns.factorplot(
        x = "backrub_steps", y = "rotamer_bin_probability",
        # hue = "state",
        col = "state",
        data = df, kind = "box",
        size = 4, aspect = .7
    )

    # g = sns.FacetGrid( df, col = "backrub_steps", hue = 'state', margin_titles = True )
    # g = g.map(sns.boxplot, "rotamer_bin_probability")
    fig = g.fig
    fig.savefig( os.path.join( output_dir, run_name + '-rotamer_bin.svg' ) )

    ##########################################################################################

    df = input_df.copy()[[
        'case_name', 'state', 'backrub_steps', 'struct_num', 'residue_number',
        'phi', 'psi', 'chi1', 'chi2',
    ]]
    df = df.loc[ df['backrub_steps'] >= -1 ]

    for col_name in ['phi', 'psi', 'chi1', 'chi2']:
        df[col_name] = df[col_name].apply( lambda x : None if np.abs(x) < 0.0000001 else x )

    dtor = np.pi / 180.0
    rtod = 180.0 / np.pi

    df = df.merge(
        df.loc[ df['backrub_steps'] == -1 ][ ['case_name', 'struct_num', 'residue_number', 'state'] + ['phi', 'psi', 'chi1', 'chi2'] ],
        how = 'left',        on = ['case_name', 'struct_num', 'residue_number'],
    )

    for col_name in ['phi', 'psi', 'chi1', 'chi2']:
        df[col_name] = rtod * np.arcsin( np.sin( dtor * df[col_name + '_x'] - dtor * df[col_name + '_y'] ) )

    df.rename( columns = { 'state_x' : 'state' }, inplace = True )
    df = df[ ['case_name', 'struct_num', 'residue_number', 'state', 'backrub_steps'] + ['phi', 'psi', 'chi1', 'chi2'] ]

    df['chi1_chi2'] = df['chi1'] + df['chi2']

    df = df.loc[ df['backrub_steps'] > 0 ]

    melt_df = pd.melt(
        df,
        id_vars = ['case_name', 'state', 'backrub_steps', 'struct_num'],
        value_vars = ['phi', 'psi', 'chi1', 'chi2', 'chi1_chi2'],
        value_name = 'torsion_delta',
        var_name = 'torsion_type',
    )

    g = sns.FacetGrid( melt_df, col = "backrub_steps", row = "torsion_type", margin_titles = True )
    g = g.map(plt.hist, "torsion_delta")
    fig = g.fig
    fig.savefig( os.path.join( output_dir, run_name + 'facet_grid.pdf' ) )

    # Facet grid 2
    g = sns.factorplot(
        x = "backrub_steps", y = "torsion_delta",
        row = "state",
        col = "torsion_type",
        data = melt_df, kind = "box",
        margin_titles = True,
        size = 4, aspect = .7
    )

    # g = sns.FacetGrid( df, col = "backrub_steps", hue = 'state', margin_titles = True )
    # g = g.map(sns.boxplot, "rotamer_bin_probability")
    fig = g.fig
    fig.savefig( os.path.join( output_dir, run_name + 'facet_grid2.pdf' ) )

    # Enumerate from facet grid 2
    for torsion_type in melt_df['torsion_type'].drop_duplicates():
        for state in melt_df['state'].drop_duplicates():
            melt_df_subset = melt_df.loc[ (melt_df['torsion_type'] == torsion_type) & (melt_df['state'] == state) ]
            fig_path = os.path.join( output_dir, '%s-%s.svg' % (state, torsion_type) )

            fig, ax = plt.subplots( figsize = (4.5,4.5), dpi = 200 )

            sns.boxplot( x = "backrub_steps", y = "torsion_delta", data = melt_df_subset, ax = ax )
            # sns.violinplot( x = "backrub_steps", y = "torsion_delta", data = melt_df_subset, ax = ax )
            if torsion_type in ['phi', 'psi']:
                ax.set_ylim( -110, 110 )
            else:
                ax.set_ylim( -160, 160 )
            ax.set_xlabel( 'Backrub steps' )
            ax.set_ylabel( '%s delta (from starting conformation)' % torsion_type.replace( 'chi1_chi2', 'chi1+chi2' ) )
            ax.set_title( sub_name )
            plt.tight_layout()
            fig.savefig( fig_path )

            # Histogram
            # sns.distplot(x, kde=False, rug=True)


def main(run_name):

    df = pd.read_csv( os.path.expanduser( '~/gits/interface_ddg/interesting.csv' ) )
    top_ids = df.iloc[:top_x]['ID']
    assert( len(top_ids) == top_x )
    bottom_ids = df.iloc[-top_x:]['ID']
    assert( len(bottom_ids) == top_x )

    torsions( run_name, top_ids, 'top%d' % top_x )
    torsions( run_name, bottom_ids, 'bottom%d' % top_x )

if __name__ == '__main__':
    main( '180115-kyleb_zemu_1.2-60000_struct-t14' )
