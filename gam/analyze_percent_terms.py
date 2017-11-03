import pandas as pd
import numpy as np
import os

for fpath in ['control_GAM_terms.csv', 'tal_GAM_terms.csv', 'ref_GAM_terms.csv']: #sorted([x for x in os.listdir(os.getcwd()) if x.endswith('_terms.csv')]):
    df = pd.read_csv(fpath)
    print fpath

    reg_columns = [x for x in df.columns if not x.endswith('_GAM')]
    gam_columns = [x for x in df.columns if x.endswith('_GAM')]
    df['total'] = df[reg_columns].sum( axis = 1 )
    df['total_GAM'] = df[gam_columns].sum( axis = 1 )

    percent_df = pd.DataFrame()
    for col in [x for x in df.columns if (not x.startswith('total') and x != 'exp_data')]:
        if col.endswith('_GAM'):
            total_col = 'total_GAM'
        else:
            total_col = 'total'
        percent_df[col] = df[col] / df[total_col]
    # print df.head()
    means_percent_df = percent_df.mean()
    # print means_percent_df.head( n = 10 )

    means = np.array([value for index, value in sorted(means_percent_df.iteritems()) if not index.endswith('_GAM')])
    means_relative = means / np.max(np.abs(means))
    means_gam = [value for index, value in sorted(means_percent_df.iteritems()) if index.endswith('_GAM')]
    means_gam_relative = means_gam / np.max(np.abs(means_gam))
    relative_delta = means_gam_relative - means_relative
    columns = [index for index, value in sorted(means_percent_df.iteritems()) if not index.endswith('_GAM')]
    means_by_term_indices = ['mean_pct', 'mean_pct_gam', 'mean_pct_relative', 'mean_pct_gam_relative', 'relative_delta']
    means_by_term = pd.DataFrame.from_records([means, means_gam, means_relative, means_gam_relative, relative_delta], index =  means_by_term_indices, columns = columns )
    means_by_term_columns = sorted(means_by_term.columns)
    means_by_term['sum'] = means_by_term[means_by_term_columns].sum( axis = 1 )
    # means_by_term['absmax'] = means_by_term[means_by_term_columns].abs().max( axis = 1 )
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    print means_by_term.applymap("{0:.3f}".format)
    # print means_by_term.applymap(lambda x: x*100).applymap("{0:.0f}".format)

    print
    for index, row in means_by_term.loc[['mean_pct', 'mean_pct_gam']].iterrows():
        for col_name, value in row.sort_values( ascending = False ).iteritems():
            if col_name == 'sum':
                continue
            print index, col_name, '%.3f' % value
        print

    print '\n\n'
