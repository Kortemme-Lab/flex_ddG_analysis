import pandas as pd
import numpy as np
import os
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

total_score_abs_cutoff = 1.0
exp_score_abs_cutoff = 1.0
error_abs_cutoff = 2.0

for fpath in ['control_GAM_terms.csv', 'tal_GAM_terms.csv', 'ref_GAM_terms.csv']: #sorted([x for x in os.listdir(os.getcwd()) if x.endswith('_terms.csv')]):
    df = pd.read_csv(fpath)
    print fpath

    reg_columns = [x for x in df.columns if (not x.endswith('_GAM') and x != 'exp_data')]
    gam_columns = [x for x in df.columns if x.endswith('_GAM')]
    df['total'] = df[reg_columns].sum( axis = 1 )
    df['total_GAM'] = df[gam_columns].sum( axis = 1 )
    df['error'] = df['total'] - df['exp_data']
    df['error_GAM'] = df['total_GAM'] - df['exp_data']
    print 'MAE: %.3f' % np.mean(np.abs(df['error']))
    print 'MAE GAM: %.3f' % np.mean(np.abs(df['error_GAM']))
    print df[['total', 'total_GAM', 'exp_data']].corr()
    df = df.loc[ (np.abs(df['total']) >= total_score_abs_cutoff) & (np.abs(df['exp_data']) >= exp_score_abs_cutoff) & (np.abs(df['error']) >= error_abs_cutoff) ]
    print 'MAE (outliers): %.3f' % np.mean(np.abs(df['error']))
    print 'MAE (outliers) GAM: %.3f' % np.mean(np.abs(df['error_GAM']))
    print df[['total', 'total_GAM', 'exp_data']].corr()
    print

    new_rows = []
    for index, row in df.iterrows():
        for col_name, value in row.iteritems():
            if col_name == 'exp_data':
                prediction_type = 'exp_data'
                score_term_name = 'total'
            elif col_name.endswith('_GAM'):
                prediction_type = 'GAM'
                score_term_name = col_name[:-4]
            else:
                prediction_type = 'Rosetta'
                score_term_name = col_name
            new_rows.append( (score_term_name, prediction_type, value) )

    categorical_df =  pd.DataFrame.from_records(data = new_rows, columns = ['score_term', 'prediction_type', 'score'] )
    print categorical_df.head()
    fig = plt.figure(figsize=(16, 8.5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(
        'abs(total) >= %.1f ; abs(exp) >= %.1f ; abs(error) >= %.1f; GAM MAE %.2f ; Rosetta MAE %.2f' % (
            total_score_abs_cutoff, exp_score_abs_cutoff, error_abs_cutoff,
            np.mean(np.abs(df['error_GAM'])),
            np.mean(np.abs(df['error'])),
        )
    )
    ax.yaxis.grid(True)
    sns.swarmplot(
        x = "score_term", y = "score", hue = "prediction_type",
        data = categorical_df,
        # jitter = 0.2,
        # palette = "Set2",
        dodge = True,
        ax = ax,
        size = 3.5,
        alpha = 0.95,
    )
    fig.savefig( fpath + '-mpl.png' )

    continue

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
    print means
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
