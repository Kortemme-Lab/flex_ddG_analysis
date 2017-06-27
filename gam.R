library(mgcv)
library(data.table)

df = fread('zcat ~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz', header = T, sep = ',')
# df = read.csv( gzfile('~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz','rt') )
# df = read.csv( '~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv' )
df = df[ PredictionRunName == 'zemu_1.2-60000_rscript_validated-t14' & MutType == 'complete' & ScoreMethodID == 30000 ]

head( df, n=10 )

df$total_diff = ( df$fa_atr + df$fa_dun + df$fa_elec + df$fa_intra_rep + df$fa_rep + df$fa_sol + df$hbond_bb_sc + df$hbond_lr_bb + df$hbond_sc ) - df$total

print( 'Dataframe loaded successfully\n' )
print( paste0( 'Number of rows where total score sum does not match within tolerance: ', sum( abs(df$total_diff) >= 0.01 ) ) )

gamobj <- gam( ExperimentalDDG ~ s(fa_atr, fx=TRUE, k=7) + s(fa_elec, fx=TRUE, k=7) + s(fa_rep, fx=TRUE, k=7) + s(fa_sol, fx=TRUE, k=7) + s(hbond_bb_sc, fx=TRUE, k=7) + s(hbond_lr_bb, fx=TRUE, k=7) + s(hbond_sc, fx=TRUE, k=7),
              # family=gaussian(link=identity),
              data = df,
              control = gam.control(nthreads = 4)
              )

summary(gamobj)
pdf("gam_results.pdf")
par(mfrow=c(3,3))
plot(gamobj)
dev.off()
