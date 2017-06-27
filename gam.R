library(mgcv)
library(data.table)

all_data = fread('zcat ~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz', header = T, sep = ',')
# all_data = read.csv( gzfile('~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz','rt') )
# all_data = read.csv( '~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv' )

# all_data = all_data[ MutType == 'complete' & (ScoreMethodID == 2500 | ScoreMethodID == 5000)] # Shorter dataset for testing

all_data$total_diff = ( all_data$fa_atr + all_data$fa_dun + all_data$fa_elec + all_data$fa_intra_rep + all_data$fa_rep + all_data$fa_sol + all_data$hbond_bb_sc + all_data$hbond_lr_bb + all_data$hbond_sc ) - all_data$total

print( 'Dataframe loaded successfully\n' )
print( paste0( 'Number of rows where total score sum does not match within tolerance: ', sum( abs(all_data$total_diff) >= 0.01 ) ) )

output_dir = "output_R"
dir.create(output_dir, showWarnings = FALSE)

calc_gam <- function(df,by_labels) {
    unique_name = paste(by_labels$PredictionRunName, by_labels$MutType, by_labels$ScoreMethodID, sep="-")
    dir.create(file.path(output_dir, unique_name), showWarnings = FALSE)
    print( unique_name )

    zz = file( file.path( file.path(output_dir, unique_name), 'gam_summary.csv'), open = "wt")
    sink(zz)

    gamobj <- gam( ExperimentalDDG ~ s(fa_atr, fx=TRUE, k=-1, bs="cs") + s(fa_elec, fx=TRUE, k=-1, bs="cs") + s(fa_rep, fx=TRUE, k=-1, bs="cs") + s(fa_sol, fx=TRUE, k=-1, bs="cs") + s(hbond_bb_sc, fx=TRUE, k=-1, bs="cs") + s(hbond_lr_bb, fx=TRUE, k=-1, bs="cs") + s(hbond_sc, fx=TRUE, k=-1, bs="cs"),
                  # family=gaussian(link=identity),
                  data = df,
                  control = gam.control(nthreads = 4)
                  )

    gamsum = summary(gamobj)
    print(gamsum)
    pdf( file.path( file.path(output_dir, unique_name), sprintf('gam_results-%s.pdf', unique_name)) )
    par(mfrow=c(3,3))
    plot(
        gamobj,
        scheme=1
    )
    frame()
    legend("topleft", legend=" ",
           title=sprintf("GAM R: %.2f", sqrt(gamsum$r.sq) ),
           bty='n')
    dev.off()
    sink()
    close(zz)
    return( sqrt(gamsum$r.sq) )
}

corr_summary = all_data[,.( gamR=calc_gam(.SD, .BY), R=cor(total, ExperimentalDDG) ), by=.(PredictionRunName,ScoreMethodID,MutType)] # .SD is subset for each group by, .BY is group labels

write.csv( corr_summary, file.path( output_dir, "corr_summary") )

print( corr_summary )
