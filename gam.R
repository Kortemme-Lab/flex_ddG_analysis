library(mgcv)
library(data.table)

all_data = fread('zcat ~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz', header = T, sep = ',')
# all_data = read.csv( gzfile('~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz','rt') )
# all_data = read.csv( '~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv' )

all_data = all_data[ MutType == 'complete' & (ScoreMethodID == 2500 | ScoreMethodID == 5000)] # Shorter dataset for testing

all_data$total_diff = ( all_data$fa_atr + all_data$fa_dun + all_data$fa_elec + all_data$fa_intra_rep + all_data$fa_rep + all_data$fa_sol + all_data$hbond_bb_sc + all_data$hbond_lr_bb + all_data$hbond_sc ) - all_data$total

print( 'Dataframe loaded successfully\n' )
print( paste0( 'Number of rows where total score sum does not match within tolerance: ', sum( abs(all_data$total_diff) >= 0.01 ) ) )

output_dir = "output_R"
dir.create(output_dir, showWarnings = FALSE)

calc_gam <- function(df, by_labels, img_type) {
    df_args <- c(by_labels, sep="-")
    unique_name = do.call(paste, df_args)
    print( unique_name )
    dir.create(file.path(output_dir, unique_name), showWarnings = FALSE)

    zz = file( file.path( file.path(output_dir, unique_name), 'gam_summary'), open = "wt")
    sink(zz)

    gamobj <- gam( ExperimentalDDG ~ s(fa_atr, fx=TRUE, k=-1, bs="cs") + s(fa_elec, fx=TRUE, k=-1, bs="cs") + s(fa_rep, fx=TRUE, k=-1, bs="cs") + s(fa_sol, fx=TRUE, k=-1, bs="cs") + s(hbond_bb_sc, fx=TRUE, k=-1, bs="cs") + s(hbond_lr_bb, fx=TRUE, k=-1, bs="cs") + s(hbond_sc, fx=TRUE, k=-1, bs="cs"),
                  # family=gaussian(link=identity),
                  data = df,
                  control = gam.control(nthreads = 4)
                  )

    gamsum = summary(gamobj)
    print(gamsum)

    if( img_type == 'pdf' ) {
        pdf( file.path( file.path(output_dir, unique_name), sprintf('gam_results-%s.pdf', unique_name)) )
    } else {
        png(
            file.path( file.path(output_dir, unique_name),
                      sprintf('gam_results-%s.png', unique_name)),
            width=20, height=20, units="cm", res=400
        )
    }
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

    if( img_type == 'pdf' ) {
        pdf( file.path( file.path(output_dir, unique_name), sprintf('gam_results-residuals-%s.pdf', unique_name)) )
    } else {
        png(
            file.path( file.path(output_dir, unique_name),
                      sprintf('gam_results-residuals-%s.png', unique_name)),
            width=20, height=20, units="cm", res=400
        )
    }
    par(mfrow=c(3,3))
    plot(
        gamobj,
        scheme=1,
        residuals=TRUE
    )
    frame()
    legend("topleft", legend=" ",
           title=sprintf("GAM R: %.2f", sqrt(gamsum$r.sq) ),
           bty='n')
    dev.off()

    sink()
    close(zz)


    # Save lpmatrix for later fitting
    step = 0.01
    boundary = 40
    newd <- data.frame(
        fa_atr = ( (-boundary/step):(boundary/step) ) * step,
        fa_elec = ( (-boundary/step):(boundary/step) ) * step,
        fa_rep = ( (-boundary/step):(boundary/step) ) * step,
        fa_sol = ( (-boundary/step):(boundary/step) ) * step,
        hbond_bb_sc = ( (-boundary/step):(boundary/step) ) * step,
        hbond_lr_bb = ( (-boundary/step):(boundary/step) ) * step,
        hbond_sc = ( (-boundary/step):(boundary/step) ) * step
    )
    pred = predict.gam(gamobj, newd)
    print( "coef" )
    print( coef(gamobj) )
    Xp <- predict(gamobj, newd, type="lpmatrix")
    xn <- c(.341, .122, .476, .981, -2, -8, 4) ## want prediction at these values

    ncols <- length(colnames(Xp))
    nsmoothterms <- length(xn)


    x0 <- 1         ## intercept column
    dx <- step      ## covariate spacing in `newd'
    for (j in 0:2) { ## loop through smooth terms
        cols <- 1+j*9 +1:9      ## relevant cols of Xp
        i <- floor(xn[j+1]/step)  ## find relevant rows of Xp
        w1 <- (xn[j+1]-i*dx)/dx ## interpolation weights
        ## find approx. predict matrix row portion, by interpolation
        x0 <- c(x0,Xp[i+2,cols]*w1 + Xp[i+1,cols]*(1-w1))
    }
    dim(x0) <- c(1,ncols)
    fv <- x0 %*% coef(gamobj) + xn[4]    ## evaluate and add offset
    ## compare to normal prediction
    print( predict(b,newdata=data.frame(x0=xn[1],x1=xn[2],
                                        x2=xn[3],x3=xn[4]),se=FALSE) )
    print( fv )

    return( sqrt(gamsum$r.sq) )
}

corr_summary = all_data[,.( gamR=calc_gam(.SD, .BY, 'pdf'), R=cor(total, ExperimentalDDG) ), by=.(PredictionRunName,ScoreMethodID,MutType)] # .SD is subset for each group by, .BY is group labels

write.csv( corr_summary, file.path( output_dir, "corr_summary.csv") )

print( corr_summary )

print( "GAM on all backrub steps" )

steps_corr_summary = all_data[,.( gamR=calc_gam(.SD, .BY, 'png'), R=cor(total, ExperimentalDDG) ), by=.(PredictionRunName,MutType)] # .SD is subset for each group by, .BY is group labels
print( steps_corr_summary )
