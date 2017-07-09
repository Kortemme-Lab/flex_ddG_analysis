library(mgcv)
library(data.table)
library(car)

all_data = fread('zcat ~/data/ddg/interface_ddg_paper/control_and_60k-score_terms.csv.gz', header = T, sep = ',')

all_data = all_data[ ( ScoreMethodID == 40 | ScoreMethodID == 10000 | ScoreMethodID == 20000 | ScoreMethodID == 30000 | ScoreMethodID == 40000 | ScoreMethodID == 50000 | ScoreMethodID == 60000 ) ] # 10k intervals

### all_data = all_data[ (MutType == 'complete' | MutType == 's2l') & (ScoreMethodID == 2500 | ScoreMethodID == 5000)] # Shorter dataset for testing

all_data$total_diff = ( all_data$fa_atr + all_data$fa_dun + all_data$fa_elec + all_data$fa_intra_rep + all_data$fa_rep + all_data$fa_sol + all_data$hbond_bb_sc + all_data$hbond_lr_bb + all_data$hbond_sc ) - all_data$total

print( 'Dataframe loaded successfully\n' )
print( paste0( 'Number of rows where total score sum does not match within tolerance: ', sum( abs(all_data$total_diff) >= 0.01 ) ) )

output_dir = "output_R"
unlink(output_dir, recursive=TRUE)
dir.create(output_dir, showWarnings = FALSE)

## save column names of prediction columns for later use
pred_col_names = c()

calc_gam <- function(df, by_labels, img_type) {
    df_args <- c(by_labels, sep="-")
    unique_name = do.call(paste, df_args)
    print( unique_name )
    dir.create(file.path(output_dir, unique_name), showWarnings = FALSE)

    zz = file( file.path( file.path(output_dir, unique_name), 'gam_summary'), open = "wt")
    sink(zz)

    gamobj <- gam( ExperimentalDDG ~ s(fa_atr, fx=TRUE, k=-1, bs="cs") + s(fa_elec, fx=TRUE, k=-1, bs="cs") + s(fa_rep, fx=TRUE, k=-1, bs="cs") + s(fa_sol, fx=TRUE, k=-1, bs="cs") + s(hbond_bb_sc, fx=TRUE, k=-1, bs="cs") + s(hbond_lr_bb, fx=TRUE, k=-1, bs="cs") + s(hbond_sc, fx=TRUE, k=-1, bs="cs"),
                  ## family=gaussian(link=identity),
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

    ## Save lpmatrix for later fitting
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
    Xp <- predict(gamobj, newd, type="lpmatrix")

    df$gamTotal = predict(gamobj, df)

    ## Save lpmatrix, coeff, step, boundary; df with predictions
    out_path = file.path( file.path(output_dir, unique_name),
                         'lpmatrix.csv' )
    write.table( Xp, out_path )
    system( paste0("gzip ", out_path) )

    out_path = file.path( file.path(output_dir, unique_name),
                         'coeff' )
    write.table( coef(gamobj), out_path )

    out_path = file.path( file.path(output_dir, unique_name),
                         'predictions.csv' )
    write.table( df, out_path )
    system( paste0("gzip ", out_path) )

    out_path = file.path( file.path(output_dir, unique_name),
                         'gamobj.Rdata.gz' )
    save(
        gamobj,
        file = out_path,
        compress = "gzip"
        )

    zz = file( file.path( file.path(output_dir, unique_name), 'step-boundary'), open = "wt")
    sink(zz)
    print(step)
    print(boundary)
    sink()
    close(zz)

    ## Predict values on all_data
    pred_col_name = paste0( unique_name, "-gamTotal" )
    pred_col_names <- c( pred_col_names, pred_col_name )
    assign( 'pred_col_names', pred_col_names, envir=.GlobalEnv )
    all_data$new_pred = predict(gamobj, all_data)
    names(all_data)[names(all_data) == "new_pred"] <- pred_col_name
    assign( 'all_data', all_data, envir=.GlobalEnv )

    ## Check calculating prediction manually
    ## xn_vals <- c(3.7054, 0.1415, -0.446, -0.8832, -0.0129, 0.1304, 0.664) ## want prediction at these values
    ## xn <- (xn_vals + boundary) / ( 2 * boundary ) # Convert to percentage in boundary range

    ## ncols <- length(colnames(Xp))
    ## nsmoothterms <- length(xn)
    ## nsections <- ( ncols - 1 ) / nsmoothterms
    ## print( paste0( 'nsmoothterms ', nsmoothterms, ' nsections ', nsections ) )

    ## x0 <- 1         ## intercept column
    ## for (j in 0:(nsmoothterms-1)) { ## loop through smooth terms
    ##     cols <- 1+j*nsections + 1:nsections      ## relevant cols of Xp
    ##     i <- floor( xn[j+1] * nrow(Xp) )  ## find relevant rows of Xp
    ##     w1 <- xn_vals[j+1] %% step / step
    ##     ## find approx. predict matrix row portion, by interpolation
    ##     x0 <- c(x0,Xp[i+2,cols]*w1 + Xp[i+1,cols]*(1-w1))
    ## }
    ## dim(x0) <- c(1,ncols)
    ## fv <- x0 %*% coef(gamobj) ## evaluate
    ## ## compare to normal prediction
    ## print( predict(gamobj, newdata=data.frame(
    ##                            fa_atr=xn_vals[1],
    ##                            fa_elec=xn_vals[2],
    ##                            fa_rep=xn_vals[3],
    ##                            fa_sol=xn_vals[4],
    ##                            hbond_bb_sc=xn_vals[5],
    ##                            hbond_lr_bb=xn_vals[6],
    ##                            hbond_sc=xn_vals[7]
    ##                        ),
    ##                se=FALSE) )

    return( sqrt(gamsum$r.sq) )
}

corr_summary = all_data[,.( gamR=calc_gam(.SD, .BY, 'pdf'), R=cor(total, ExperimentalDDG) ), by=.(PredictionRunName,ScoreMethodID,MutType)] # .SD is subset for each group by, .BY is group labels

write.csv( corr_summary, file.path( output_dir, "corr_summary.csv") )

print( corr_summary )

print( "GAM on all backrub steps" )

steps_corr_summary = all_data[,.( gamR=calc_gam(.SD, .BY, 'png'), R=cor(total, ExperimentalDDG) ), by=.(PredictionRunName,MutType)] # .SD is subset for each group by, .BY is group labels
print( steps_corr_summary )

write.csv( steps_corr_summary, file.path( output_dir, "steps_corr_summary.csv") )

## Save all data
print( "Saving all data" )
out_path = file.path( output_dir, "all_data.csv")
write.csv( all_data, out_path )
system( paste0("gzip ", out_path) )

## Cross correlation of gam prediction columns
print( "Cross correlation" )
selection_cols = c( pred_col_names, c("ExperimentalDDG", "PredictionRunName", "MutType", "ScoreMethodID") )
complete_data = all_data[ , selection_cols, with=FALSE ]

alt_gam_r = NULL
for ( pred_col in pred_col_names ) {
    print( pred_col )
    print( head(complete_data[,pred_col,with=FALSE]) )
    print( "ExperimentalDDG" )
    print( head(complete_data[,"ExperimentalDDG",with=FALSE]) )
    print( cor(complete_data[,pred_col,with=FALSE], complete_data[,"ExperimentalDDG",with=FALSE]) )
    alt_gam_r_inner = complete_data[,.( R=cor(get(pred_col), ExperimentalDDG, use = "complete.obs", method = "pearson") ), by=.(PredictionRunName,MutType,ScoreMethodID)] # .SD is subset for each group by, .BY is group labels
    names(alt_gam_r_inner)[names(alt_gam_r_inner)=="R"] <- paste0(pred_col, "-expR")
    if( is.null(alt_gam_r) ) {
        alt_gam_r = alt_gam_r_inner
    } else {
        alt_gam_r = merge(alt_gam_r, alt_gam_r_inner, all=TRUE)
    }
}

print( alt_gam_r )
out_path = file.path( output_dir, "alt_gam_r.csv")
write.csv( alt_gam_r, out_path )
