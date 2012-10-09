
# {{{1 Imports/boilerplate
library(grt)

library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

library(foreach)
library(doMC)
registerDoMC()
# }}}1

#
# {{{1 General statistics
bic.general <- function(negloglik, n, k) log( n )*k - 2*negloglik
# }}}1

# {{{1 Translation functions
radiens2degrees <- function( rad ) 180 * rad / pi

xy2angle <- function(xyvec) {
        # For a unit vector of x/y coordinates, finds their angle 
        # TODO: catch the case where xyvec is not a unit vector.
        if (xyvec[2] > 0) { 
                    theta.rad <- acos( xyvec[1] )  }
    else { 
                theta.rad <- 2*pi - acos( xyvec[1] ) }
        radiens2degrees( theta.rad )
}
# }}}1

# {{{1 GLC statistics
nobs.glc <- function(fitobj) {
    nrow(fitobj$model)
}

bic.model.glc <- function(fitobj ) {
    loglik <- logLik(fitobj)
    nparams <- attr(loglik, "df")
    this.nobs <- nrow(fitobj$model)
    bic.general(-c(loglik), this.nobs, nparams )
}

anova.glc <- function(lik0, lik1) {
    df <- attr( lik1, "df" ) - attr( lik0, "df" )
    if ( df < 1 ) stop( "First object must be the reduced model.")
    # if ( lik0 == 0 ) { return( 1 ) }
    # if ( lik1 == 0 ) { return( 0 ) }
    chisq <- -2 * (c(lik0) - c(lik1))
    # stopifnot ( lik0 <= lik1 )
    stopifnot ( chisq >= 0 )
    pval <- 1 - pchisq( chisq, df=df )
    pval
}
# }}}1

# {{{1 GRT model functions
#resps <- c("CH1", "CH2")
resps <- c(T,F)
labels <- c( "Unimodal", "Bimodal", "2D", "Null" ) 
oned.labels <- c( "Unimodal", "Bimodal", "Null" ) 
alpha <- .01

compare.models <- function(loglik0, loglik1) {
    # Returns true if the full model is better than the reduced model.
    if( attr(loglik0, 'df') > attr(loglik1, 'df') )  stop( "First object must be the reduced model.")
    if( attr(loglik0, 'df') == attr(loglik1, 'df') )  return (loglik0 < loglik1)
    if (loglik0 > loglik1 ) {
        warning( "Restricted model better fit than full model. Could be a numerical error." )
        return( F )
    }
    anova.glc( loglik0, loglik1 ) < alpha
}


null.model.loglik <- function( df ) {
    # code the responses
    n <- nrow(df)
    # Assuming ch1 == False
    nch1 <- sum(!df$resp)
    bias <- nch1 / n
    loglLik <- log(bias)*nch1 + log(1-bias) * (n-nch1)
    attr( loglLik, "df" ) <- 1
    loglLik
}

fits.table <- function(df, withplot=F) {
    # Metadata:
    if (unique(df$alllab) == "true") unlab <- "sham"
    else unlab <- "unlab"
    nlab <- unique(df$nlab)
    order <- unique(df$order)
    
    labcond <- paste( unlab, nlab, order, sep='.' )
    
    if (length( unique ( df$resp ) ) < 2) {
        warning( paste("always gave the same response.") )
        return (
            data.frame(BestFit="Null", BestFitOne="Null", Used="Null",
                        unimodcoeff=0, bimodcoeff=0, twodAngle=-999,
                        twodBimod=0, twodNoise=0, twodBias=1,
                        labcond=labcond, unlab=unlab,
                        nlab=nlab, order=order )
            )
    }
    model.formulas <- list(Unimodal = resp ~ unimod,
                           Bimodal = resp ~ bimod,
                           "2D" = resp ~ bimod + unimod)
    fit.glc <- function(form) {
        model.fit <- glc( form, data=df )
        if ( model.fit$par$noise == 10 ) {
            df <- attr(model.fit$logLik, "df")
            model.fit$logLik <- -Inf 
            attr(model.fit$logLik, "df") <- df
        }
        model.fit
    }
    glc.fits <- llply( model.formulas, fit.glc )
    
    # Record oned coeffs # TODO may not be working
    unimodcoeff <- glc.fits$Unimodal$par$coeffs
    bimodcoeff <- glc.fits$Bimodal$par$coeffs
    
    # Record the twod coeffs
    twodcoeffs <- glc.fits$"2D"$par$coeffs
    twod.angle <- xy2angle( twodcoeffs )
    twod.bimod <- twodcoeffs[1] 
    twod.noise <- glc.fits$"2D"$par$noise
    twod.bias <- glc.fits$"2D"$par$bias
    
    # If we had to pick a 1D rule, which is it?
    if ( compare.models( logLik(glc.fits$Bimodal), logLik(glc.fits$Unimodal)) ) {
        used.index <- 1
    } else { used.index <- 2 }
    used <- labels[[used.index]]
    
    # Null model
    nullloglik <- null.model.loglik( df )
    # null2loglik <- log(.5) * nrow(df)
    # attr( null2loglik, "df" ) <- 0
    
    # Choose between null or the better 1D model:
    if ( compare.models(nullloglik, logLik(glc.fits[[used]])) ) {
        onedwinner <- used
        logLik.onedwinner <- logLik(glc.fits[[used]])
        onedwinner.coeff <- c( unimodcoeff, bimodcoeff, 0 )[used.index]
    } else {
        onedwinner <- "Null"
        logLik.onedwinner <- nullloglik
        onedwinner.coeff <- 0
    }
    
    # Choose between that winner and the 2D model
    if ( compare.models( logLik.onedwinner, logLik(glc.fits$"2D")) ) {
        winner <- "2D"
        winner.coeff <- twodcoeffs
    } else {
        winner <- onedwinner
        winner.coeff <- onedwinner.coeff
    }
    
    data.frame(BestFit=winner, BestFitOne=onedwinner, Used=used,
               unimodcoeff=unimodcoeff, bimodcoeff=bimodcoeff,
               twodAngle=twod.angle, twodBimod=twod.bimod,
               twodNoise=twod.noise, twodBias=twod.bias, labcond=labcond,
               unlab=unlab, nlab=nlab, order=order )
}
# }}}1

# {{{1 Reading in sims
# {{{2 Params
# }}}2
#
# Softmax, chooses true or false given a fixed probably of false
softmax_binomial <- function(pfalse, tau) {
  exponentiated <- exp(c(pfalse, 1-pfalse) / tau)
  odds_false <- exponentiated[1] / sum(exponentiated)
  runif(1)>odds_false
}

run_anderson_once <- function(...) {
  output_columns <- c("type","bimod","unimod","label","clust")
  command <- paste("./testanderson", ...)
  read.csv(pipe(command), col.names=output_columns)
}
# }}}1

# {{{1 Simulation code
simulate_anderson <- function(cparam=1, nlab=16, order="interspersed", encoding="actual", tau=.05, plotting=F, echo=F) {
  resps <- run_anderson_once(cparam, nlab, order, encoding)
  inferences <- subset(resps, type=="INFER")
  inferences$resp <- apply(matrix(inferences$label), 1, function(p) softmax_binomial(p, tau))
  if (echo) { print( inferences) }
  if (plotting) {
    fit <- glc(resp ~ bimod + unimod, data=inferences )
    theplot <- ggplot(inferences)
    theplot <- theplot  + geom_point(aes(x=bimod, y=unimod, colour=label))
    intercept <- -fit$par$bias / fit$par$coeff[2]
    slope <- -fit$par$coeff[1]/fit$par$coeff[2] 
    #plot( unimod ~ bimod, data=subset(inferences, label==F), col='blue' )
    #points( unimod ~ bimod, data=subset(inferences, label==T), col='red' )
    #abline( a=intercept, b=slope )
    theplot <- theplot + geom_abline(intercept=intercept, slope=slope)
    print(theplot)
  }
  
  # Prepping for use in fits.table
  inferences$alllab <- "false"
  inferences$order <- order
  inferences$nlab <- nlab
  inferences$encoding <- encoding
  
  fits.table(inferences)
}

plot_anderson <- function(cparam=1, nlab=16, order="interspersed", encoding="actual", tau=.05, plotting=F, echo=F) {
  resps <- run_anderson_once(cparam, nlab, order, encoding)
  stims <- subset(resps, type=="STIM")
  stims <- subset(stims, bimod>0 & unimod>0 ) # Remove unlabeld
  if (length(unique(stims$label)) == 3) labels <- c("unlabeled", "ch1", "ch2")
  else if (length(unique(stims$label)) == 2) labels <- c("ch1", "ch2")
  else if (length(unique(stims$label)) == 1) labels <- c("unlabeled")
  stims$label <- factor(stims$label, labels=labels)
  ggplot(stims) + geom_point( aes(y=unimod, x=bimod, colour=label, size=2) ) + geom_path(aes(x=bimod, y=unimod, alpha=1), arrow=arrow(type="closed", length=unit(.25,"inches"), angle=15))
}

plot_anderson(order="interspersed") # Sanity check

run_sims <- function(runs, nreps) {
  #ntotal <- nrow(runs) * nreps
  #pb <- txtProgressBar(min=0, max=ntotal, style=3)
  #count <- 0
  
  all_sims <- ddply(runs, names(runs), function(df) {
            attach(df[1,]) # Yes yes I know attaching is a dirty dirty thing to do.
            count <- 0
            ret <- data.frame()
            while ( count < nreps ) {
              this_sim <- try(simulate_anderson(cparam, nlab, order, encoding, tau, plotting=F))
              if (class(this_sim) == "try-error") next # Sometimes the call fails.
              count <- count + 1
              ret <- rbind( ret, this_sim )
            }
            ret
            #count <- count + nreps
            #setTxtProgressBar(pb, count)
                      }, .parallel=T)
  
  write.csv(all_sims, file="sims.csv")
  all_sims
}

# {{{2 Try one run
simulate_anderson(1, -1, plotting=T, echo=T)
# }}}2
# }}}1

# {{{1 Read the sim results
#sims <- read.csv("sims.csv")
runs <- expand.grid(nlab=c(-1), 
                    order=c("interspersed", "labfirst", "lablast"),
                    cparam=c(1, 1.8, .7/.3),
                    tau=c(.05), 
                    encoding=c("actual"))
runs <- subset( runs, ! ((encoding=="guess" | encoding=="softguess" ) & order!="interspersed"))

nreps <- 300
sims <- run_sims(runs, nreps)
print(summary(sims$BestFit))

counts <- ddply( sims, .(nlab, cparam, tau, order, encoding), function(x) summary(x$BestFit) )
counts$ratio <- counts$Bimodal/counts$"2D"
counts$params <- do.call(paste, c(counts[c("cparam", "tau")], sep = ":"))

counts$twod <- counts$"2D"
#print(ggplot(counts) + geom_line(aes(x=nlab, y=twod, group=params, linetype=factor(tau), colour=factor(cparam))) + facet_wrap(~order))

counts.melted <- melt(counts,
                      id=c("nlab", "cparam", "tau", "order", "encoding"),
                      measure.vars=c("twod", "Bimodal", "Unimodal", "Null"),
                      variable.name="bestfit",
                      value.name="count")

counts.melted$groupingparam <- do.call(paste, c(counts.melted[c("bestfit", "encoding")], sep = ":"))

print(ggplot(counts.melted) +
      geom_line(aes(x=factor(nlab), y=count, group=groupingparam, linetype=encoding, colour=bestfit)) + 
      facet_grid(cparam~order, labeller=label_both) +
      labs(title=paste(nreps, "simulated runs of the Anderson Rational model")))
# }}}1

# vim: foldmethod=marker
