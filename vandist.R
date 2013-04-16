
library(data.table)
#library(ggplot2)
source("simfunctions.R")

# {{{1 Functions
run_vandist_once <- function(...) {
  x <- run_anderson_once(task="vandist", ...)
  x$label <- x$bimod - x$unimod > 0
  x
}

# By default, sigma is .5 (pooled sd / 3)
#x <- run_anderson_once(task="vandist", a0=10, proplab=.5)
#x
#names(x)

plot_clusters <- function(...) {
  simresult <- run_anderson_once(task="vandist", ...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  print(stims)
  print(names(stims))
  ggplot(stims) + geom_point(aes(x=bimod, y=unimod, colour=factor(clust)))
}

plot_label <- function(...) {
  simresult <- run_vandist_once(...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  ggplot(stims) + geom_point(aes(x=bimod, y=unimod, colour=factor(label)))
}

plot_guesses <- function(...) {
  simresult <- run_anderson_once(task="vandist", ...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  print(stims)
  print(names(stims))
  ggplot(stims) + geom_point(aes(x=bimod, y=unimod, colour=guess))
}

run_vandist_sims <- function(runs, nreps) {
  do.run <- function(run) {
    tau <- run$tau
    bias_sd <- run$bias_sd
    run['bias_sd'] <- NULL
    run['tau'] <- NULL
    prop <- unique(run$prop)
    count <- 0
    ret <- data.frame()
    while (count<nreps) {
      if (bias_sd==0) { run$bias <- 0 }
      else run$bias <- rnorm(1, mean=0, sd=bias_sd)
      # Removing unusable items:
      output <- do.call(run_vandist_once, run)
      this_sim <- subset(output, type=="STIM")
      this_sim$proplab <- prop
      this_sim$trialnum <- 1:nrow(this_sim)
      this_sim$bias <- run$bias
      this_sim$hit <- apply(matrix(this_sim$guess), 1, function(p) softmax_binomial(p, tau))!=this_sim$label
      count <- count + 1
      ret <- rbind(ret, this_sim)
    }
    ret
  }
  ret <- ddply(runs, names(runs), do.run, .parallel=T)
  ret$block <- floor((ret$trialnum-1)/80)+1
  ret
}


vandist_et_al_sim <- function(bias_sd=0, n=1, tau=.05, ...) {
  do.run <- function(prop) {
    bias <- bias_sd
    output <- run_vandist_once(bias=bias, proplab=prop, ...)
    ret <- subset(output, type=="STIM")
    ret$proplab <- prop
    ret$trialnum <- 1:nrow(ret)
    ret
  }
  runs <- c(replicate(n, .5), replicate(n, 1))
  stims <- adply(runs, 1, do.run, .parallel=F)
  stims$block <- floor((stims$trialnum-1)/80)+1
  stims$hit <- softmax_binomial(stims$guess, tau)==stims$label
  ddply(stims, .(proplab,block), summarise, acc=mean(hit))
  
}
# }}}1

# {{{1 Test above functions
accuracies <- vandist_et_al_sim(n=100, a0=15, alpha=2.333, sigma0=.125, bias_sd=1)
ggplot(accuracies) + geom_line(aes(x=factor(block), y=acc, group=factor(proplab), colour=factor(proplab)))

plot_label(alpha=2.3333, proplab=.5, a0=10, bias=2)
plot_clusters(alpha=2.3333, proplab=.5, a0=10)
plot_guesses(alpha=2.3333, proplab=.5, a0=10, bias=rnorm(1))
plot_guesses(alpha=2.3333, proplab=.5, a0=10, bias=-2)
plot_guesses(alpha=2.3333, proplab=.25, a0=10)
plot_guesses(alpha=1, proplab=.25, a0=10)
plot_clusters(alpha=1, a0=10)
# }}}1

# {{{1 Simulations
# {{{2 Run sims
# Here we want to check for the effect of labeled items.
runs <- expand.grid(proplab=c(.5, 1),
                    lambda0=c(.5, 1),
                    a0=c(1, 5, 10, 15),
                    bias_sd=c(0,1,2),
                    tau=c(0, .05, .15),
                    sigma0=c(.05, .125, .25, .5),
                    alpha=c(1, .7/.3),
                    numtrials=c(800))

runs[runs$proplab==1,]$numtrials <- 400
nreps <- 100
sims <- as.data.table(run_vandist_sims(runs, nreps))
sims[proplab==1, block:=floor((trialnum-1)/40)+1]

accuracies <- sims[,list(acc=mean(hit)), by=c(names(runs), 'block')
#accuracies <- ddply(sims.reindexed, c(names(runs), 'block'), summarise, acc=mean(hit))
accuracies <- read.csv("vandistsims.csv")
write.csv(accuracies, file="vandistsims.csv")
# }}}2

# {{{2 Error function
expone <- c(.73, .80, .79, .87, .87, .87, .90, .88, .93, .91)
exptwo <- c(.72, .81, .8,  .89, .92, .89, .94, .93, .95, .94)

vandist.sqerr <- function(block) {
  ifelse (block['proplab']==1,
          block['acc']-exptwo[block['block']],
          block['acc']-expone[block['block']]) ^ 2
}

accuracies$vandisterr <- apply(accuracies, 1, vandist.sqerr)

params <- names(runs)[2:(length(names(runs))-1)]
subjfits <- ddply(accuracies, params, summarise, err=sum(vandisterr))

head(subjfits[with(subjfits, order(err)),])

# }}}2

# {{{2 Plot sims
ggplot(subset(accuracies, alpha==1)) + geom_line(aes(x=factor(block), y=acc, group=proplab, colour=factor(proplab))) + facet_wrap(~ lambda0+ a0+ bias_sd+ tau+ sigma0)
ggsave("figs/vandistsims.pdf")
# }}}2

#
# vi: foldmethod=marker
#
