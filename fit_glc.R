
source("simfunctions.R")

# {{{1 Try single runs

simulate_anderson(alpha=2.333, nlab=-1, bias=100, a0=10, lambda0=1, sigma0=0.15, tau=0.05, plotting=T, echo=F)
plot_anderson(alpha=1, lambda0=1, a0=10, nlab=16, bias=2)
# It looks like negative bias is weaker than positive.
#b <- rnorm(1, mean=0, sd=4)
simulate_anderson(alpha=2.333, nlab=-1, bias=7, sigma0=0.35, tau=0.05, plotting=T, echo=T)
simulate_anderson(alpha=2.333, nlab=-1, bias=7, tau=0.05, plotting=T, echo=T)
simulate_anderson(alpha=1, nlab=0, order="interspersed", encoding="encodeactual", bias=0, tau=0.05, plotting=F)
simulate_anderson(alpha=2.333, nlab=16, bias=0, plotting=T)
simulate_anderson(alpha=1, nlab=16, bias=0)

# }}}1

# {{{1 Run sims across conditions.
runs <- expand.grid(nlab=c(0,4,16,-1), 
                    order=c("interspersed", "labeledfirst", "labeledlast"),
                    sigma0=c(.15),
                    a0=c(10),
                    lambda0=c(1),
                    alpha=c(1, .7/.3),
                    tau=c(.05), 
                    bias_sd=c(0, 1, 2, 5),
                    encoding=c("encodeactual"))
runs <- subset(runs, ! ((alpha==1 & order!="interspersed") | (alpha==1 & nlab==4)))
nrow(runs)

nreps <- 500
ofile <- "bias.csv"
sims <- read.csv(ofile)
#sims <- run_sims(runs, nreps, ofile)
sims$nlab[sims$nlab==-1] <- Inf

counts <- ddply(sims, .(nlab, alpha, sigma0, a0, lambda0, tau, order, bias_sd, encoding), function(x) summary(x$BestFit))
#counts$ratio <- counts$Bimodal/counts$"2D"
#counts$params <- do.call(paste, c(counts[c("alpha", "tau")], sep = ":"))
counts$twod <- counts$"2D"
#print(ggplot(counts) + geom_line(aes(x=nlab, y=twod, group=params, linetype=factor(tau), colour=factor(alpha))) + facet_wrap(~order))
# }}}1

# {{{1 taking a look
semisup <- subset(counts, nlab %in% c(0,16) & bias_sd==0 & order=="interspersed")
countsquantified <- ddply(semisup, .(alpha, sigma0, a0, lambda0, tau, order, bias_sd, encoding), function(df) {
      unlab <- subset(df, nlab==0)
      lab <- subset(df, nlab==16)
      denom <- function(x) x$twod + x$Unimodal + x$Bimodal #+ x$Null
      unlabpercent <- unlab$Bimodal / denom(unlab)
      labpercent <- lab$Bimodal / denom(lab)
      data.frame(lab=labpercent, unlab=unlabpercent)})
ggplot(countsquantified) + geom_point(aes(x=lab, y=unlab, colour=factor(sigma0))) + geom_abline() + facet_grid(lambda0~a0)

# }}}1

# {{{1 Plots
# {{{2 Setup
counts.melted <- melt(counts,
                      id=c("nlab", "alpha", "tau", "order", "encoding", "bias_sd"),
                      measure.vars=c("2D", "Bimodal", "Unimodal", "Null"),
                      variable.name="bestfit",
                      value.name="count")

#counts.melted$groupingparam <- do.call(paste, c(counts.melted[c("bestfit", "encoding")], sep = ":"))
#counts.melted

highalpha <- subset(counts.melted, alpha>1)[1,'alpha']
lowalpha <- 1
# }}}2

# {{{2 Experiment 1
plotdata <- subset(counts.melted, alpha==highalpha & order=="interspersed")

ggplot(plotdata) + geom_line(aes(x=factor(nlab), y=count, group=bestfit))  + facet_wrap(bestfit~bias_sd)
# }}}2

# {{{2 Experiment 2
plotdata <- subset(counts.melted, nlab==16 & order=="interspersed")

ggplot(plotdata) + geom_line(aes(x=factor(alpha), y=count, group=bestfit))  + facet_wrap(bestfit~bias_sd)
# }}}2

# {{{2 Experiment 3
plotdata <- subset(counts.melted, nlab==16 & alpha==highalpha)

ggplot(plotdata) + geom_line(aes(x=factor(order), y=count, group=bestfit))  + facet_wrap(bestfit~bias_sd)
# }}}2
# }}}1
