
source("simfunctions.R")

# {{{1 Try single runs

simulate_anderson(task="tvtask", alpha=2.333, nlab=-1, bias=0, a0=10, alab=1, lambda0=1, sigma0=0.15, tau=0.05, plotting=T, order="interspersed", echo=F)

# }}}1
#
# {{{1 Just simulate exp 1
runs <- expand.grid(task="tvtask",
                    nlab=c(0,4,16,-1), 
                    order=c("interspersed"),
                    sigma0=c(.25),
                    a0=c(16),
                    alab=c(.5),
                    lambda0=c(1),
                    alpha=c(.7/.3),
                    tau=c(.05), 
                    bias_sd=c(0),
                    encoding=c("encodeactual"))
runs <- subset(runs, ! ((alpha==1 & order!="interspersed") | (alpha==1 & nlab==4)))
nrow(runs)

nreps <- 300
ofile <- "/dev/null"
sims <- run_sims(runs, nreps, ofile)
sims$nlab[sims$nlab==-1] <- Inf

counts <- ddply(sims, .(nlab, alpha, sigma0, a0, alab, lambda0, tau, order, bias_sd, encoding), function(x) summary(x$BestFit))
counts$twod <- counts$"2D"
counts
# }}}1

# {{{1 Run sims across conditions.
runs <- expand.grid(task="tvtask",
                    nlab=c(0,4,16,-1), 
                    order=c("interspersed", "labeledfirst", "labeledlast"),
                    sigma0=c(.125),
                    a0=c(30),
                    alab=c(1),
                    lambda0=c(1),
                    alpha=c(1, .7/.3),
                    tau=c(.05), 
                    bias_sd=c(0, 1.5),
                    encoding=c("encodeactual"))
runs <- subset(runs, ! ((alpha==1 & order!="interspersed") | (alpha==1 & nlab==4) | (order=="labeledfirst" & nlab!=16) | (order=="labeledlast" & nlab==4)))
nrow(runs)

nreps <- 101
ofile <- "data/forpub"
sims <- read.csv(ofile)
#sims <- run_sims(runs, nreps, ofile)
sims$nlab[sims$nlab==-1] <- Inf

counts <- ddply(sims, .(nlab, alpha, sigma0, a0, alab, lambda0, tau, order, bias_sd, encoding), function(x) summary(x$BestFit))
counts$twod <- counts$"2D"
counts
# }}}1

# {{{1 Plots
# {{{2 Setup
counts.melted <- melt(counts,
                      id=c("nlab", "alpha", "tau", "order", "encoding", "bias_sd", "sigma0", "a0"),
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
ggplot(plotdata) + geom_line(aes(x=factor(nlab), y=count, group=bias_sd, linetype=factor(bias_sd)))  + facet_grid(~bestfit)
ggsave("exp1.pdf")
# }}}2

# {{{2 Experiment 2
identify.exp2.cond <- function(row) {
    if (row["order"]!="interspersed") return(NA)
    nlab <- as.numeric(row["nlab"]) 
    if (is.infinite(nlab)) nlab <- "All"
    alpha <- as.numeric(row["alpha"]) 
    if (alpha==1) {
        return(paste(nlab, "-lab-resp", sep=""))
    } else if (nlab==16) {
        return(paste(nlab, "-lab-noresp", sep=""))
    } else return(NA)
}

plot.data <- subset(counts.melted, order=="interspersed")
plot.data$cond <- factor(apply(plot.data, 1, identify.exp2.cond))
plot.data$cond <- factor(plot.data$cond, levels=levels(plot.data$cond)[c(1,3,2,4)])
plot.data <- subset(plot.data, ! is.na(cond) )

ggplot(plot.data) + geom_line(aes(x=cond, y=count, group=bias_sd, linetype=factor(bias_sd)))  + facet_grid(~bestfit)
ggsave("exp2.pdf")
# }}}2

# {{{2 Experiment 3
identify.exp3.cond <- function(row) {
    nlab <- as.numeric(row["nlab"]) 
    order <- row["order"]
    if (nlab==0 & order=="interspersed") return("US")
    if (order=="interspersed") return(NA)
    alpha <- as.numeric(row["alpha"]) 
    if (alpha==1) return(NA)
    if (is.infinite(nlab) & order=="labeledlast") return("FS-lablast")
    if (nlab == 16) {
      if (order=="labeledfirst") return("SS-labfirst")
      if (order=="labeledlast") return("SS-lablast")
    }
    return(NA)
}
plot.data <- subset(counts.melted, alpha==highalpha)
plot.data$cond <- factor(apply(plot.data, 1, identify.exp3.cond))
plot.data$cond <- factor(plot.data$cond, levels=levels(plot.data$cond)[c(4,3,2,1)])
plot.data <- subset(plot.data, ! is.na(cond))

ggplot(plot.data) + geom_line(aes(x=cond, y=count, group=bias_sd, linetype=factor(bias_sd)))  + facet_grid(~bestfit)
ggsave("exp3.pdf")
# }}}2
# }}}1
