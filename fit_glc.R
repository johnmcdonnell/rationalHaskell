
source("simfunctions.R")

# {{{1 Try single runs

test.one.run <- function() {
  simulate_anderson(task="tvtask", alpha=2.333, nlab=-1, bias=0, a0=10, alab=1, lambda0=1, sigma0=0.15, tau=0.05, plotting=T, order="interspersed", echo=F)
}

# }}}1
#
# {{{1 Just simulate unsupervised cond of Exp 1
runs <- expand.grid(task="tvtask",
                    nlab=c(0), 
                    order=c("interspersed"),
                    sigma0=c(.125, .25),
                    a0=c(1, 10, 30),
                    alab=c(1),
                    lambda0=c(1),
                    alpha=c(.7/.3),
                    tau=c(0.05), 
                    bias_sd=c(0, 1.5),
                    encoding=c("encodeactual"))
runs <- subset(runs, ! ((alpha==1 & order!="interspersed") | (alpha==1 & nlab==4)))
nrow(runs)
nreps <- 1000
ofile <- "data/testfirstcond_tau05.csv"

sims <- read.csv(ofile)
#sims <- run_sims(runs, nreps, ofile)
sims$nlab[sims$nlab==-1] <- Inf

counts <- ddply(sims, .(nlab, alpha, sigma0, a0, alab, lambda0, tau, order, bias_sd, encoding), function(x) summary(x$BestFit))
counts$twod <- counts$"2D"
counts

cond1.empirical.proportions <- c(twod=.1714286, Unimodal=0.257142, Bimodal=0.5714)
cond1Error <- function(counts) {
  denom <- counts$twod + counts$Unimodal + counts$Bimodal# + counts$Null
  (((counts$twod / denom) - cond1.empirical.proportions['twod'])^2 +
    ((counts$Unimodal / denom) - cond1.empirical.proportions['Unimodal'])^2 +
    ((counts$Bimodal / denom) - cond1.empirical.proportions['Bimodal'])^2 +
    ((counts$Null / denom) - 0)^2
    )
}


counts$err <- cond1Error(counts)

# What is the best fit without bias?
head(subset(counts[with(counts, order(err)),], bias_sd>0))

# What is the best fit with bias?
head(subset(counts[with(counts, order(err)),], bias_sd==0))

counts[with(counts, order(-err)),]
# }}}1

# {{{1 Run sims across conditions.
runs <- expand.grid(task="tvtask",
                    nlab=c(0,4,16,-1), 
                    order=c("interspersed", "labeledfirst", "labeledlast"),
                    sigma0=c(.125, .25),
                    a0=c(1, 10, 30),
                    alab=c(1),
                    lambda0=c(1),
                    alpha=c(1, .7/.3),
                    tau=c(0.05), 
                    bias_sd=c(0, 1.5),
                    encoding=c("encodeactual"))
runs <- subset(runs, ! (((nlab==4 | alpha==1) & order!="interspersed") | (alpha==1 & nlab==4) | (order=="labeledfirst" & nlab!=16)))
nrow(runs)

nreps <- 1000
ofile <- "data/newgrid_tau.csv"
ofile <- "data/forpubthree.csv"
#sims <- read.csv(ofile)

simulate.all <- function(runs, nreps, ofile) {
  sims <- run_sims(runs, nreps, ofile)
  sims$nlab[sims$nlab==-1] <- Inf
  counts <- ddply(sims, .(nlab, alpha, sigma0, a0, alab, lambda0, tau, order, bias_sd, encoding), function(x) summary(x$BestFit))
  counts$twod <- counts$"2D"
  counts
}

counts <- simulate.all(runs, nreps, ofile)

# }}}1

# {{{1 Fitting model fits to data
# {{{2 Experiment 1
# {{{3 Functions
get.expone.proportions <- function(df) {
      denom <- function(x) x$twod + x$Unimodal + x$Bimodal #+ x$Null
      get.condition.proportions <- function(cond.df) {
        denom <- denom(cond.df)
        data.frame(
          bimod = cond.df$Bimodal / denom,
          unimod = cond.df$Unimodal / denom,
          twod = cond.df$twod / denom)
      }
      
      ddply(df, .(nlab), get.condition.proportions)
}
score.run <- function(df) {
  bimod.props <- c(0.5714,.442857, .3, .3)
  unimod.props <- c(0.257142,.22857, .3, .2)
  twod.props <- c(.1714286, .32857, .4, .5)
  ssqerr <- sum(c(df$bimod[1]-bimod.props[1], df$unimod[1]-unimod.props[1], df$twod[1]-twod.props[1]) ^ 2)
  data.frame(sqerr=ssqerr)
}
rate.bimodality.effect <- function(df) {
  fourlab <- subset(df, nlab==4)
  sixteenlab <- subset(df, nlab==16)
  data.frame(bimodeffect=fourlab$bimod - sixteenlab$bimod)
}
# }}}3

params <- .(sigma0, a0, alab, lambda0, tau, order, bias_sd, encoding)
expone <- subset(counts, order=="interspersed" & alpha>1)
expone.condition.proportions <- ddply(expone, params, get.expone.proportions)

model.fits.to.data <- ddply(expone.condition.proportions, params, score.run)
model.fits.to.data[with(model.fits.to.data, order(-sqerr)),]

head(model.fits.to.data[with(model.fits.to.data, order(sqerr)),])
head(subset(model.fits.to.data[with(model.fits.to.data, order(sqerr)),], bias_sd==0))

bimod.effect <- ddply(expone.condition.proportions, params, rate.bimodality.effect)
subset(bimod.effect, bias_sd==1.5)
ddply(bimod.effect, .(sigma0), summarise, mean(bimodeffect))
ddply(bimod.effect, .(bias_sd), summarise, mean(bimodeffect))
# }}}2
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

# vi: foldmethod=marker
