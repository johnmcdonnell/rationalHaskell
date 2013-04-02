
#library(ggplot2)
source("simfunctions.R")

# By default, sigma is .5 (pooled sd / 3)
x <- run_anderson_once(task="vandist", a0=10, proplab=.5)
x

names(run_anderson_once())

plot_clusters <- function(...) {
  simresult <- run_anderson_once(task="vandist", ...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  print(stims)
  print(names(stims))
  ggplot(stims) + geom_point(aes(x=bimod, y=unimod, colour=factor(clust)))
}

plot_guesses <- function(...) {
  simresult <- run_anderson_once(task="vandist", ...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  print(stims)
  print(names(stims))
  ggplot(stims) + geom_point(aes(x=bimod, y=unimod, colour=guess))
}

plot_clusters(alpha=2.3333, proplab=.5, a0=10, bias=2)
plot_clusters(alpha=2.3333, proplab=.5, a0=10)
plot_guesses(alpha=2.3333, proplab=.5, a0=10, bias=rnorm(1))
plot_guesses(alpha=2.3333, proplab=.5, a0=10, bias=-2)
plot_guesses(alpha=2.3333, proplab=.25, a0=10)
plot_guesses(alpha=1, proplab=.25, a0=10)
plot_clusters(alpha=1, a0=10)

plot(3)
