
library(ggplot2)

run_anderson_once <- function(...) {
  output_columns <- c("type","x","y","label","clust")
  command <- paste("./testanderson", ...)
  read.csv(pipe(command), col.names=output_columns, header=F)
}

names(run_anderson_once())

plot_clusters <- function(...) {
  simresult <- run_anderson_once(...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  ggplot(stims) + geom_point(aes(x=x, y=y, colour=factor(clust)))
}

plot_clusters(1, 0.1)

plot(3)
