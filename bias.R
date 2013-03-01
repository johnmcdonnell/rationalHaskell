
library(ggplot2)

run_anderson_once <- function(...) {
  output_columns <- c("type","x","y","label","clust")
  command <- paste("./randombias", ...)
  read.csv(pipe(command), col.names=output_columns, header=F)
}

names(run_anderson_once())
run_anderson_once()


plot_clusters <- function(...) {
  simresult <- run_anderson_once(...)
  stims <- subset(simresult, type=="STIM")[1:100,]
  ggplot(stims) + geom_point(aes(x=x, y=y, colour=factor(clust)))
}
plot_clusters()

plot_inference <- function(...) {
  simresult <- run_anderson_once(...)
  stims <- subset(simresult, type=="INFER")[1:100,]
  ggplot(stims) + geom_point(aes(x=x, y=y, colour=label))
}
plot_inference(2.3333)


