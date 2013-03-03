
library(ggplot2)

run_anderson_once <- function(argstring="") {
  output_columns <- c("type","x","y","label","clust")
  command <- paste("./testanderson", argstring)
  read.csv(pipe(command), col.names=output_columns, header=F)
}
run_anderson_once()

plot_clusters <- function(...) {
  simresult <- run_anderson_once(...)
  stims <- subset(simresult, type=="STIM")
  ggplot(stims) + geom_point(aes(x=x, y=y, colour=factor(clust)))
}
plot_clusters("--alphaparam=2.333 --bias=-1")

plot_inference <- function(...) {
  simresult <- run_anderson_once(...)
  stims <- subset(simresult, type=="INFER")
  ggplot(stims) + geom_point(aes(x=x, y=y, colour=label))
}
plot_inference("--alphaparam=2.333 --bias=-1")



