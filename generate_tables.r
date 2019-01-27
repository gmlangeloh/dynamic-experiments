library(xtable)

files <- list.files(path="results", pattern="*.out", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  t <- read.table(x, header=TRUE)
  print(xtable(t), include.rownames=F)
})
