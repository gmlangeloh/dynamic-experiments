lteqgt <- function(col1, col2) {

  print(sum(col1 < col2))
  print(sum(col1 == col2))
  print(sum(col1 > col2))

  print(sum(col1 - col2))

}

addcol <- function(dt, name, col1, col2) {
  dt[paste(name)] <- col1 / col2
  return(dt)
}

instances <- function() {
  i <- read.table("results/instances.out", header=T)
  names(i)[names(i) == 'name'] <- 'V1'
  return(i)
}

compare <- function(path1, path2) {

  #Merge the tables to make comparing stats easier
  df1 <- read.table(path1, header=F)
  df2 <- read.table(path2, header=F)
  mf <- merge(df1, df2, by='V1')
  #Keep only stuff that don't timeout in either
  mf <- mf[which(mf$V4.x != Inf & mf$V4.y != Inf), ]
  mf <- merge(mf, instances(), by='V1')

  mf <- addcol(mf, "t_ratio", mf$V2.x, mf$V2.y)
  mf <- addcol(mf, "polys_ratio", mf$V4.x, mf$V4.y)
  mf <- addcol(mf, "mons_ratio", mf$V5.x, mf$V5.y)
  mf <- addcol(mf, "deg_ratio", mf$V6.x, mf$V6.y)
  mf <- addcol(mf, "reds_ratio", mf$V7.x, mf$V7.y)

  #print(summary(mf))

  return(mf)
}

algorithmname <- function(path) {
  dirs <- strsplit(path, "/", fixed=T)
  filename <- dirs[[1]][length(dirs[[1]])]
  algorithm <- strsplit(filename, ".", fixed=T)
  return(algorithm[[1]][1])
}

meanratios <- function() {
  files <- list.files(path="raw-results", pattern="*.out", full.names=TRUE, recursive=FALSE)
  df <- data.frame(matrix(ncol=8, nrow=0))
  colnames(df) <- c("Algorithm1, Algorithm2, no_t, t, polys, mons, deg, reds")
  for (i in 1:length(files)) {
    for (j in 1:i) {
      if (i != j) {
        mf <- compare(files[[i]], files[[j]])
        name1 <- algorithmname(files[[i]])
	name2 <- algorithmname(files[[j]])
        nf <- data.frame("Algorithm1"=name1, "Algorithm2"=name2, "no_t"=nrow(mf),"t"=mean(mf$t_ratio), "polys"=mean(mf$polys_ratio), "mons"=mean(mf$mons_ratio), "deg"=mean(mf$deg_ratio), "reds"=mean(mf$reds_ratio))
  	df <- rbind(df, nf)
      }
    }
  }
  library(xtable)
  print(xtable(df), include.rownames=F)
  return(df)
}

bettivsmixed <- function() {
  betti <- list.files(path="raw-results", pattern="*betti.out", full.names=TRUE, recursive=FALSE)
  mixed <- list.files(path="raw-results", pattern="*mixed.out", full.names=TRUE, recursive=FALSE)
  df <- data.frame(matrix(ncol=8, nrow=0))
  colnames(df) <- c("Algorithm1, Algorithm2, no_t, t, polys, mons, deg, reds")
  for (i in 1:length(betti)) {
    for (j in 1:length(mixed)) {
      mf <- compare(betti[[i]], mixed[[j]])
      name1 <- algorithmname(betti[[i]])
      name2 <- algorithmname(mixed[[j]])
      nf <- data.frame("Algorithm1"=name1, "Algorithm2"=name2, "no_t"=nrow(mf), "t"=mean(mf$t_ratio), "polys"=mean(mf$polys_ratio), "mons"=mean(mf$mons_ratio), "deg"=mean(mf$deg_ratio), "reds"=mean(mf$reds_ratio))
      df <- rbind(df, nf)
    }
  }
  library(xtable)
  print(xtable(df), include.rownames=F)
  return(df)
}

hilbertvsmixed <- function() {
  hilbert <- list.files(path="raw-results", pattern="*hilbert.out", full.names=TRUE, recursive=FALSE)
  mixed <- list.files(path="raw-results", pattern="*mixed.out", full.names=TRUE, recursive=FALSE)
  df <- data.frame(matrix(ncol=8, nrow=0))
  colnames(df) <- c("Algorithm1, Algorithm2, no_t, t, polys, mons, deg, reds")
  for (i in 1:length(hilbert)) {
    for (j in 1:length(mixed)) {
      mf <- compare(hilbert[[i]], mixed[[j]])
      name1 <- algorithmname(hilbert[[i]])
      name2 <- algorithmname(mixed[[j]])
      nf <- data.frame("Algorithm1"=name1, "Algorithm2"=name2, "no_t"=nrow(mf), "t"=mean(mf$t_ratio), "polys"=mean(mf$polys_ratio), "mons"=mean(mf$mons_ratio), "deg"=mean(mf$deg_ratio), "reds"=mean(mf$reds_ratio))
      df <- rbind(df, nf)
    }
  }
  library(xtable)
  print(xtable(df), include.rownames=F)
  return(df)
}

min_max_sizes <- function() {
  files <- list.files(path="raw-results", pattern="*.out", full.names=TRUE, recursive=FALSE)
  m <- read.table("raw-results/static.out", header=FALSE)
  m$V2 <- NULL
  m$V3 <- NULL
  m$V5 <- NULL
  m$V6 <- NULL
  m$V7 <- NULL
  names(m)[names(m) == 'V4'] <- 'raw-results/static.out'
  lapply(files, function(x) {
    t <- read.table(x, header=FALSE)
    t$V2 <- NULL
    t$V3 <- NULL
    t$V5 <- NULL
    t$V6 <- NULL
    t$V7 <- NULL
    names(t)[names(t) == 'V4'] <- x
    m <<- merge(m, t, by='V1')
  })
  V1 <- m$V1
  m$V1 <- NULL
  m$minbasis <- apply(m, 1, FUN=min)
  m$maxbasis <- apply(m, 1, FUN=max)
  m$minarg <- colnames(m)[apply(m, 1, which.min)]
  m$V1 <- V1
  return(m[, c('V1', 'minbasis', 'minarg', 'maxbasis')])
}
