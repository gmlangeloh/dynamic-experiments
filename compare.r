gmean <- function(data) {
  return(exp(mean(log(data))))
}

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

algorithmname_noheuristic <- function(path) {
  name_heuristic <- algorithmname(path)
  name <- strsplit(name_heuristic, "-", fixed=T)
  l <- length(name[[1]])
  fullname <- name[[1]][1:(l-1)]
  return(paste(fullname, collapse='-'))
}

algorithmnames_noheuristic <- function(paths) {
  lapply(paths, algorithmname_noheuristic)
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

compare_heuristics <- function(heuristic1, heuristic2) {
  patternh1 <- paste(c('*', heuristic1, '.out'), collapse='')
  patternh2 <- paste(c('*', heuristic2, '.out'), collapse='')
  filesh1 <- list.files(path="raw-results", pattern=patternh1, full.names=TRUE, recursive=FALSE)
  filesh2 <- list.files(path="raw-results", pattern=patternh2, full.names=TRUE, recursive=FALSE)
  df <- data.frame(matrix(ncol=8, nrow=0))
  colnames(df) <- c("Algorithm1, Algorithm2, no_t, t, polys, mons, deg, reds")
  for (i in 1:length(filesh1)) {
    for (j in 1:length(filesh2)) {
      if (algorithmname_noheuristic(filesh1[[i]]) == algorithmname_noheuristic(filesh2[[j]])) {
        mf <- compare(filesh1[[i]], filesh2[[j]])
        name1 <- algorithmname(filesh1[[i]])
        name2 <- algorithmname(filesh2[[j]])
        nf <- data.frame("Algorithm1"=name1, "Algorithm2"=name2, "no_t"=nrow(mf), "t"=gmean(mf$t_ratio), "polys"=gmean(mf$polys_ratio), "mons"=gmean(mf$mons_ratio), "deg"=gmean(mf$deg_ratio), "reds"=gmean(mf$reds_ratio))
        df <- rbind(df, nf)
      }
    }
  }
  library(xtable)
  print(xtable(df), include.rownames=F)
  return(df)
}

compare_all_heuristics <- function() {
  compare_heuristics('betti', 'mixed')
  compare_heuristics('hilbert', 'betti')
  compare_heuristics('hilbert', 'mixed')
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

count_str <- function(col1, col2) {
  smaller <- sum(col1 < col2)
  equals <- sum(col1 == col2)
  larger <- sum(col1 > col2)
  paste(c(toString(smaller), '/', toString(equals), '/', toString(larger)), collapse='')
}

compare_all_wrt_param <- function(parameter) {
  static <- "raw-results/static.out"
  hilbert <- c(static, list.files(path="raw-results", pattern="*hilbert.out", full.names=TRUE, recursive=FALSE))
  mixed <- c(static, list.files(path="raw-results", pattern="*mixed.out", full.names=TRUE, recursive=FALSE))
  total_algorithms <- length(hilbert)

  #Ratio matrix for this parameter --- row algorithm / col algorithm
  m <- matrix(ncol=total_algorithms, nrow=total_algorithms)
  mcount <- matrix(ncol=total_algorithms, nrow=total_algorithms)
  param1 <- paste(parameter, '.x', sep='')
  param2 <- paste(parameter, '.y', sep='')

  for (i in 1:total_algorithms) {
    for (j in 1:total_algorithms) {
      if (i < j) { #Hilbert
        algo1 <- read.table(hilbert[[i]], header=F)
        algo2 <- read.table(hilbert[[j]], header=F)
      }
      else { #Mixed
        algo2 <- read.table(mixed[[i]], header=F)
        algo1 <- read.table(mixed[[j]], header=F)
      }
      mf <- merge(algo1, algo2, by='V1')
      mf <- mf[which(mf$V4.x != Inf & mf$V4.y != Inf), ]
      col1 <- mf[[param1]]
      col2 <- mf[[param2]]
      res <- col1 / col2
      res[is.nan(res)] <- 1
      m[i, j] <- gmean(res)
      mcount[i, j] <- count_str(col1, col2)
    }
  }

  #Make a nice data frame with named rows/columns and print
  names <- algorithmnames_noheuristic(hilbert)
  algorithms <- lapply(names, function(x) paste('row-', x, sep=''))
  df <- data.frame(m)
  colnames(df) <- c(names)
  rownames(df) <- c(algorithms)
  print(xtable(df, caption=parameter, digits=2), include.rownames=T)

  df2 <- data.frame(mcount)
  colnames(df2) <- c(names)
  rownames(df2) <- c(algorithms)
  print(xtable(df2, caption=parameter), include.rownames=T)
  return(df)
}

compare_all <- function() {
  compare_all_wrt_param('V2') #Time
  compare_all_wrt_param('V4') #Basis size
  compare_all_wrt_param('V5') #Monomials in basis
  compare_all_wrt_param('V6') #Degrees
  compare_all_wrt_param('V7') #s-reductions
  return()
}

count_timeouts <- function() {
  static <- "raw-results/static.out"
  static_res <- read.table(static, header=F)
  cat("static", nrow(static_res[which(static_res$V4 != Inf), ]), "\n")

  hilbert <- c(static, list.files(path="raw-results", pattern="*hilbert.out", full.names=TRUE, recursive=FALSE))
  total_algorithms <- length(hilbert)

  for (i in 1:total_algorithms) {
    algo <- read.table(hilbert[[i]], header=F)
    cat(hilbert[[i]], nrow(algo[which(algo$V4 != Inf), ]), "\n")
  }
}

compare_all_size_degree <- function() {
  static <- "raw-results/static.out"
  hilbert <- c(static, list.files(path="raw-results", pattern="*hilbert.out", full.names=TRUE, recursive=FALSE))
  total_algorithms <- length(hilbert)

  #Ratio matrix for this parameter --- row algorithm / col algorithm
  m <- matrix(ncol=total_algorithms, nrow=total_algorithms)
  #mcount <- matrix(ncol=total_algorithms, nrow=total_algorithms)

  for (i in 1:total_algorithms) {
    for (j in 1:total_algorithms) {
      #Read results for each algorithm
      if (i < j) { #Basis size
          param1 <- 'V4.x'
          param2 <- 'V4.y'
      } else {
          param1 <- 'V6.x'
          param2 <- 'V6.y'
      }
      algo1 <- read.table(hilbert[[i]], header=F)
      algo2 <- read.table(hilbert[[j]], header=F)
      mf <- merge(algo1, algo2, by='V1')
      mf <- mf[which(mf$V4.x != Inf & mf$V4.y != Inf), ]
      col1 <- mf[[param1]]
      col2 <- mf[[param2]]
      res <- col1 / col2
      res[is.nan(res)] <- 1
      m[i, j] <- gmean(res)
      #mcount[i, j] <- count_str(col1, col2)
    }
  }

  #Make a nice data frame with named rows/columns and print
  names <- algorithmnames_noheuristic(hilbert)
  algorithms <- lapply(names, function(x) paste('row-', x, sep=''))
  df <- data.frame(m)
  colnames(df) <- c(names)
  rownames(df) <- c(algorithms)
  print(xtable(df, caption="size and degree", digits=2), include.rownames=T)

  #df2 <- data.frame(mcount)
  #colnames(df2) <- c(names)
  #rownames(df2) <- c(algorithms)
  #print(xtable(df2, caption=parameter), include.rownames=T)
  return(df)
}
