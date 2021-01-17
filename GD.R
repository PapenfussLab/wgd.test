# Genome doubling test
#
# Author: Lachlan Mcintosh
# Conceived by: Lachlan Mcintosh and Tony Papenfuss

parent_dir <- "."
all_files <- list.files(parent_dir)
all_files <- all_files[[1]]
all_files
library(ggplot2)
library(gridExtra)

# a function to set the number of xticks when using ggplot
number_ticks <- function(n) {
  function(limits)
    pretty(limits, n)
}

load_centromeres <- function(centromeres_file) {
  centromeres <- read.delim(centromeres_file, header = FALSE)
  colnames(centromeres) <- c("chromosome", "x")
  centromeres["x"] <- centromeres["x"] * 10^6
  return(centromeres)
}

centromeres <- load_centromeres("centromeres.txt")

get_CN_track <- function(mat, xmax, centromeres = centromeres) {
  break.from <- 0
  break.to <- 3 * 10^9
  break.by <- 50 * 10^6
  number.of.intervals <- 20
  
  g <- ggplot(mat) +
    geom_segment(aes(x = start.pos, xend = end.pos, y = A, yend = A,
                     col = chromosome), size = 2) +
    geom_segment(aes(x = start.pos, xend = end.pos, y = B, yend = B), size = 1) +
    geom_vline(data = centromeres, aes(xintercept = x, col = chromosome)) +
    
    # from ?formula
    # There are two special interpretations of . in a formula.
    # The usual one is in the context of a data argument
    # of model fitting functions and means
    # ‘all columns not otherwise in the formula’: see terms.formula.
    facet_grid(. ~ chromosome, scales = "free_x", space = "free") +
    
    scale_x_continuous(breaks = seq(break.from, break.to, break.by)) +
    scale_y_continuous(breaks = number_ticks(number.of.intervals)) +
    coord_cartesian(ylim = c(0, xmax)) +
    
    theme_bw() +
    theme(axis.text.x = element_blank(), legend.position = "none") +
    xlab("50 MB ticks") +
    ylab("major/minor CN") +
    ggtitle("CN track")
  
  return(g)
}

likelihoodGDany <-
  function(alpha,
           simpledata,
           half = FALSE,
           FA = TRUE,
           maxN = 6,
           maxM = 3) {
    # likelihoodhalfGD <- function(alpha,simpledata,half = FALSE){
    likelihoods <- data.frame('N' = numeric(),
                              'M' = numeric(),
                              'likelihood' = numeric())
    mylikelihood <-
      lapply(simpledata[1, c("A", "B")], function(x)
        getmylikelihood(x, 0,-1, FA = F))
    mylikelihood_nogd <-
      lapply(simpledata[1, c("A", "B")], function(x)
        getmylikelihood(x, 0,-1, FA = F))
    
    for (N in 0:maxN) {
      lowest <- -1
      if (half)
        lowest <- 0
      for (M in lowest:maxM) {
        # should allow a setting to free this up:
        # a and b are not used, safe to delete?
        a <- alpha
        b <- 1 - 2 * alpha
        c <- alpha
        
        N <- as.character(N)
        M <- as.character(M)
        
        for (i in seq_len(nrow(simpledata))) {
          if (sum(is.na(simpledata[i, c("A", "B")])) > 0) {
            #print(paste("FAILED",N,M,alpha,simpledata[i,]))
            new <- c(-Inf,-Inf,-Inf,-Inf)
            if (i == 1)
              new_likelihood <- new
            if (i > 1)
              new_likelihood <- rbind(new_likelihood, new)
            next
          }
          
          if (half) {
            if (min(simpledata[i, c("A", "B")]) == 0) {
              mylikelihood$A <-
                getmylikelihood(min(simpledata[i, c("A", "B")]), N, M, FA =
                                  F)
              mylikelihood$B <-
                getmylikelihood(max(simpledata[i, c("A", "B")]), N, M, FA =
                                  FA)
              mylikelihood_nogd$A <-
                getmylikelihood(min(simpledata[i, c("A", "B")]), N, "-1", FA =
                                  F)
              mylikelihood_nogd$B <-
                getmylikelihood(max(simpledata[i, c("A", "B")]), N, "-1", FA =
                                  FA)
            } else{
              mylikelihood <- lapply(simpledata[i, c("A", "B")], function(x)
                getmylikelihood(x, N, M, FA = F))
              mylikelihood_nogd <-
                lapply(simpledata[i, c("A", "B")], function(x)
                  getmylikelihood(x, N, "-1", FA = F))
            }
          } else{
            if (min(simpledata[i, c("A", "B")]) == 0) {
              mylikelihood$A <-
                getmylikelihood(min(simpledata[i, c("A", "B")]), N, M, FA =
                                  F)
              mylikelihood$B <-
                getmylikelihood(max(simpledata[i, c("A", "B")]), N, M, FA =
                                  FA)
              mylikelihood_nogd <- mylikelihood
            } else{
              mylikelihood <- lapply(simpledata[i, c("A", "B")], function(x)
                getmylikelihood(x, N, M, FA = F))
              mylikelihood_nogd <- mylikelihood
            }
          }
          if (is.null(mylikelihood$A) |
              is.null(mylikelihood$B) |
              is.null(mylikelihood_nogd$A) |
              is.null(mylikelihood_nogd$B)) {
            #print(paste("FAILED",N,M,alpha,simpledata[i,]))
            new <- c(-Inf,-Inf,-Inf,-Inf)
            if (i == 1)
              new_likelihood <- new
            if (i > 1)
              new_likelihood <- rbind(new_likelihood, new)
            next
          }
          val <- sapply(mylikelihood, function(x)
            log(as.numeric(eval(
              parse(text = x)
            ))))
          val_nogd <- sapply(mylikelihood_nogd, function(x)
            log(as.numeric(eval(
              parse(text = x)
            ))))
          if (i == 1)
            new_likelihood <- c(val, val_nogd)
          if (i > 1)
            new_likelihood <-
            rbind(new_likelihood, c(val, val_nogd))
        }
        
        if (half) {
          optA <- new_likelihood[, 1] + new_likelihood[, 4]
          optB <- new_likelihood[, 2] + new_likelihood[, 3]
          best <- optB
          best[which(optA > optB)] <- optA[which(optA > optB)]
        } else{
          best <- new_likelihood[, 1] + new_likelihood[, 2]
        }
        likelihood <- sum(best)
        likelihoods <-
          rbind(likelihoods,
                data.frame(
                  'N' = N,
                  'M' = M,
                  'likelihood' = likelihood
                ))
      }
    }
    LIKES <- matrix(rep(NA), nrow = (maxN + 1), ncol = (maxM + 2))
    for (N in 0:(maxN)) {
      for (M in-1:maxM) {
        temp <- likelihoods[which(likelihoods$N == N),]
        temp <- temp[which(temp$M == M),]
        if (nrow(temp) > 0) {
          LIKES[N + 1, M + 2] <- as.numeric(temp[, "likelihood"])
        }
      }
    }
    # print(LIKES)
    return(LIKES)
  }

getmylikelihood <- function(x, N, M, FA,
                            terms.path = file.path("GD", "terms")) {
  if (FA) {
    FA <- "True"
  } else {
    FA <- "False"
  }
  
  dir <- paste0("BTrue_FA", FA, "_N", N, "_M", M)
  file <- paste0("c", x, "_12_dec")
  filename <- file.path(terms.path, dir, file)
  
  # should this be an error or just a warning?
  if (!file.exists(filename)) {
    stop(filename, " can not be opened.")
  }
  
  return(readLines(filename))
}

outputdf <- data.frame(
  "name" = character(),
  "lmax" = numeric(),
  "alphamax" = numeric(),
  "indexmaxN" = numeric(),
  "indexmaxM" = numeric(),
  "lmax_half" = numeric(),
  "alphamax_half" = numeric(),
  "indexmaxN_half" = numeric(),
  "indexmaxM_half" = numeric(),
  "lmax_nogd" = numeric(),
  "alphamax_nogd" = numeric(),
  "indexmax_nogd" = numeric()
)

# what is this file and is it still / will be used?
# all_files = "/Users/lmcintosh/Downloads/CA004-8.CNV.txt"

i <- 0
for (file in all_files) {
  try({
    print(file)
    i <- i + 1
    print(i)
    # UNCOMMENT IN GENERAL
    # data = try(read.delim(pCaste(parent_dir,file,sep="/")))
    data <- read.delim(paste(parent_dir, file, sep = "/"))
    # data = read.delim(file)
    if (is.null(nrow(data))) {
      print("Could not load this file")
      next
    }
    chr_levels <- as.character(sort(unique(as.numeric(
      data$chromosome
    ))))
    chr_levels <- c(chr_levels[1:(length(chr_levels) - 1)], "X")
    #assumes that the only non numeric chromosome is "X"
    
    data$chromosome <- factor(data$chromosome, levels = chr_levels)
    
    # now we want to make a simplified version of the data where every chromosomal arm is independent!
    data$A <- data$Minor.Copy.Number
    data$B <- data$Major.Copy.Number
    data$start.pos <- data$start
    data$end.pos <- data$end
    
    simple_data <-
      data[0, c("chromosome", "A", "B", "start.pos", "end.pos")]
    
    for (chr in chr_levels) {
      chr_data <- data[which(data$chromosome == chr),]
      chr_data_p <-
        chr_data[which(chr_data$end.pos < centromeres[which(centromeres$chromosome == chr), "x"]),]
      chr_data_q <-
        chr_data[which(chr_data$start.pos > centromeres[which(centromeres$chromosome == chr), "x"]),]
      chr_data_middle <-
        chr_data[which(chr_data$start.pos < centromeres[which(centromeres$chromosome == chr), "x"] &
                         chr_data$end.pos > centromeres[which(centromeres$chromosome == chr), "x"]),]
      
      if (nrow(chr_data_middle) > 0) {
        # that is there is a segment that overlaps the centromere:
        chr_data_middle_p <- chr_data_middle
        chr_data_middle_q <- chr_data_middle
        chr_data_middle_p$end.pos <-
          centromeres[which(as.character(centromeres$chromosome) == as.character(chr)), "x"]
        chr_data_middle_q$start.pos <-
          centromeres[which(centromeres$chromosome == chr), "x"]
        chr_data_p <- rbind(chr_data_p, chr_data_middle_p)
        chr_data_q <- rbind(chr_data_q, chr_data_middle_q)
      }
      
      # calculate the summary statistics for each arm:
      meanAp <-
        sum((chr_data_p$end - chr_data_p$start) * chr_data_p$A) / sum(chr_data_p$end -
                                                                        chr_data_p$start)
      meanBp <-
        sum((chr_data_p$end - chr_data_p$start) * chr_data_p$B) / sum(chr_data_p$end -
                                                                        chr_data_p$start)
      startp <- min(chr_data_p$start.pos, na.rm = TRUE)
      endp <- max(chr_data_p$end.pos, na.rm = TRUE)
      
      meanAq <-
        sum((chr_data_q$end - chr_data_q$start) * chr_data_q$A) / sum(chr_data_q$end -
                                                                        chr_data_q$start)
      meanBq <-
        sum((chr_data_q$end - chr_data_q$start) * chr_data_q$B) / sum(chr_data_q$end -
                                                                        chr_data_q$start)
      startq <- min(chr_data_q$start.pos, na.rm = TRUE)
      endq <- max(chr_data_q$end.pos, na.rm = TRUE)
      
      # put the learnt information in a data frame:
      simple_data <-
        rbind(
          simple_data,
          data.frame(
            "chromosome" = chr,
            "A" = round(meanAp, 0),
            "B" = round(meanBp, 0),
            "start.pos" = startp,
            "end.pos" = endp
          )
        )
      simple_data <-
        rbind(
          simple_data,
          data.frame(
            "chromosome" = chr,
            "A" = round(meanAq, 0),
            "B" = round(meanBq, 0),
            "start.pos" = startq,
            "end.pos" = endq
          )
        )
    }
    g1 <- get_CN_track(data, max(data[, c("A", "B")], na.rm = TRUE))
    simple_data <-
      simple_data[which(!is.nan(simple_data$A)),] # what is the purpose of this line?
    g2 <-
      get_CN_track(simple_data, max(simple_data[, c("A", "B")], na.rm = TRUE))
    print(grid.arrange(g1, g2))
    
    # simple_data$zero = apply(simple_data[,c("A","B")],1,min)
    #
    # zero_data = simple_data[which(simple_data$zero == 0),]
    # non_zero_data = simple_data[which(simple_data$zero > 0),]
    #
    # X = c(non_zero_data[,"A"],non_zero_data[,"B"],apply(zero_data[,c("A","B")],1,min))
    #
    # sizeX = length(X)
    # X <- table(X) #[which(X>0)])
    # X <- X/sum(X)
    
    X <- table(c(simple_data[, c("A", "B")])) #[which(X>0)])
    X <- X / sum(X)
    X#
    # # Xzero are the ones that are in their final state protected from zero.
    # Xzero = apply(zero_data[,c("A","B")],1,max)
    # sizeXzero = length(Xzero)
    # Xzero <- table(X) #[which(X>0)])
    # Xzero <- X/sum(Xzero)
    
    alphamax_nogd <- 0
    lmax_nogd <- -Inf
    
    alphamax <- 0
    lmax <- -Inf
    
    alphamax_half <- 0
    lmax_half <- -Inf
    
    for (alpha in 1:10 / 50) {
      # likes = likelihoodGD(alpha,X,sizeX,Xzero,sizeXzero)
      # likes_half = likelihoodhalfGD(alpha,simple_data)
      
      likes <- likelihoodGDany(alpha, simple_data, half = F, FA = F)
      likes_half <-
        likelihoodGDany(alpha, simple_data, half = T, FA = F)
      
      l <- max(likes[, 2:ncol(likes)], na.rm = TRUE)
      l_nogd <- max(likes[, 1], na.rm = TRUE)
      
      l_half <- max(likes_half, na.rm = TRUE)
      if (l > lmax) {
        lmax <- l
        alphamax <- alpha
        indexmax <-
          which(likes == max(likes[, 2:ncol(likes)], na.rm = TRUE), arr.ind = TRUE)
        # the first row is for no GD....
      }
      if (l_half > lmax_half) {
        lmax_half <- l_half
        alphamax_half <- alpha
        indexmax_half <-
          which(likes_half == max(likes_half, na.rm = TRUE), arr.ind = TRUE)
      }
      if (l_nogd > lmax_nogd) {
        lmax_nogd <- l_nogd
        alphamax_nogd <- alpha
        indexmax_nogd <-
          which(likes == max(likes[, 1], na.rm = TRUE), arr.ind = TRUE)
      }
      print(c(lmax, alphamax, indexmax[1] - 1, indexmax[2] - 2))
      print(c(
        lmax_half,
        alphamax_half,
        indexmax_half[1] - 1,
        indexmax_half[2] - 2
      ))
      print(c(lmax_nogd, alphamax_nogd, indexmax_nogd[1] - 1))
      # print(likes)
      # print(likes_half)
      print(alpha)
      print(X)
    }
    
    print(c(lmax, alphamax, indexmax[1] - 1, indexmax[2] - 1))
    print(c(
      lmax_half,
      alphamax_half,
      indexmax_half[1] - 1,
      indexmax_half[2] - 2
    ))
    print(c(lmax_nogd, alphamax_nogd, indexmax_nogd[1] - 1))
    print(X)
    
    outputdf <- rbind(
      outputdf,
      data.frame(
        "name" = file,
        "lmax" = lmax,
        "alphamax" = alphamax,
        "indexmaxN" = indexmax[1] - 1,
        "indexmaxM" = indexmax[2] - 2,
        "lmax_half" = lmax_half,
        "alphamax_half" = alphamax_half,
        "indexmaxN_half" = indexmax_half[1] -
          1,
        "indexmaxM_half" = indexmax_half[2] - 2,
        "lmax_nogd" = lmax_nogd,
        "alphamax_nogd" = alphamax_nogd,
        "indexmax_nogd" = indexmax_nogd[1] - 1
      )
    )
    
    if (lmax > lmax_nogd & lmax > lmax_half) {
      print("the GD")
    } else if (lmax_half > lmax_nogd & lmax_half > lmax) {
      print("half the GD")
    } else{
      print("no GD")
    }
  })
}
outputdf <- outputdf[, 1:12]
colnames(outputdf)
outputdf$GD <-
  outputdf$lmax > outputdf$lmax_nogd# & outputdf$lmax > outputdf$lmax_half
# outputdf$GD_half <- outputdf$lmax_half > outputdf$lmax_nogd & outputdf$lmax_half > outputdf$lmax
outputdf$no_GD <-
  outputdf$lmax_nogd >= outputdf$lmax #& outputdf$lmax_nogd > outputdf$lmax_half

c(sum(outputdf$GD), sum(outputdf$GD_half), sum(outputdf$no_GD))

outputdf[, c("lmax", "lmax_nogd", "GD", "no_GD")]

k <- 3
n <- 92
outputdf$AIC_full <- 2 * k - 2 * outputdf$lmax
# outputdf$AIC_half <- 2*(k+23)-2*outputdf$lmax_half
outputdf$AIC_null <- 2 * (k - 1) - 2 * outputdf$lmax_nogd

outputdf$GD_AIC <- outputdf$AIC_full <  outputdf$AIC_null
# outputdf$GD_half_AIC <- outputdf$AIC_half < outputdf$AIC_full & outputdf$AIC_half < outputdf$AIC_null
outputdf$no_GD_AIC <-  outputdf$AIC_null <= outputdf$AIC_full

# AIC is founded on information theory:
# it offers a relative estimate of the information lost
# when a given model is used to represent the process that
# generates the data. In doing so, it deals with the trade-off
# between the goodness of fit of the model and the complexity of
# the model. We choose the candidate model that minimized
# the information loss.

outputdf$bestAIC <-
  apply(outputdf[, c("AIC_full", "AIC_null")], 1, min)
outputdf$rl_GD <- exp((outputdf$bestAIC - outputdf$AIC_full) / 2)
outputdf$rl_noGD <- exp((outputdf$bestAIC - outputdf$AIC_null) / 2)

write.table(outputdf, file = paste0("/Users/lmcintosh", "/outputdf.txt"))

all_files[which(!(all_files %in% outputdf$name))]
