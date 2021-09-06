
library(tidyverse)

# I removed a check for the object to save time. Maybe revert to this version.
cluster_addProcessHistory <- function(object, ph) {
  if (!inherits(ph, "ProcessHistory"))
    stop("Argument 'ph' has to be of type 'ProcessHistory' or a class ",
         "extending it!")
  object@.processHistory[[(length(object@.processHistory) + 1)]] <- ph
  object
}


# Removed the update_feature_definintion check to make it run faster. Revert to this version if the script becomes too slow.
cluster_adjustRtimePeakGroups <- function(object, param = PeakGroupsParam(),
                                  msLevel = 1L) {
  subs <- integer(0)
  if (!length(subs)) {subs <- seq_along(fileNames(object))}
  
  nSamples <- length(subs)
  pkGrp <- .getPeakGroupsRtMatrix(
    peaks = chromPeaks(object, msLevel = msLevel),
    # peakIndex = .peakIndex(
    #         .update_feature_definitions(
    #             featureDefinitions(object), rownames(chromPeaks(object)),
    #             rownames(chromPeaks(object, msLevel = msLevel))))
    peakIndex = .peakIndex(featureDefinitions(object)),
    sampleIndex = subs,
    missingSample = nSamples - (nSamples * minFraction(param)),
    extraPeaks = extraPeaks(param)
  )
  colnames(pkGrp) <- basename(fileNames(object))[subs]
  pkGrp
}

cluster_do_adjustRtime_peakGroups <-
  function(peaks, peakIndex, rtime, minFraction = 0.9, extraPeaks = 1,
           smooth = c("loess", "linear"), span = 0.2,
           family = c("gaussian", "symmetric"),
           peakGroupsMatrix = matrix(ncol = 0, nrow = 0),
           subset = integer(), subsetAdjust = c("average", "previous")){
    
    cat("please be aware that clusterXCMS does not support subsets (e.g., pooled QCs). \nAny processing including subset information must be done on post-xcms tabular data")
    ## Translate minFraction to number of allowed missing samples.
    nSamples      <- length(subset)
    missingSample <- nSamples - (nSamples * minFraction)
    if (nrow(peakGroupsMatrix)) {
      rt <- peakGroupsMatrix
    } else
      rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, subset,
                                   missingSample)

    message("Performing retention time correction using ", nrow(rt),
            " peak groups.")
    
    ## Calculate the deviation of each peak group in each sample from its
    ## median
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)
    
    if (smooth == "loess") {
      mingroups <- min(colSums(!is.na(rt)))
      if (mingroups < 4) {
        smooth <- "linear"
        warning("Too few peak groups for 'loess', reverting to linear",
                " method")
      } else if (mingroups * span < 4) {
        span <- 4 / mingroups
        warning("Span too small for 'loess' and the available number of ",
                "peak groups, resetting to ", round(span, 2))
      }
    }
    
    rtdevsmo <- vector("list", nSamples)
    
    ## Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE
    warn.tweak.rt <- FALSE
    
    rtime_adj <- rtime
    ## Adjust samples in subset.
    for (i in seq_along(subset)) {
      i_all <- subset[i]              # Index of sample in whole dataset.
      pts <- na.omit(data.frame(rt = rt[, i], rtdev = rtdev[, i]))
      
      ## order the data.frame such that rt and rtdev are increasingly ordered.
      pk_idx <- order(pts$rt, pts$rtdev)
      pts <- pts[pk_idx, ]
      if (smooth == "loess") {
        lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span,
                                     degree = 1, family = family))
        
        rtdevsmo[[i]] <- xcms:::na.flatfill(
          predict(lo, data.frame(rt = rtime[[i_all]])))
        ## Remove singularities from the loess function
        rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                        quantile(abs(rtdevsmo[[i]]), 0.9,
                                 na.rm = TRUE) * 2] <- NA
        if (length(naidx <- which(is.na(rtdevsmo[[i]])))){
          rtdevsmo[[i]][naidx] <- suppressWarnings(
            approx(na.omit(data.frame(rtime[[i_all]], rtdevsmo[[i]])),
                   xout = rtime[[i_all]][naidx], rule = 2)$y
          )
          }
        
        ## Check if there are adjusted retention times that are not ordered
        ## increasingly. If there are, search for each first unordered rt
        ## the next rt that is larger and linearly interpolate the values
        ## in between (see issue #146 for an illustration).
        while (length(decidx <- which(diff(rtime[[i_all]] - rtdevsmo[[i]]) < 0))) {
          warn.tweak.rt <- TRUE  ## Warn that we had to tweak the rts.
          rtadj <- rtime[[i_all]] - rtdevsmo[[i]]
          rtadj_start <- rtadj[decidx[1]] ## start interpolating from here
          ## Define the
          next_larger <- which(rtadj > rtadj[decidx[1]])
          if (length(next_larger) == 0) {
            ## Fix if there is no larger adjusted rt up to the end.
            next_larger <- length(rtadj) + 1
            rtadj_end <- rtadj_start
          } else {
            next_larger <- min(next_larger)
            rtadj_end <- rtadj[next_larger]
          }
          ## linearly interpolate the values in between.
          adj_idxs <- (decidx[1] + 1):(next_larger - 1)
          incr <- (rtadj_end - rtadj_start) / length(adj_idxs)
          rtdevsmo[[i]][adj_idxs] <- rtime[[i_all]][adj_idxs] -
            (rtadj_start + (1:length(adj_idxs)) * incr)
        }
        
        rtdevsmorange <- range(rtdevsmo[[i]])
        if (any(rtdevsmorange / rtdevrange > 2))
          warn.overcorrect <- TRUE
      } else {
        if (nrow(pts) < 2) {
          stop("Not enough peak groups even for linear smoothing ",
               "available!")
        }
        ## Use lm instead?
        fit <- lsfit(pts$rt, pts$rtdev)
        rtdevsmo[[i]] <- rtime[[i_all]] * fit$coef[2] + fit$coef[1]
        ptsrange <- range(pts$rt)
        minidx <- rtime[[i_all]] < ptsrange[1]
        maxidx <- rtime[[i_all]] > ptsrange[2]
        rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
        rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
      }
      ## Finally applying the correction
      rtime_adj[[i_all]] <- rtime[[i_all]] - rtdevsmo[[i]]
    }
    ## Adjust the remaining samples.
    rtime_adj <- adjustRtimeSubset(rtime, rtime_adj, subset = subset,
                                   method = subsetAdjust)
    rtime_adj
  }

cluster_adjustRtime <- 
  function(object, param, msLevel = 1L) {
            if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
              object <- updateObject(object)
            startDate <- date()
            ## If param does contain a peakGroupsMatrix extract that one,
            ## otherwise generate it.
            if (nrow(peakGroupsMatrix(param))){
              pkGrpMat <- peakGroupsMatrix(param)
              } else {
              pkGrpMat <- cluster_adjustRtimePeakGroups(object, param = param)
              }
            message("Data must only contain one(!) MS level")
            res <- cluster_do_adjustRtime_peakGroups(
              peaks = chromPeaks(object, msLevel = msLevel),
              # OBS: I removed .update_feature_definitions, since we only use one MS level!
              peakIndex = featureDefinitions(object)$peakidx,
              rtime = rtime(object, bySample = TRUE),
              minFraction = minFraction(param),
              extraPeaks = extraPeaks(param),
              smooth = "loess",
              span = 0.8,
              family = "gaussian",
              peakGroupsMatrix = pkGrpMat,
              subset = 1:length(fileNames(object)),
              subsetAdjust = subsetAdjust(param)
            )
            ## Add the pkGrpMat that's being used to the param object.
            peakGroupsMatrix(param) <- pkGrpMat
            ## Dropping the peak groups but don't remove its process history
            ## step.
            ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
            object <- dropFeatureDefinitions(object)
            ## Add the results. adjustedRtime<- should also fix the retention
            ## times for the peaks! Want to keep also the latest alignment
            ## information
            adjustedRtime(object) <- res
            if (length(ph)) {
              object <- cluster_addProcessHistory(object, ph[[length(ph)]])
            }
            ## Add the process history step, get the msLevel from the peak
            ## detection step.
            ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
            xph <- XProcessHistory(param = param, date. = startDate,
                                   type. = .PROCSTEP.RTIME.CORRECTION,
                                   fileIndex = 1:length(fileNames(object)),
                                   msLevel = msLevel)
            object <- cluster_addProcessHistory(object, xph)
            validObject(object)
            object
          }


###################################### How to use ##############################################

msLevel <- 1
param <- PeakGroupsParam(smooth = "loess",
                         minFraction	= 0.5,
                         span = 0.25,
                         extraPeaks = 500,
                         family = "gaussian")

xset_aligned <- cluster_adjustRtime(object, param)
xset_aligned <- applyAdjustedRtime(xset_aligned)


save(xset_aligned, file = "./tmp/xset_aligned.Rdata")

