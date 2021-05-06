#' @title Apply majority rule
#'
#' @name majority_list
#'
#' @description Each point within a cluster is assigned to the ground_truth class of
#' highest frequency in the cluster. Function defined in ``AnalysisTools.R``.
#'
#' @param ground_truth, list of ground truth labels
#' @param clustering_labels, output from a clustering algorithm
#'
#' @export
majority_list <- function(ground_truth, clustering_labels){
  # To each clusters label assign the class in the ground truth having more
  # occurrences in the clusters
  rule <- apply(table(ground_truth, clustering_labels), 2, which.max)
  labels <- names(rule)
  # The next line extend the computation for cases where
  # the labels in the ground_truth are not consecutives
  rule <- sort(unique(ground_truth))[rule]
  mm <- match(clustering_labels, labels)
  rule[mm]
}

#' @title Compute Normalized Mutual Information score
#'
#' @name getscore_NMI
#'
#' @description Function defined in ``AnalysisTools.R``.
#'
#' @param ground_truth, list of ground truth labels
#' @param labels, output from a clustering algorithm after majority assignment
#'
#' Note: required Python package scikitlearn
#'
#' @export
getscore_NMI <- function(ground_truth, labels){
  nmi <- skl$normalized_mutual_info_score(ground_truth,labels)
  nmi
}

#' @title Compute grid search for external parameters
#'
#' @name grid_search
#'
#' @description Function defined in ``AnalysisTools.R``.
#'
#' @param method, clustering method 'DPAclustering'
#' @param min, minimum value for the external parameter to explore
#' @param max, maximum value for the external parameter to explore
#' @param step, increment of the sequence
#' @param dataset, input data set
#' @param ground_truth, list of ground truth labels
#' @param ukwn, default NULL. If `ukwn` is set to the label in the ground truth corresponding to
#'              an NA or UKNOWN classification the UKNOWN cell will be excluded by the score
#'              calculation and performance evaluation. Otherwise set `uknw`=NULL to include
#'              those cells as well.
#' @param fileprefix, default NULL. If `fileprefix` is provided it will read clustering results
#'              from file, otherwise the clustering algorithm will be run again.
#' @param dir, set the path to the directory where to save clustering outputs to.
#'
#' @export
grid_search <- function(method, min, max, step, dataset, ground_truth, ukwn=NULL,
                        fileprefix=NULL, dir){
  modellist <- list()
  if(!is.null(ukwn)){
    mm <- ground_truth != ukwn
    ground_truth <- ground_truth[mm]
  }
  switch(
    EXPR = method,
    "DPAclustering" = {
      for (param in seq(min, max, step)){
        if (is.null(fileprefix)){
          DPAobj <- DPAclustR::runDPAclustering(dataset, Z=param, dir=dir)
          clustering_labels <- DPAobj@labels
        } else {
          clustering_labels <- read.csv(paste0(dir,fileprefix,toString(param)))[,1]
        }
        if(!is.null(ukwn)){
          clustering_labels <- clustering_labels[mm]
        }
        mjr <- majority_list(ground_truth, clustering_labels)
        nmi <- getscore_NMI(ground_truth, mjr)
        key <- toString(param)
        modellist[[key]] <- nmi
      }
    },
    stop("Invalid method. Please choose 'DPAclustering' ")
  )
  results <- (modellist)
}
