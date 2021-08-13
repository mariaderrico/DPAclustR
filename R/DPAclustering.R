# For interfacing with python and "DPA" clustering
#' @importFrom reticulate py_module_available import_from_path
NULL


#' @title Run Density Peak Advanced clustering
#'
#' @name runDPAclustering
#'
#' @description Function defined in ``DPAclustering.R``.
#'
#' @param data, input data after transformation or normalization
#' @param Z, the number of standard deviations fixing the level of
#'            statistical confidence at which one decides to consider
#'            a cluster meaningful. Default value is set to 1.
#' @param metric, string, valid metrics are "precomputed", "euclidean".
#'            If metric is "precomputed", data is assumed to be a distance matrix.
#'            Default is "euclidean".
#' @param maxpop, clusters with population above maxpop are discarded.
#'            If maxpop is euqal to NULL, no filter is applied.
#' @param dir, path to the directory where to save clustering results. Default is NULL.
#'
#' @export
runDPAclustering <- function(data, Z=1, metric="euclidean", maxpop=NULL, dir=NULL){
  message(paste("Running DPA clustering for Z=",Z))
  # ReadInput requires a data.matrix input
  if (inherits(data, "data.frame"))
    data <- data.matrix(data)
  message("Input data loaded.")
  est <- DPA$DensityPeakAdvanced(Z=Z, metric=metric)
  res_dpa <- est$fit(data)
  message("DPA clustering done.")
  labels <- res_dpa$labels_
  # Create topography data frame
  b_topography <- as.data.frame(as.matrix(res_dpa$topography_))
  L <- length(unlist(b_topography$V1))
  if(L>0){
    topography <- data.frame(unlist(b_topography$V1)[seq(1,L,4)], unlist(b_topography$V1)[seq(2,L+1,4)],
                             unlist(b_topography$V1)[seq(3,L+1,4)], unlist(b_topography$V1)[seq(4,L+1,4)])
    colnames(topography) <- c("source","target","border","border_err")
    # Increase by one the labels so they start from 1
    topography$source <- topography$source+1
    topography$target <- topography$target+1
    # Create a distances measure between clusters as max(all_densities)-border_density
    topography$value <- c(lapply(topography$border, function(x) if (x>0) {max(res_dpa$densities_)-x} else{x}))
  }else{
    topography <- data.frame()
  }
  # Select clusters with population higher than maxpop
  if(!is.null(maxpop)){
    C <- 1:max(labels)
    subset <- which(as.numeric(table(labels))>maxpop)
    C_temp <- rapply(list(C), function(x) ifelse(x%in%subset, C[x]-1, max(subset)), how = "replace")
    C_subset <- which(unlist(C_temp)<max(subset))-1
    labels_maxpop <- unlist(rapply(list(labels),function(x) ifelse(!(x%in%C_subset),max(C_subset)+1,x), how = "replace"))
  }
  else {
    labels_maxpop <- labels
  }
  # Labels must be starting from 1
  # labels_maxpop+1
  if(!is.null(dir)){
    write.csv(labels_maxpop+1,paste0(dir,"cell_clustering_DPA_Z",toString(Z)), row.names = FALSE)
    write.csv(res_dpa$densities_,paste0(dir,"densities_DPA_Z",toString(Z)), row.names = FALSE)
    write.csv(res_dpa$centers_,paste0(dir,"centers_DPA_Z",toString(Z)), row.names = FALSE)
    if(L>0){
      topography$value <- as.numeric(topography$value)
      topography$source <- as.numeric(topography$source)
      topography$target <- as.numeric(topography$target)
      topography$border <- as.numeric(topography$border)
      topography$border_err <- as.numeric(topography$border_err)
      write.csv(topography,paste0(dir,"topography_DPA_Z",toString(Z)), row.names = FALSE)
    }
  }
  #return(new("DPAresults",
  return(list(labels=labels_maxpop+1,
             density=res_dpa$densities_,
             centers=res_dpa$centers_,
             topography=topography))
}


#setClass(Class="DPAresults",
#         representation(
#           labels="array",
#           density="numeric",
#           centers="integer",
#           topography="list"
#         )
#)

