#' @title Create network graph for clustering results
#'
#' @name plot_networkD3
#'
#' @description Function defined in ``Plotting.R``. It visualizes the clustering results as
#' a network.
#'
#' @param ground_truth, list of ground_truth labels if available
#' @param clustering_labels, list of clustering labels from clustering
#' @param topography, list of connections between clusters as provided by the
#'                    DPA clustering or as returned by the function build_topography
#' @param popmin, default is 100. Clusters with population lower than popmin are excluded from the visualization.
#' @param ukwn, default is NULL. It can be set to the ground-truth label identified as uknown
#' @param opacity, default equals to 0.8. It set the opacity of nodes in the graph
#' @param ColourScale, set the colours of each node based on the most abundance ground_truth class.
#'
#' Note: required networkD3 package
#'
#' @importForm networkD3 forceNetwork
#' @export
plot_networkD3 <- function(ground_truth, clustering_labels, topography,
                           popmin = 100, ukwn = NULL, opacity = 0.8,
                           ColourScale = 'd3.scaleOrdinal()
                           .domain(["ADPr-Peptide", "randomPeptide"])
                           .range(["#FF6900", "#694489"]);'){
  # Remove uknown types in the ground_truth
  if(!is.null(ukwn)){
    mm <- ground_truth != ukwn
    ground_truth <- ground_truth[mm]
    clustering_labels <- clustering_labels[mm]
  }
  # Create the data structure for nodes in the graph
  # First column will be the clustering lable, then the ground_truth label
  # assigned according to majority rule.
  rule <- apply(table(ground_truth, clustering_labels), 2, which.max)
  labels <- names(rule)
  rule <- sort(unique(ground_truth))[rule]
  pop <- table(clustering_labels)
  nodes <- data.frame(rule, pop)
  colnames(nodes) <- c("group", "name", "size")

  # Filter out cluster with population lower than pop:
  listnames <- nodes[nodes$size>popmin,"name"]
  listukwn <- nodes[nodes$group==ukwn,"name"]
  # networkd3 requires consecutive labels: rename after pop cuts or on case of missing labels
  nodes <- nodes[nodes$size>popmin,]
  rownames(nodes) <- 1:nrow(nodes)
  nodes$new_name <- rownames(nodes)

  # Create the data structure for edges in the graph:
  # Remove links which has zero value, for visualization
  topography <- topography[topography$value>0,]
  # Remove links to nodes which have population lower than popmin
  topography <- topography[topography$source %in% listnames & topography$target %in% listnames
                           & !(topography$source %in% listukwn)
                           & !(topography$target %in% listukwn),]
  #topography$value <- as.numeric(topography$value)*0.01
  topography$source <- c(lapply(topography$source, function(x) nodes[nodes$name==x,"new_name"]))
  topography$target <- c(lapply(topography$target, function(x) nodes[nodes$name==x,"new_name"]))

  # networkd3 seems taking labels from zero (?) : shift labels by one:
  nodes$name <- as.numeric(nodes$name)-1
  topography$source <- as.numeric(topography$source)-1
  topography$target <- as.numeric(topography$target)-1

  networkD3::forceNetwork(Links = topography, Nodes = nodes,
                          Source = "source", Target = "target",
                          Value = "value", NodeID = "new_name",
                          Nodesize = "size",
                          Group = "group", opacity = opacity,
                          linkDistance = networkD3::JS("function(d) { return 10*d.value; }"),
                          linkWidth = networkD3::JS("function(d) { return d.value/5; }"),
                          zoom = FALSE, fontSize = 20,
                          colourScale = networkD3::JS(ColourScale), opacityNoHover = 1)

}

#' @title Create dendrogram for clustering results
#'
#' @name plot_dendrogram
#'
#' @description Function defined in ``Plotting.R``. It visualizes the clustering results as
#' a dendrogram.
#'
#' @param clustering_labels, list of clustering labels from clustering
#' @param topography, list of connections between clusters as provided by the
#'                    DPA clustering or as returned by the function build_topography
#' @param maxD, normalization value
#' @param popmin, default is 0. Clusters with population lower than popmin are excluded from the visualization.
#' @param method, default "average". It it a hclust parameter, see help(hclust) for other options.
#'
#' Note: required hclust package
#'
#' @importFrom tidyr spread
#' @export

plot_dendrogram <- function(clustering_labels, topography, maxD, popmin=0, method="average"){
  # Check that the topography is not empty
  if(nrow(topography)==0){
    stop("No dendrogram can be created: DPA returned a single cluster and an empty topography")
  }
  message("Start creating dendrogram")

  # Create the data structure for leaves in the tree
  pop <- table(clustering_labels)
  nodes <- data.frame(pop)
  colnames(nodes) <- c("name", "size")

  # Filter out cluster with population lower than pop:
  listnames <- nodes[nodes$size>popmin,"name"]
  # networkd3 requires consecutive labels: rename after pop cuts or on case of missing labels
  nodes <- nodes[nodes$size>popmin,]
  rownames(nodes) <- 1:nrow(nodes)
  nodes$new_name <- rownames(nodes)

  ### Create the data structure for the dendrogram
  # convert similarities into distance normalizing to the maximum density in the data set
  topography$temp_value <- lapply(topography$value, function(x){ (maxD-x)/maxD})

  # Remove leaves which have population lower than popmin
  topography <- topography[topography$source %in% listnames & topography$target %in% listnames,]
  topography$source <- c(lapply(topography$source, function(x) {nodes[nodes$name==x,"new_name"]}))
  topography$target <- c(lapply(topography$target, function(x) {nodes[nodes$name==x,"new_name"]}))

  # create the similarity matrix
  data_wide <- tidyr::spread(topography[,c("source","target","temp_value")], target, temp_value, fill=NA)
  data_wide <- as.data.frame(t(data_wide[, !(names(data_wide) %in% "source")]))
  Nclus <- nrow(nodes)-1
  names(data_wide) <- 1:Nclus
  # add first row and last column to adjust to NclusxNclus similarity matrix format
  data_wide <- as.matrix(data_wide)
  data_wide <- rbind(c(rep(0,Nclus)), data_wide)
  data_wide <- cbind(data_wide, c(rep(0,Nclus)))
  data_wide <- as.data.frame(data_wide)
  names(data_wide) <- 1:(Nclus+1)
  rownames(data_wide) <- 1:(Nclus+1)
  data_wide <- as.matrix(data_wide)
  # create a lower diagonal matrix with zeros elsewhere
  data_wide[upper.tri(data_wide,diag=TRUE)] <- 0

  # Build the dist object from the similarity matrix
  d<-as.dist(data_wide, diag=TRUE)
  hc <- hclust(d, method=method)

  # Set dendrogram labels to original clustering labels
  #plot(x = hc, labels = nodes$name, cex = 0.5)
  # Return the dendrogram for plotting
  message("Dendrogram done.")
  return(hc)

}
