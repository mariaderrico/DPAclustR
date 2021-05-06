context("Applying DPA clustering algorithm")

# Only test if DPA is available as python package


skip_if_no_dpa <- function() {
  have_dpa <- reticulate::py_module_available("DPA")
  if (!have_dpa)
    skip("DPA not available for testing")
}

# this "try" block is necessary because:
# a system that does not have python at all will generate a warning,
# which can generate a NOTE during R CMD check (winbuilder)
has.dpa = FALSE
try({
  has.dpa = reticulate::py_module_available("DPA")
}, silent=TRUE)



#if (has.dpa) { 
test_that("the DPA clustering run properly", {
  skip_if_no_dpa()
  setwd("./benchmarkdata")
  df_norm <- read.csv("df_norm.csv")
  clustering_bmk <- read.csv("DPA_results/cell_clustering_DPA_Z1", header=TRUE)$x

  # dir=NULL if you don't want to save the results
  DPAresult <- runDPAclustering(df_norm, Z=1, dir=NULL)
  cell_clustering_DPA <- DPAresult@labels
  expect_equal(all.equal(as.numeric(cell_clustering_DPA),as.numeric(clustering_bmk)), TRUE)  

})

