Density Peaks Advanced clustering in R
======================================

The `DPAclustR` package is a wrapper to the Python library `DPA` implementing the Density Peaks Advanced (DPA) clustering algorithm as introduced in the paper "Automatic topography of high-dimensional data sets by non-parametric Density Peak clustering", published on `M. d'Errico, E. Facco, A. Laio, A. Rodriguez, Information Sciences, Volume 560, June 2021, 476-492`_  (also available on `arXiv`_).

The package offers the following features:

* Advanced version of the `DP` clustering algorithm, published in the `Clustering by fast search and find of density peaks`_ paper, which includes an automatic search of cluster centers and assessment of statistical significance of the clusters.

Source files
------------

The source R codes are inside the ``R`` folder::

    .
    |-- ...
    |-- R/
    |   |-- Plotting.R                   # Visualization tools:
    |   |                                # dendrogram and network visualizations.
    |   |
    |   |-- AnalysisTools.R              # Normalized Mutual Information score,
    |   |                                # grid search for choosing of
    |   |                                # the external parameters
    |   |
    |   |-- DPAclustering.R              # Function running the Density Peak Advanced clustering.
    |
    |
    |-- ...

Getting started
===============

The source code of DPAclustR is on `github DPAclustR repository`_.


We suggest you to use `renv`_ as project environment for R.
The ``renv.lock`` file can be used to restore a
working status of the DPAclustR library by running the command ``renv::restore()``.
Alternatively, a new project environment can be initialized with the command ``renv::init()``.
See Installation section for more details.

The ``DPA`` Python package is required to use the ``DPAclustR::runDPAclustering`` function. Please see intruction on `github DPA repository`_ for installing it.
For more details on the ``runDPAclustering`` and its usage see Quickstart section below.

Installation
============

Assuming you already have R or RStudio installed on your machine,
run the following commands to install the SingleCellAnalysis package from github::

    renv::init() # suggested
    install.packages("devtools")
    devtools::install_github("mariaderrico/DPAclustR")


For development, you can clone the DPAclustR source code::

     git clone https://github.com/mariaderrico/DPAclustR.git

The ``renv.lock`` file in the package can be used to restore a working status of the
DPAclustR library by running the following commands::

     setwd("yourlocalpath/DPAclustR")
     renv::restore()
     library(DPAclustR)

Quickstart
----------

A use-case example of the analysis workflow is provided as jupyter notebook in ``Analysis_example.ipynb``. The same analysis workflow is provided as R script in ``Analysis_example.R`` that can be easily run within RStudio.

To run the clustering analysis using the DPA method, as described under section ``Analysis tools``
of the use-case example, the DPA Python package has to be loaded. Assuming the DPA package has been installed following the instructions available in the `github DPA repository`_, using virtualenv, the following commands have to be run::

    library(reticulate)
    use_virtualenv("path_toEnv/venvdpa/", required=TRUE)
    setwd("path_toDPA/DPA")
    DPA <- import_from_path("DPA", path="src/Pipeline/")

.. References

.. _`github DPAclustR repository`: https://github.com/mariaderrico/DPAclustR.git
.. _`renv`: https://blog.rstudio.com/2019/11/06/renv-project-environments-for-r/
.. _`github DPA repository`: https://github.com/mariaderrico/DPA.git
.. _`M. d'Errico, E. Facco, A. Laio, A. Rodriguez, Information Sciences, Volume 560, June 2021, 476-492`: https://www.sciencedirect.com/science/article/pii/S0020025521000116?dgcid=author
.. _`arXiv`: https://arxiv.org/abs/1802.10549v2
.. _`Computing the free energy without collective variables`: https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00916
.. _`Estimating the intrinsic dimension of datasets by a minimal neighborhood information`: https://export.arxiv.org/pdf/1803.06992
.. _`Clustering by fast search and find of density peaks`: http://science.sciencemag.org/content/344/6191/1492.full.pdf

