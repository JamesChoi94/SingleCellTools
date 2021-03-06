% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StandardLR.R
\name{StandardLR}
\alias{StandardLR}
\title{Ligand-Receptor scoring - mean of mean expressions}
\usage{
StandardLR(
  object,
  ref.path = NULL,
  lr.ref = NULL,
  split.by = NULL,
  min.pct = 0.1,
  assay = "RNA",
  slot = "data",
  resample = 1000,
  adjust.pval = FALSE
)
}
\arguments{
\item{object}{Seurat object containing RNA expression data.}

\item{ref.path}{Character string of path to ligand-receptor pair reference
list. LR pair reference must contain a column labeled "Pair.Name" with values
for each LR pair. If \code{NULL}, \code{lr.ref} must be provided.}

\item{lr.ref}{Data.frame of ligand-receptor pair reference. Must contain a
column labeled "Pair.Name" with values for each LR pair. If \code{NULL},
\code{ref.path} must be provided.}

\item{split.by}{Character string of a column name in
\code{slot(tmp,'meta.data')} by which to split cells before calculating LR
scores e.g. across multiple conditions or time-points.}

\item{min.pct}{Numeric minimum percentage a ligand or receptor gene must be
expressed in any cell cluster to be retained for LR scoring. Note: scores
will still be calculated between pairs where one cluster expresses at 15\% and
the other at 0\%. These scores should be filtered out at visualization
(see: \code{\link{PlotLR}}).}

\item{assay}{Character string to select which assay slot of Seurat object to
use.}

\item{resample}{Numeric number of times to sample cells for permutation test.}

\item{adjust.pval}{Logical determining whether to perform max-T p-value
adjustment.}
}
\value{
A data.frame containing results of the standard ligand-receptor
analysis. The columns of the data.frame as as follow:
\itemize{
\item{'Pair_name'} : Name of ligand gene followed and separated by name of
receptor gene.
\item{'Score'} : Average of mean ligand expression in ligand cluster and
mean receptor expression in receptor cluster.
\item{'pval'} : Estimated p-value using a permutation test. See description
for more method details.
\item{'adj_pval'} : Adjusted p-value using a max-T adjustment. See https://statistics.berkeley.edu/sites/default/files/tech-reports/633.pdf equation 11 for more details.
\item{'Ligand_cell'} : Ligand-expressing cell cluster.
\item{'Receptor_cell'} : Receptor-expressing cell cluster.
}
}
\description{
Ligand-Receptor scoring - mean of mean expressions
}
