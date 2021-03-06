#' Ligand-Receptor scoring - mean of mean expressions
#'
#'
#'
#'
#' @import dplyr
#' @import tibble
#' @importFrom rlang "!!"
#'
#' @param object Seurat object containing RNA expression data.
#' @param ref.path Character string of path to ligand-receptor pair reference
#' list. LR pair reference must contain a column labeled "Pair.Name" with values
#' for each LR pair. If \code{NULL}, \code{lr.ref} must be provided.
#' @param lr.ref Data.frame of ligand-receptor pair reference. Must contain a
#' column labeled "Pair.Name" with values for each LR pair. If \code{NULL},
#' \code{ref.path} must be provided.
#' @param split.by Character string of a column name in
#' \code{slot(tmp,'meta.data')} by which to split cells before calculating LR
#' scores e.g. across multiple conditions or time-points.
#' @param min.pct Numeric minimum percentage a ligand or receptor gene must be
#' expressed in any cell cluster to be retained for LR scoring. Note: scores
#' will still be calculated between pairs where one cluster expresses at 15% and
#'  the other at 0%. These scores should be filtered out at visualization
#' (see: \code{\link{PlotLR}}).
#' @param assay Character string to select which assay slot of Seurat object to
#' use.
#' @param resample Numeric number of times to sample cells for permutation test.
#' @param adjust.pval Logical determining whether to perform max-T p-value
#' adjustment.
#'
#' @return A data.frame containing results of the standard ligand-receptor
#' analysis. The columns of the data.frame as as follow:
#' \itemize{
#'  \item{'Pair_name'} : Name of ligand gene followed and separated by name of
#'  receptor gene.
#'  \item{'Score'} : Average of mean ligand expression in ligand cluster and
#'  mean receptor expression in receptor cluster.
#'  \item{'pval'} : Estimated p-value using a permutation test. See description
#'  for more method details.
#'  \item{'adj_pval'} : Adjusted p-value using a max-T adjustment. See https://statistics.berkeley.edu/sites/default/files/tech-reports/633.pdf equation 11 for more details.
#'  \item{'Ligand_cell'} : Ligand-expressing cell cluster.
#'  \item{'Receptor_cell'} : Receptor-expressing cell cluster.
#' }
#'
#' @export
#'


StandardLR <- function(
  object,
  ref.path = NULL,
  lr.ref = NULL,
  split.by = NULL,
  min.pct = 0.1,
  assay = "RNA",
  slot = "data",
  resample = 1000,
  adjust.pval = FALSE
) {

  # 'object' input check
  if (class(object) != 'Seurat') {
    stop('\"object\" must be of class Seurat')
  }
  # 'ref.path' or 'lr.ref' check
  if (is.null(ref.path) && is.null(lr.ref)) {
    stop(strwrap('Must provide either \"ref.path\" (path of LR reference) or
    \"lr.ref\" (data.frame of imported LR reference).', prefix = ' '))
  }

  # Load reference csv.
  if (!is.null(ref.path)) {
    lr.ref <- read.csv(file = ref.path, stringsAsFactors = FALSE)
  }

  # Check for Pair.Name column. Error if not present.
  if (!any(colnames(lr.ref) == 'Pair.Name')) {
    stop(strwrap('Ligand-receptor reference list requires a column titled
    \"Pair.Name\". Entries in \"Pair.Name\" should be formatted as: [Ligand
    gene]_[Receptor gene]. E.g. Apoe_Lrp1"', prefix = ' '))
  } else {
    pair_column <- which(colnames(lr.ref) == 'Pair.Name')
  }
  message('LR reference has ', nrow(lr.ref), ' rows (LR pairs).')

  # Split Pair.Name into ligand and receptor names
  tmp <- strsplit(x = lr.ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)

  # Retain all ligand-receptor pairs where both gene names are present in the
  # gene expression dataset
  all_genes <- rownames(slot(object = object[[assay]], 'counts'))
  lr.ref <- lr.ref[ligand_names %in% all_genes & receptor_names %in% all_genes,]
  if (nrow(lr.ref) == 0) {
    stop('No LR pairs were detected in provided Seurat data.')
  }

  # Get ligands and receptors that were detected in data
  tmp <- strsplit(x = lr.ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)

  # Message regarding use of seurat identities
  tmp_idents <- paste(unique(slot(object = object, name = 'active.ident')),
                      collapse = ', ')
  message(paste0('Using active identities: ', tmp_idents))

  # New vector of all genes to retrieve data
  retrieve_genes <- union(ligand_names, receptor_names)
  active_idents <- slot(
    object = object,
    name = 'active.ident'
  )

  # Extract data
  Seurat::DefaultAssay(object) <- assay
  lr_data <- Seurat::FetchData(
    object = object,
    vars = c(split.by, retrieve_genes),
    slot = slot
  )
  active_idents <- slot(
    object = object,
    name = 'active.ident'
  )
  lr_data <- cbind(active_idents, lr_data)
  message('Using expression values for ', ncol(lr_data), ' genes across ',
          nrow(lr_data), ' cells.')

  # Calculate average/percent expression for each gene, by "split.by" if
  # provided.
  exp_avg <- lr_data %>%
    group_by(active_idents, .add = TRUE)
  exp_pct <- lr_data %>%
    group_by(active_idents, .add = TRUE)
  if (!is.null(split.by)) {
    split.by_name <- as.name(split.by)
    exp_avg <- exp_avg %>%
      group_by(!!split.by_name, .add = TRUE)
    exp_pct <- exp_pct %>%
      group_by(!!split.by_name, .add = TRUE)
  }
  exp_avg <- exp_avg %>%
    summarise(across(where(is.numeric), .fns = mean))
  exp_pct <- exp_pct %>%
    summarise(across(where(is.numeric),
                     .fns = function(x) round(mean(x > 0), 3)))

  # Determine which genes meet minimum percent detection threshold
  minpct_genes <- sapply(
    X = exp_pct[sapply(exp_pct, is.numeric)],
    FUN = function(x) any(x > 0)
  )
  minpct_genes <- names(minpct_genes)[minpct_genes]
  lr.ref <- lr.ref[ligand_names %in% minpct_genes &
                     receptor_names %in% minpct_genes,]
  lr_genes <- sort(unique(minpct_genes))

  # Table of cell counts, further split if "split.by" provided.
  if (!is.null(split.by)) {
    cell_counts <- table(
      slot(object = object, name = 'active.ident'),
      slot(object = object, name = 'meta.data')[[split.by]]
    )
  } else {
    cell_counts <- table(slot(object = object, name = 'active.ident'))
  }

  # Cell-level expression matrix for genes in LR reference
  exp_mat <- as.matrix(lr_data[sapply(lr_data, is.numeric)])

  # Extract ligand/receptor gene names + expression matrices
  tmp <- strsplit(x = lr.ref[['Pair.Name']], split = '_')
  ligand_names <- sapply(X = tmp, FUN = `[`, 1)
  receptor_names <- sapply(X = tmp, FUN = `[`, 2)
  exp_mat_ligands <- exp_mat[,match(ligand_names, colnames(exp_mat))]
  exp_mat_receptors <- exp_mat[,match(receptor_names, colnames(exp_mat))]


  # Data.frame with results and null distribution score values
  if (class(lr_data[['active_idents']]) != 'factor') {
    stop('Seurat identities must be of class factor.')
  }
  if (!is.null(split.by)) {
    var_set <- expand.grid(levels(lr_data[['active_idents']]),
                           levels(lr_data[['active_idents']]),
                           levels(lr_data[[split.by]]))
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell', 'split.by')
  } else {
    var_set <- expand.grid(levels(lr_data[['active_idents']]),
                           levels(lr_data[['active_idents']]))
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell')
  }

  null_scores <- vector(mode = 'list', length = nrow(var_set))
  names(null_scores) <- apply(
    X = var_set,
    MARGIN = 1,
    FUN = paste,
    collapse = '_'
  )

  # Use maxT method for multiple hypotheses p-value adjustment
  # (Benjamini-Hochberg appears too hard for 1000x permutations)
  if (adjust.pval) {
    maxT_scores <- matrix(0, nrow = resample, ncol = nrow(lr.ref))
    colnames(maxT_scores) <- paste(colnames(exp_mat_ligands),
                                   colnames(exp_mat_receptors),
                                   sep = '_')
  }

  # progress bar
  pb = txtProgressBar(min = 0, max = nrow(var_set), initial = 0, style = 3)

  # calculate null distribution of randomly permuted ligand-receptor score values
  for(i in 1:nrow(var_set)) {
    # these are factors, but cell_counts is sorted by level
    cellx <- as.character(var_set[['Ligand_cell']][i])
    celly <- as.character(var_set[['Receptor_cell']][i])
    if (!is.null(split.by)) {
      split.byz <- as.character(var_set[['split.by']][i])
      countx <- cell_counts[cellx, split.byz]
      county <- cell_counts[celly, split.byz]
    } else {
      countx <- cell_counts[cellx]
      county <- cell_counts[celly]
    }
    if (countx == 0 | county == 0) {
      next()
    }

    # Util function to randomly select n_c cells (n = #, c = cell-type) and
    # extract element positions.
    cell_sample <- function(x) {
      tmp <- rep(FALSE, nrow(exp_mat))
      tmp[sample(nrow(exp_mat), size = x, replace = FALSE)] <- TRUE
      return(tmp)
    }
    null_index_x <- t(replicate(n = resample, expr = cell_sample(countx)))
    null_index_y <- t(replicate(n = resample, expr = cell_sample(county)))

    # Calculate average expression of all ligands and receptors for the n_c
    # cells
    null_avg_lig <- (null_index_x %*% exp_mat_ligands) / countx
    null_avg_rec <- (null_index_y %*% exp_mat_receptors) / county

    # Calculate null LR scores
    tmp_scores <- 1/2 * (null_avg_lig + null_avg_rec)
    colnames(tmp_scores) <- paste(colnames(null_avg_lig),
                                  colnames(null_avg_rec),
                                  sep = '_')

    # For max-T method, take element-wise maximum values (this iterates
    # nrow(var_set) times)
    null_scores[[i]] <- tmp_scores
    if (adjust.pval) {
      change_score <- tmp_scores > maxT_scores
      maxT_scores[change_score] <- tmp_scores[change_score]
    }

    # progress bar
    setTxtProgressBar(pb, i)
  }

  # Calculate average expression matrix
  exp_names <- exp_avg[['active_idents']]
  if (!is.null(split.by)) {
    exp_names <- paste(exp_names, exp_avg[[split.by]], sep = '_')
  }

  # Get indices from expression matrix for ligands/receptors
  index_x <- var_set[['Ligand_cell']]
  index_y <- var_set[['Receptor_cell']]
  if (!is.null(split.by)) {
    index_x <- paste(index_x, var_set[['split.by']], sep = '_')
    index_y <- paste(index_y, var_set[['split.by']], sep = '_')
  }
  index_x <- match(x = index_x, table = exp_names)
  index_y <- match(x = index_y, table = exp_names)
  index_l <- match(ligand_names, colnames(exp_avg))
  index_r <- match(receptor_names, colnames(exp_avg))
  exp_avg_lig <- as.matrix(exp_avg[index_x, index_l])
  exp_avg_rec <- as.matrix(exp_avg[index_y, index_r])

  # Calculate LR scores (for data, not nulls)
  exp_scores <- 1/2 * (exp_avg_lig + exp_avg_rec)

  # NOTE: Percent threshold step is removed here and left to the visualization
  # step for removal of data points.
  # # Determine which cells do not have at least 10% expression of their gene
  # exp_pct_l <- exp_pct[index_x, index_l] < 0.1
  # exp_pct_r <- exp_pct[index_y, index_r] < 0.1
  # exp_pct_lr <- exp_pct_l | exp_pct_r

  # Replace values that don't meet 10% threshold
  # exp_scores[exp_pct_lr] <- NA
  rownames(exp_scores) <- apply(
    X = var_set, MARGIN = 1,
    FUN = paste,
    collapse = '_'
  )
  colnames(exp_scores) <- paste(
    colnames(exp_avg_lig),
    colnames(exp_avg_rec),
    sep = '_'
  )

  # Calculate p-values from ecdf using null scores for each LR-pair.
  # NOTE: A true permutation test will test all possible permutations of cell-
  # sampling. Since this is computational impractical, we estimate with 1000
  # permutations. Thus if LR-scores ("effects") are large, then p-values can be
  # zero.
  message('\nCalculating p-values...')
  pvals <- matrix(NA, nrow = nrow(var_set), ncol = ncol(exp_scores))

  # For raw p-values:
  for(i in 1:nrow(pvals)) {
    tmp_nulls <- null_scores[[i]]
    tmp_scores <- exp_scores[i,]
    for(j in 1:length(tmp_scores)) {
      # Skip when cell counts were zero
      if (is.null(tmp_nulls)) {
        next()
      } else {
        pvals[i,j] <- 1-ecdf(tmp_nulls[,j])(tmp_scores[j])
      }
    }
  }
  colnames(pvals) <- paste(
    colnames(exp_avg_lig),
    colnames(exp_avg_rec),
    sep = '_'
  )
  rownames(pvals) <- apply(
    X = var_set,
    MARGIN = 1,
    FUN = paste,
    collapse = '_'
  )

  # For adjusted p-values:
  if (adjust.pval) {
    adj_pvals <- matrix(1, nrow = nrow(var_set), ncol = nrow(lr.ref))
    for(i in 1:ncol(adj_pvals)) {
      adj_pvals[,i] <- sapply(X = exp_scores[,i],
                              nulls = maxT_scores[,i],
                              FUN = function(x, nulls) {1-ecdf(nulls)(x)})
    }
    colnames(adj_pvals) <- paste(
      colnames(exp_avg_lig),
      colnames(exp_avg_rec),
      sep = '_')
    rownames(adj_pvals) <- apply(
      X = var_set,
      MARGIN = 1,
      FUN = paste,
      collapse = '_'
    )
  }

  # Long-form data.table of results
  tmp <- expand.grid(rownames(exp_scores),
                     colnames(exp_scores),
                     stringsAsFactors = FALSE)
  tmp_idents <- strsplit(x = tmp[[1]], split = '_')
  tmp_pairs <- strsplit(x = tmp[[2]], split = '_')
  ligand_cell <- sapply(X = tmp_idents, FUN = `[[`, 1)
  receptor_cell <- sapply(X = tmp_idents, FUN = `[[`, 2)
  ligand <- sapply(X = tmp_pairs, FUN = `[[`, 1)
  receptor <- sapply(X = tmp_pairs, FUN = `[[`, 2)
  cell_pair <- paste(ligand_cell, receptor_cell, sep = '_')
  lr_pair <- paste(ligand, receptor, sep = '_')
  tmp_scores <- exp_scores %>% reshape2::melt() %>% .[['value']]
  tmp_avg_l <- c(exp_avg_lig)
  tmp_avg_r <- c(exp_avg_rec)
  tmp_pvals <- pvals %>% reshape2::melt() %>% .[['value']]
  if (adjust.pval) {
    tmp_adj_pvals <- adj_pvals %>% reshape2::melt() %>% .[['value']]
  }
  tmp_pct <- exp_pct %>% as.matrix()
  tmp_pct_l <- c(tmp_pct[index_x, index_l])
  tmp_pct_r <- c(tmp_pct[index_y, index_r])
  split_var <- NA
  if (!is.null(split.by)) {
    split_var <- sapply(X = tmp_idents, FUN = `[[`, 3)
    split_var <- factor(x = split_var, levels = levels(lr_data[[split.by]]))
    cell_prop <- round(prop.table(cell_counts, margin = 1), 3)*100
    tmp_count_l <- mapply(
      FUN = function(x,y) cell_prop[x,y],
      ligand_cell,
      split_var
    )
    tmp_count_r <- mapply(
      FUN = function(x,y) cell_prop[x,y],
      receptor_cell,
      split_var
    )
  }

  # Final result compilation
  results <- data.frame(
    'Pair_name' = lr_pair,
    'Score' = tmp_scores,
    'pval' = tmp_pvals,
    'adj_pval' = ifelse(test = exists('tmp_adj_pvals'),
                        yes = tmp_adj_pvals,
                        no = NA),
    'Ligand_cell' = ligand_cell,
    'Receptor_cell' = receptor_cell,
    'split.by' = split_var,
    'Ligand' = ligand,
    'Receptor' = receptor,
    'Ligand_avgExp' = tmp_avg_l,
    'Receptor_avgExp' = tmp_avg_r,
    'Ligand_pct' = tmp_pct_l,
    'Receptor_pct' = tmp_pct_r,
    'Cell_pair' = cell_pair,
    'LR_pair' = lr_pair,
    stringsAsFactors = FALSE
  )
  if (!is.null(split.by)) {
    results[['Ligand_cell_pct']] <- tmp_count_l
    results[['Receptor_cell_pct']] <- tmp_count_r
  }
  message('Done!')
  return(results)
}
