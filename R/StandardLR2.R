#' Ligand-Receptor scoring - mean of mean expressions
#'
#'
#'
#'
#' @import dplyr
#' @import tibble
#' @import data.table
#' @importFrom rlang "!!"
#' @importFrom BiocParallel bpmapply
#'
#' @param seurat.object Seurat object containing RNA expression data.
#' @param lr.ref Data.frame of ligand-receptor pair reference. Must contain a
#' column labeled "Pair.Name" with values for each LR pair. If \code{NULL},
#' \code{lr.ref.path} must be provided.
#' @param lr.ref.path Character string of path to ligand-receptor pair reference
#' list. LR pair reference must contain a column labeled "Pair.Name" with values
#' for each LR pair. If \code{NULL}, \code{lr.ref} must be provided.
#' @param split.by Character string of a column name in
#' \code{slot(tmp,'meta.data')} by which to split cells before calculating LR
#' scores e.g. across multiple conditions or time-points. If \code{NULL}, no
#' splitting is done.
#' @param min.pct Numeric minimum percentage a ligand or receptor gene must be
#' expressed in any cell cluster to be retained for LR scoring. Note: scores
#' will still be calculated between pairs where one cluster expresses at 15% and
#'  the other at 0%. These scores should be filtered out at visualization
#' (see: \code{\link{PlotLR}}).
#' @param assay Character string to select which assay slot of Seurat object to
#' use. Default is "RNA" slot.
#' @param slot Character string to select which slot of the provided to assay
#' to extract values (e.g. counts, data, scaled.data). Default is "data" slot.
#' @param resample Numeric number of times to sample cells for permutation test.
#' Default is 1000.
#' @param adjust.pval Logical determining whether to perform max-T p-value
#' adjustment. Default is TRUE.
#' @param ligand.cells Character vector of cell-types to scared as the ligand
#' expressing cell-type e.g. c("Macrophage","Microglia"). If \code{NULL}, use
#' all cell-types.
#' @param receptor.cells Character vector of cell-types to scared as the
#' receptor expressing cell-type e.g. c("Macrophage","Microglia"). If
#' \code{NULL}, use all cell-types.
#' @param split.vars Character vector of subset of "split.by" variable e.g. if
#' split.by = "time", then split.vars = c("Uninjured","1dpi'). If \code{NULL},
#' use all values of "split.by".
#' @param BPPARAM Instance of \code{BiocParallelParam} class for back-ends.
#' Default is first entry of \code{BiocParallel::registered()}.
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

StandardLR2 <- function(
  seurat.object,
  lr.ref = NULL,
  lr.ref.path = NULL,
  split.by = NULL,
  min.pct = 0.1,
  assay = "RNA",
  slot = "data",
  resample = 1000,
  adjust.pval = FALSE,
  ligand.cells = NULL,
  receptor.cells = NULL,
  split.vars = NULL,
  BPPARAM = bpparam()
) {

  # 'seurat.object' input check
  if (class(seurat.object) != 'Seurat') {
    stop('\"seurat.object\" must be of class Seurat.')
  }

  # Check that provided assay is in seurat.object.
  if (!(assay %in% names(slot(object = seurat.object, name = 'assays')))) {
    stop('assay \"', assay, '\" cannot be found in any slots of Seurat object.')
  }

  # 'lr.ref.path' or 'lr.ref' check
  if (is.null(lr.ref.path) && is.null(lr.ref)) {
    stop(strwrap('Must provide either \"lr.ref.path\" (path of LR reference) or
    \"lr.ref\" (data.frame of imported LR reference).', prefix = ' '))
  }

  # Load reference csv.
  if (!is.null(lr.ref)) {
    lr_ref <- lr.ref
  }
  if(!is.null(lr.ref.path)) {
    lr_ref <- read.csv(file = lr.ref.path, stringsAsFactors = FALSE)
  }

  # Check for Pair.Name column. Error if not present.
  if (!any(colnames(lr_ref) == 'Pair.Name')) {
    stop(strwrap('Ligand-receptor reference list requires a column titled
    \"Pair.Name\". Entries in \"Pair.Name\" should be formatted as: [Ligand
    gene]_[Receptor gene]. E.g. Apoe_Lrp1"', prefix = ' '))
  } else {
    message(paste0('LR reference has ', nrow(lr_ref), ' rows (LR pairs).'))
  }

  # Message regarding use of seurat identities
  active_idents <- as.character(x = slot(
    object = seurat.object,
    name = 'active.ident'
  ))
  if (any(grepl(' |\\.|\\*|\\&', x = unique(active_idents)))) {
    stop(strwrap('Active identities contain special characters (such as
                    spaces). Replace with "_".', prefix = ' '))
  }
  tmp_idents <- paste(unique(active_idents), collapse = ', ')
  message(paste0('Using active identities: ', tmp_idents))
  rm(tmp_idents)

  # Check if ligand_cells and receptor_cells are in active.idents.
  if (!is.null(ligand.cells)) {
    if (!all(ligand.cells %in% active_idents)) {
      stop('Not all provided \"ligand.cells\" are in active.idents.')
    }
  }
  if (!is.null(receptor.cells)) {
    if (!all(receptor.cells %in% active_idents)) {
      stop('Not all provided \"receptor.cells\" are in active.idents.')
    }
  }

  # Split Pair.Name into ligand and receptor names. Retain those that are
  # present in the expression data.
  all_genes <- rownames(slot(object = seurat.object[[assay]], name = 'data'))
  lr_ref <- lr_ref[!duplicated(lr_ref[['Pair.Name']]),]
  tmp <- strsplit(x = lr_ref[['Pair.Name']], split = '_')
  ligand_genes <- sapply(X = tmp, FUN = `[`, 1)
  receptor_genes <- sapply(X = tmp, FUN = `[`, 2)
  ligand_present <- ligand_genes %in% all_genes
  receptor_present <- receptor_genes %in% all_genes
  lr_ref <- lr_ref[ligand_present & receptor_present, ]
  if (nrow(lr_ref) == 0) {
    stop('No LR pairs were detected in provided Seurat data.')
  } else {
    ligand_genes <- ligand_genes[ligand_present & receptor_present]
    receptor_genes <- receptor_genes[ligand_present & receptor_present]
    message(paste0('Genes detected for ', nrow(lr_ref), ' LR pairs'))
  }
  suppressWarnings(rm(tmp, ligand_present, receptor_present, all_genes))

  # Time-stamping
  t1 <- Sys.time()

  # Extract data. Convert factor-type column variables to character-type.
  Seurat::DefaultAssay(seurat.object) <- assay
  lr_data <- Seurat::FetchData(
    object = seurat.object,
    vars = unique(c(split.by, ligand_genes, receptor_genes)),
    slot = slot
  )
  is_factor <- which(sapply(X = lr_data, FUN = class) == 'factor')
  lr_data[is_factor] <- lapply(
    X = lr_data[is_factor],
    FUN = as.character
  )
  lr_data <- cbind(active_idents, lr_data)
  message(paste0('Using expression values for ', ncol(lr_data)-1, ' genes across ',
          nrow(lr_data), ' cells.'))
  suppressWarnings(rm(active_idents, is_factor))

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
    summarise(
      across(.cols = where(is.numeric),
             .fns = mean
      )
    )
  exp_pct <- exp_pct %>%
    summarise(
      across(.cols = where(is.numeric),
             .fns = function(x) round(mean(x > 0), 3)
      )
    )

  # Determine which genes meet minimum percent detection threshold
  minpct_genes <- sapply(
    FUN = function(x, min.pct) {any(x > min.pct)},
    X = exp_pct[sapply(exp_pct, is.numeric)],
    min.pct = min.pct
  )
  minpct_genes <- names(which(minpct_genes))
  ligand_detected <- ligand_genes %in% minpct_genes
  receptor_detected <- receptor_genes %in% minpct_genes
  ligand_genes <- ligand_genes[ligand_detected & receptor_detected]
  receptor_genes <- receptor_genes[ligand_detected & receptor_detected]
  lr_ref <- lr_ref[ligand_detected & receptor_detected,]
  suppressWarnings(rm(minpct_genes, ligand_detected, receptor_detected))

  # Generate all possible combinations of cell-type pairs (and split.by). Create
  # columns with cell-type counts for permutation test sampling.
  if (is.null(ligand.cells)) {
    var1 <- sort(unique(lr_data[['active_idents']]))
  } else {
    var1 <- ligand.cells
  }
  if (is.null(receptor.cells)) {
    var2 <- sort(unique(lr_data[['active_idents']]))
  } else {
    var2 <- receptor.cells
  }
  if (!is.null(split.by)) {
    if (is.null(split.vars)) {
      var3 <- sort(unique(lr_data[[split.by]]))
    } else {
      var3 <- split.vars
    }
    var_set <- expand.grid(var1, var2, var3, stringsAsFactors = FALSE)
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell', 'split.by')
    cell_counts <- table(
      slot(object = seurat.object, name = 'active.ident'),
      slot(object = seurat.object, name = 'meta.data')[[split.by]]
    )
  } else {
    var_set <- expand.grid(var1, var2, stringsAsFactors = FALSE)
    colnames(var_set) <- c('Ligand_cell', 'Receptor_cell')
    cell_counts <- table(slot(object = seurat.object, name = 'active.ident'))
  }
  var_set$Ligand_cell_count <- NA
  var_set$Receptor_cell_count <- NA
  for (i in 1:nrow(var_set)) {
    tmp_row_lig <- which(rownames(cell_counts) == var_set$Ligand_cell[i])
    tmp_row_rec <- which(rownames(cell_counts) == var_set$Receptor_cell[i])
    tmp_col <- which(colnames(cell_counts) == var_set$split.by[i])
    var_set$Ligand_cell_count[i] <- cell_counts[tmp_row_lig, tmp_col]
    var_set$Receptor_cell_count[i] <- cell_counts[tmp_row_rec, tmp_col]
  }
  suppressWarnings(rm(tmp_row_lig, tmp_row_rec, tmp_col, var1, var2, var3, i))

  # Calculate LR scores (for data, not nulls). Convert to list to vectorize.
  index_x <- var_set[['Ligand_cell']]
  index_y <- var_set[['Receptor_cell']]
  exp_names <- exp_avg[['active_idents']]
  if (!is.null(split.by)) {
    index_x <- paste(index_x, var_set[['split.by']], sep = '_')
    index_y <- paste(index_y, var_set[['split.by']], sep = '_')
    exp_names <- paste(exp_names, exp_avg[[split.by]], sep = '_')
  }
  index_x <- match(x = index_x, table = exp_names)
  index_y <- match(x = index_y, table = exp_names)
  index_l <- match(ligand_genes, colnames(exp_avg))
  index_r <- match(receptor_genes, colnames(exp_avg))
  exp_avg_lig <- as.matrix(exp_avg[index_x, index_l])
  rownames(exp_avg_lig) <- var_set[['Ligand_cell']]
  exp_avg_rec <- as.matrix(exp_avg[index_y, index_r])
  rownames(exp_avg_rec) <- var_set[['Receptor_cell']]
  exp_scores <- 1/2 * (exp_avg_lig + exp_avg_rec)
  rownames(exp_scores) <- apply(
    X = var_set,
    MARGIN = 1,
    FUN = function(x) {
      if (!is.null(split.by)) {
        paste(x[c('Ligand_cell','Receptor_cell','split.by')], collapse = '_')
      } else {
        paste(x[c('Ligand_cell','Receptor_cell')], collapse = '_')
      }
    }
  )
  colnames(exp_scores) <- paste(colnames(exp_avg_lig), colnames(exp_avg_rec),
                                sep = '_')

  # Extract individual matrices for ligand/receptor gene expression values to be
  # used for random sampling of cells (via matrix multiplication). Reformat the
  # lr_data data.frame to data.table before setting matrix to allow duplicated
  # column names (gene names).
  lr_data_ligands <- as.matrix(
    x = data.table::as.data.table(lr_data)[, ligand_genes, with = FALSE]
  )
  rownames(lr_data_ligands) <- rownames(lr_data)
  lr_data_receptors <- as.matrix(
    x = data.table::as.data.table(lr_data)[, receptor_genes, with = FALSE]
  )
  rownames(lr_data_receptors) <- rownames(lr_data)

  message('Generating null distributions...')
  null_scores <- bpmapply(
    FUN = GenerateNulls,
    l.count = var_set$Ligand_cell_count,
    r.count = var_set$Receptor_cell_count,
    MoreArgs = list(
      total.count = nrow(lr_data),
      resample = resample,
      lr.data.ligands = lr_data_ligands,
      lr.data.receptors = lr_data_receptors,
      Sample_Random_Cells = Sample_Random_Cells
    ),
    BPPARAM = BPPARAM
  ) # [resample, lr_pairs, cell_pairs]

  message('Calculating p-values...')
  pvals <- vector(mode = 'list', length = nrow(var_set))
  for (i in 1:length(pvals)) {
    pvals[[i]] <- mapply(
      exp.score = exp_scores[i,],
      null.scores = as.data.frame(null_scores[,,i]),
      FUN = CalculatePvals
    )
  }
  names(pvals) <- rownames(exp_scores)

  if (adjust.pval) {
    adj_pvals <- vector(mode = 'list', length = nrow(var_set))
    max_null <- null_scores[,,1]
    for (i in 1:length(adj_pvals)) {
      max_null[null_scores[,,i] > max_null] <- null_scores[,,i][null_scores[,,i] > max_null]
    }
    for (i in 1:nrow(var_set)) {
      adj_pvals[[i]] <- mapply(
        exp.score = exp_scores[i,],
        null.scores = as.data.frame(max_null),
        FUN = CalculatePvals
      )
    }
    names(adj_pvals) <- rownames(exp_scores)
  }

  # Long-form data.table of results
  exp_scores <- data.frame(exp_scores, stringsAsFactors = FALSE) %>%
    reshape2::melt()
  pvals <- t(data.frame(pvals, stringsAsFactors = FALSE)) %>%
    reshape2::melt()
  exp_avg_lig <- exp_avg_lig %>% reshape2::melt()
  exp_avg_rec <- exp_avg_rec %>% reshape2::melt()
  exp_pct_lig <- exp_pct[index_x, index_l] %>% reshape2::melt()
  exp_pct_rec <- exp_pct[index_y, index_r] %>% reshape2::melt()

  tmp_idents <- strsplit(
    x = gsub(
      pattern = '\\.',
      replacement = '-',
      x = pvals[['Var1']]
    ),
    split = '_'
  )
  tmp_pairs <- strsplit(
    x = gsub(
      pattern = '\\.',
      replacement = '-',
      x = pvals[['Var2']]
    ),
    split = '_'
  )
  ligand_cell <- sapply(X = tmp_idents, FUN = `[[`, 1)
  receptor_cell <- sapply(X = tmp_idents, FUN = `[[`, 2)
  ligand <- sapply(X = tmp_pairs, FUN = `[[`, 1)
  receptor <- sapply(X = tmp_pairs, FUN = `[[`, 2)
  cell_pair <- paste(ligand_cell, receptor_cell, sep = '_')
  lr_pair <- paste(ligand, receptor, sep = '_')
  tmp_scores <- exp_scores[['value']]
  tmp_pvals <- pvals[['value']]
  tmp_avg_l <- exp_avg_lig[['value']]
  tmp_avg_r <- exp_avg_rec[['value']]
  tmp_pct_l <- exp_pct_lig[['value']]
  tmp_pct_r <- exp_pct_rec[['value']]
  if (adjust.pval) {
    adj_pvals <- t(data.frame(adj_pvals, stringsAsFactors = FALSE)) %>%
      reshape2::melt()
    tmp_adj_pvals <- adj_pvals[['value']]
  } else {
    tmp_adj_pvals <- NA
  }

  if (!is.null(split.by)) {
    split_var <- sapply(X = tmp_idents, FUN = `[[`, 3)
    cell_prop <- round(prop.table(cell_counts, margin = 1), 3)*100
    tmp_count_l <- mapply(
      FUN = Get_Cell_Proportion,
      cell.type = ligand_cell,
      split.var = split_var,
      MoreArgs = list(cell.table = cell_prop)
    )
    tmp_count_r <- mapply(
      FUN = Get_Cell_Proportion,
      cell.type = receptor_cell,
      split.var = split_var,
      MoreArgs = list(cell.table = cell_prop)
    )
  }

  # Final result compilation
  results <- data.frame(
    'Pair_name' = lr_pair,
    'Score' = tmp_scores,
    'pval' = tmp_pvals,
    'adj_pval' = tmp_adj_pvals,
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

  # Finish
  t1 <- Sys.time() - t1
  message('Done!')
  message('Run time:')
  print(t1)
  return(results)
}


GenerateNulls <- function(
  resample,
  total.count,
  l.count,
  r.count,
  lr.data.ligands,
  lr.data.receptors,
  Sample_Random_Cells
) {
  null_index_x <- t( # nrow = resample, ncol = nrow(lr_data)
    x = replicate(
      n = resample,
      expr = Sample_Random_Cells(
        sample.count = l.count,
        total.count = total.count
      )
    )
  )
  null_index_y <- t(
    x = replicate(
      n = resample,
      expr = Sample_Random_Cells(
        sample.count = r.count,
        total.count = total.count
      )
    )
  )
  null_avg_lig <- (null_index_x %*% lr.data.ligands) / l.count
  null_avg_rec <- (null_index_y %*% lr.data.receptors) / r.count
  null_scores <- 0.5 * (null_avg_lig + null_avg_rec) # nrow = resample, ncol = # lr pairs

  return(null_scores)
}


CalculatePvals <- function(
  null.scores,
  exp.score
) {
  p <- 1 - (sum(exp.score >= null.scores) / length(null.scores))
  return(p)
}


Sample_Random_Cells <- function(
  sample.count,
  total.count
) {
  tmp <- rep(FALSE, times = total.count)
  tmp_switch <- sample(x = total.count, size = sample.count, replace = FALSE)
  tmp[tmp_switch] <- TRUE
  return(tmp)
}


Get_Cell_Proportion <- function(
  cell.table,
  cell.type,
  split.var
) {
  prop <- cell.table[cell.type, split.var]
  return(prop)
}
