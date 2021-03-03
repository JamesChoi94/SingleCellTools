#' Ligand-Receptor plotting - mean of mean expressions
#'
#'
#'
#'
#' @import dplyr
#' @import tibble
#' @importFrom rlang "!!"
#'
#' @param results data.frame containing output from \code{\link{StandardLR}}.
#' @param ligands Character vector with specific ligands of interest to plot. If
#' \code{NULL}, all available ligand genes will be plotted.
#' @param receptors Character vector with sepcific receptors of interest to
#' plot. If \code{NULL}, all available receptor genes will be plotted.
#' @param l.cells Character vector with specific cell clusters to plot as the
#' ligand-expressing cell-type. If \code{NULL}, all available cell cluster types
#' will be plotted.
#' @param r.cells Character vector with specific cell clusters to plot as the
#' receptor-expressing cell-type. If \code{NULL}, all available cell cluster
#' types will be plotted.
#' @param split.subset Character vector with values of \code{split.by} column
#' from \code{results} to retain. All other values will be removed from
#' plotting. If \code{NULL}, all available values of \code{split.by} will be
#' plotted.
#' @param split.along.y Logical whether \code{split.by} plot facets should run
#' horizontally (\code{FALSE}) or vertically (\code{TRUE}).
#' @param use.adj.pval Logical whether to plot adjusted p-values column from
#' \code{results}. Default is \code{FALSE}.
#' @param pval.threshold Numeric value of the maximum p-value or adjusted
#' p-value to plot. All scores with p-values greater than the threshold will be
#' removed.
#' @param min.exp.percent Numeric value of the minimum percent that a cell
#' cluster must express a gene in order to be plotted. Percent expresses below
#' this value will not be plotted.
#' @param min.cell.percent Numeric value (need to expand more).
#' @param resample Numeric value of the number of resamplings used in the
#' permutation test.
#'
#' @return A ggplot2 object with the ligand-receptor plot.
#'
#' @export
#'

PlotLR <- function(
  results,
  ligands = NULL,
  receptors = NULL,
  l.cells = NULL,
  r.cells = NULL,
  split.subset = NULL,
  split.along.y = FALSE,
  use.adj.pval = FALSE,
  pval.threshold = 0.05,
  min.exp.percent = 0.1,
  min.cell.percent = 0,
  resample = 1000
) {
  # Import locally
  tmp_results <- results

  # Remove p-values for which expression percents do not meet provided threshold
  tmp_results[['pval']] <- ifelse(test = results[['Ligand_pct']] < min.exp.percent |
                                    results[['Receptor_pct']] < min.exp.percent,
                                  yes = NA,
                                  no = results[['pval']])
  if(use.adj.pval) {
    if(all(is.na(tmp_results[['adj_pval']]))) {
      stop('Adjusted p-values not calculated')
    }
    tmp_results[['adj_pval']] <- ifelse(test = results[['Ligand_pct']] < min.exp.percent |
                                          results[['Receptor_pct']] < min.exp.percent,
                                        yes = NA,
                                        no = results[['adj_pval']])
  }


  # Series of filters based on those provided by user.
  if(!is.null(ligands)) {
    tmp_results <- tmp_results %>% filter(Ligand %in% ligands)
    if(nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if(!is.null(receptors)) {
    tmp_results <- tmp_results %>% filter(Receptor %in% receptors)
    if(nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if(!is.null(l.cells)) {
    tmp_results <- tmp_results %>% filter(Ligand_cell %in% l.cells)
    if(nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if(!is.null(r.cells)) {
    tmp_results <- tmp_results %>% filter(Receptor_cell %in% r.cells)
    if(nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  if(!is.null(pval.threshold)) {
    tmp_results <- tmp_results %>% filter(pval < pval.threshold)
    if(nrow(tmp_results) == 0) {
      stop("No remaining results after filter.")
    }
  }
  # Remove rows where at least "min.cell.percent" of a cell-type is not present in a given subset of "split_by".
  if(!is.null(split.subset)) {
    if(!all(is.na(tmp_results[['split_by']]))) {
      stop('Provided \"split.subset\" but no values present in \"split_by\" column of \"results\".')
    }
    tmp_results <- tmp_results %>% filter(split_by %in% split.subset)
    if(nrow(tmp_results) == 0) {
      stop('No remaining results after filter.')
    }
  }
  tmp_results <- tmp_results %>% ungroup()
  tmp_pairs <- unique(tmp_results[['Pair_name']][!is.na(tmp_results[['pval']])])
  tmp_results <- tmp_results[tmp_results[['Pair_name']] %in% tmp_pairs,]

  # Calculate -log10 of p-values for visualization.
  min_col <- floor(min(tmp_results[['Score']])/0.5) * 0.5
  tmp_results[['log_pval']] <- -log10(tmp_results[['pval']])
  tmp_results[['log_pval']][is.infinite(tmp_results[['log_pval']])] <- log10(resample)
  if(use.adj.pval) {
    tmp_results[['log_adj_pval']] <- -log10(tmp_results[['adj_pval']])
    tmp_results[['log_adj_pval']][is.infinite(tmp_results[['log_adj_pval']])] <- log10(resample)
  }
  if(use.adj.pval) {
    tmp_pval <- 'log_adj_pval'
  } else {
    tmp_pval <- 'log_pval'
  }
  tmp_plot <- tmp_results %>%
    ggplot() +
    geom_point(mapping = aes_string(x = 'Pair_name', y = 'Receptor_cell', size = tmp_pval, fill = 'Score'),
               color = 'black', pch = 21)
  if(!all(is.na(tmp_results[['split_by']]))) {
    if(split.along.y) {
      tmp_plot <- tmp_plot + facet_grid(Ligand_cell + split_by ~ ., switch = 'y', drop = TRUE)
    } else {
      tmp_plot <- tmp_plot + facet_grid(Ligand_cell ~ split_by, switch = 'y', drop = TRUE)
    }
  } else {
    tmp_plot <- tmp_plot + facet_grid(Ligand_cell ~ ., switch = 'y')
  }
  tmp_plot <- tmp_plot +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = 'Spectral')),
                         limits = c(min_col, NA)) +
    scale_radius(limits = c(0,NA), range = c(1,6)) +
    scale_y_discrete(position = 'right') +
    theme(strip.text = element_text(size = 12, color = 'black', face = 'bold'),
          strip.background = element_rect(fill = NA, color = 'black'),
          axis.title = element_blank(),
          axis.text = element_text(size = 12, color = 'black'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.text = element_text(size = 12, color = 'black'),
          legend.title = element_text(size = 12, color = 'black', face = 'bold'),
          legend.key = element_rect(fill = NA),
          panel.background = element_rect(fill = NA, color = 'black'),
          panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'grey70')) +
    guides(fill = guide_colorbar(title = 'Score',
                                 frame.linewidth = 1,
                                 ticks.linewidth = 1,
                                 frame.colour = 'black',
                                 ticks.colour = 'black'),
           size = guide_legend(title = '-log10(p-value)',
                               override.aes = list(fill = 'black')))
  return(tmp_plot)
}
