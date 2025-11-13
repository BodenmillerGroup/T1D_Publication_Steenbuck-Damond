library(ggplot2)
library(dplyr)
library(SpatialExperiment)
library(SummarizedExperiment)

# Helper functions for T1D analysis
#' List of functions
#' - `getPaths`: get default paths for the analysis (change the values if needed).
#' - `censor_dat`: removes the outliers on the upper side by capping the values at the provided quantile (copied from bbRtools).
#' - `summarizeHeatmap`: function to return mean or median counts by cluster and by channel from a `SingleCellExperiment` object.
#' - `plotDimRed`: plot reduced dimension and color by variable.
#' - `plotDimRedchannels`: plot markers on reduced dimensions.
#' - `plotBoxes`: plot boxplots.
#' - `plotViolins`: plot violin plots.
#' - `plotDensity`: plot density plots.
#' - `calcCompo` : calculate cell type composition.
#' - `YOLOader` : there are two ways to build `CytoImageList` objects: the proper way and YOLOader.
#' - `sample_one`: sample exactly one element (copied from bbsnippets, heatscatter).
#' - `mytheme`: Plotting theme for ggplot2

#' getPaths
#' @return paths of data, git and output folders
getPaths <- function(script_name) {
  dataset_name <- "T1D"
  home <- file.path("/", "mnt", "central_nas", "projects",
                     "type1_diabetes", "nathan", "T1D_Vol")

  home_data <- file.path(home)

  home_analysis <- file.path(home, "T1D_analysis")

  cluster_home <- file.path("/", "data", "nsteen")

  script_name <- gsub(".Rmd", "", script_name)
  panel_type <- gsub("[^_]*_[^_]*_[^_]*_", "", script_name)
  object_type <- gsub(paste0("_", panel_type), "", script_name)
  object_type <- gsub("[^_]*_[^_]*_", "", object_type)

  exp_name <- paste(dataset_name, panel_type, object_type, sep = "_")

  folder_in <- file.path(home_data, "processing")
  cluster_folder_in <- file.path("/", "scratch", "nsteen", "processing")
  
  folder_out <- file.path(home_data, "temp_analysis", "analysis")
  if (!dir.exists(folder_out)) dir.create(folder_out)

  folder_script <- file.path(folder_out, script_name)
  if (!dir.exists(folder_script)) dir.create(folder_script)
  cluster_folder_script <- file.path(cluster_home, script_name)

  paths <- list(
    dataset_name = dataset_name,
    panel_type = panel_type,
    object_type = object_type,
    exp_name = exp_name,
    home = home,
    home_git = home_analysis,
    home_data = home_data,
    folder_in = folder_in,
    folder_out = folder_out,
    folder_script = folder_script,
    cluster_folder_script = cluster_folder_script,
    cluster_home = cluster_home,
    cluster_folder_in = cluster_folder_in
  )

  return(paths)
}


#' censor_dat
#' Function copied from [bbRtools](https://github.com/BodenmillerGroup/bbRtools).
#' Removes the outliers on the upper side by capping the values at the provided quantile.
#'
#' @param x values to censor
#' @param quant quantile to censor, i.e. how many percent of values are considered outliers
#' @param symmetric censor on both side. In this case the outliers are assumed to be symetric on both sides. For example if a quantile of 5\% (0.05) is choosen, in the symetric case 2.5\% (0.025) of values are censored on both sides.
#'
#' @return returns the percentile of each value of x
#' @export

censor_dat <- function(x, quant = 0.999, symmetric = FALSE) {
  if (symmetric) {
    lower_quant <- (1 - quant)/2
    quant <- quant + lower_quant
  }
  q <- stats::quantile(x, quant)
  x[x > q] <- q
  if (symmetric) {
    q <- stats::quantile(x, lower_quant)
    x[x < q] <- q
  }
  return(x)
}


#' #' summarize_heatmap
#' #' Returns median counts by cluster and by channel.
#' #'
#' #' @param x A `SingleCellExperiment` or `SpatialExperiment` object.
#' #' @param expr_values A string corresponding to an assay in x, should be in `assayNames(x)`.
#' #' @param cluster_by Name of the column containing the clusters.
#' #' @param channels Channels to include, should be in `rownames(x)`. If `NULL`, all channels will be summarized.
#' #' @param fun c("median","mean") chose if the median or the mean should be returned (default:"median") (optional).
#' #'
#' #' @return summarized data as a matrix.
#' #' @export
#'
summarize_heatmap <- function(x, expr_values, cluster_by, channels = NULL, fun = "median") {
  require(data.table)
  # Argument checks
  if (is.null(expr_values) || !(expr_values %in% SummarizedExperiment::assayNames(x))) {
    expr_values <- SummarizedExperiment::assayNames(x)[1]
    print(paste0("Warning: Assay type not provided or assay type not in 'assayNames(x)', '", expr_values, "' used."))
  }
  if (is.null(channels)) {
    channels <- rownames(x)
  }
  if (!all(channels %in% rownames(x))) {
    stop("Channel names do not correspond to the rownames of the assay")
  }
  if (is.null(cluster_by)) {
    stop("Cluster column not provided")
  }
  if (!(cluster_by %in% colnames(SummarizedExperiment::colData(x)))) {
    stop("The 'cluster_by' argument should correspond to a colData(x) column")
  }
  if (length(cluster_by) > 1) {
    stop("'cluster_by' takes only one argument")
  }

  if (length(expr_values) > 1) {
    stop("'expr_values' takes only one argument")
  }
  if (!(fun %in% c("median", "mean"))) {
    stop("'fun' takes either 'median' or 'mean' as an argument")
  }

  # Convert the data to a melted format
  if (expr_values %in% SummarizedExperiment::assayNames(x)) {
    dat <- as.data.table(t(SummarizedExperiment::assay(x, expr_values)[channels, ]))
  } else if (expr_values %in% SummarizedExperiment::reducedDimNames(x)) {
    dat <- as.data.table(t(SummarizedExperiment::assay(x, "exprs")[channels, ]))
  }
  dat[, id := colnames(x)]
  dat[, cluster := SummarizedExperiment::colData(x)[, cluster_by]]
  dat <- melt.data.table(dat,
                         id.vars = c("id", "cluster"),
                         variable.name = "channel",
                         value.name = expr_values)

  # Summarize the data
  # Get mean/median for each cluster and channel combination
  if (fun == "median") {
    dat_summary <- dat[, list(
      summarized_val = median(get(expr_values)),
      cellspercluster = .N),
      by = c("channel", "cluster")]
  } else if (fun == "mean") {
    dat_summary <- dat[, list(
      summarized_val = mean(get(expr_values)),
      cellspercluster = .N),
      by = c("channel", "cluster")]
  }
  # Decast the summarized data and convert to a matrix
  hm_cell <- dcast.data.table(dat_summary,
                              formula = "cluster ~ channel",
                              value.var = "summarized_val")
  hm_clusters <- hm_cell$cluster
  hm_cell <- as.matrix(hm_cell[, -1, with = FALSE])

  # Add rownames
  rownames(hm_cell) <- hm_clusters
  # Return the summarized values
  return(as.matrix(hm_cell))
}

#' sample_one: samples exactly one value (copied from Vito's heatscatter snippet in bbsnippets)
sample_one <- function(x) {
  # Helper function to sample exactly 1 element
  if (length(x) == 0) {
    return(NA)
  } else if (length(x) == 1) {
    return(x)
  } else {
    return(sample(x, 1))
  }
}


#' plot_dim_red
#' Plot a variable on reduced dimensions
#'
#' @param dat A `data.frame`.
#' @param dimred the name of the reduced dimensions to use.
#' @param color_by Name of the `colData(sce)`column containing the variable to color by.
#' @param sample (boolean) should the rows of `dat` be randomly sampled.
#' @param size point size.
#' @param alpha point alpha.
#' @param axes.label vector indicating the labels of the x and y axes.
#' @param palette vector indicating the color values to use.
#' @param palette_continuous (boolean) should a continuous (rather than discrete) color palette by used.
#'
#' @return a `ggplot` object
#' @export

plot_dim_red <- function(dat, dimred, color_by, sample = TRUE,
                         size = 0.1, alpha = 0.5, axes_labels = c("Reduced dim 1", "Reduced dim 2"),
                         palette = NULL, palette_continuous = FALSE) {
  if (sample == TRUE)
    dat <- dat[sample(nrow(dat)), ]

  p <- dat |>
    ggplot2::ggplot(ggplot2::aes(x = get(paste0(dimred, ".1")),
                                 y = get(paste0(dimred, ".2")))) +
    ggplot2::ggtitle(paste(dimred, color_by, sep = "-")) +
    ggplot2::guides(color = ggplot2::guide_legend(title = color_by, override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(x = axes_labels[1], y = axes_labels[2]) +
    mytheme$standard() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
    )

  nval <- nrow(unique(dat[, color_by]))
  if (is.null(nval)) { 
    nval <- length(unique(dat[, color_by]))
  }
  # nval.pal <- ifelse(is.null(palette), nval, length(palette))


  if (isFALSE(palette_continuous)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = as.factor(get(color_by))),
                          size = size, alpha = alpha)
    if (is.null(palette) && nval > 200) {
      stop("Please provide a palette with enough colors or set `palette_continuous` as `TRUE`")
    } else if (!is.null(palette)) {
      p <- p + ggplot2::scale_colour_manual(values = palette)
    } else if (nval <= 15) {
      p <- p + ggplot2::scale_colour_manual(values = palettes$colors[c(FALSE, TRUE)])
    } else if (nval <= 30) {
      p <- p + ggplot2::scale_colour_manual(values = palettes$colors)
    } else if (nval <= 50) {
      p <- p + ggplot2::scale_colour_manual(values = palettes$colors50)
    } else {
      p <- p + ggplot2::scale_colour_manual(values = palettes$colorsfull)
    }  
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = get(color_by)),
                          size = size, alpha = alpha)
    if (!is.null(palette)) {
      p <- p + ggplot2::scale_color_continuous(type = palette)
    } else {
      require(viridis)
      p <- p + scale_color_viridis(option = "viridis")
    }
  }
  return(p)
}


#' plot_dim_red_channels
#' Plot markers on reduced dimensions
#'
#' @param dat A `data.frame`.
#' @param dimred the name of the reduced dimensions to use.
#' @param expr_values the name of the assay to use.
#' @param channels a vector containing the channels to use.
#' @param censor_val counts censor value (percentile clipping).
#' @param ncol number of column in facetted plot.
#' @param size point size (only used if data is displayed as points).
#' @param alpha point alpha (only used if data is displayed as points).
#' @param axes.label vector indicating the labels of the x and y axes.
#' @param force_points (boolean) force display of points rather than summary statistics.
#'
#' @return a `ggplot` object
#' @export

plot_dim_red_channels <- function(
    dat, dimred, expr_values, channels, censor_val = 0.95,
    axes_labels = c("Reduced dim 1", "Reduced dim 2"),
    ncol = ceiling(sqrt(length(channels))), size = 0.5, alpha = 1,
    force_points = FALSE) {

  p <- dat[dat$channel %in% channels & sample(nrow(dat)), ] |>
    ggplot2::ggplot(ggplot2::aes(x = get(paste0(dimred, ".1")),
                                 y = get(paste0(dimred, ".2")))) +
    ggplot2::facet_wrap(~channel, ncol = ncol) +
    ggplot2::ggtitle(paste(dimred, expr_values, sep = "-")) +
    ggplot2::labs(x = axes_labels[1], y = axes_labels[2]) +
    mytheme$large() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
    )

  if (isTRUE(force_points)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = censor_dat(get(expr_values), censor_val)),
                                 size = size, alpha = alpha) +
      viridis::scale_color_viridis(name = expr_values, option = "viridis")
  } else {
    p <- p + ggplot2::stat_summary_2d(ggplot2::aes(z = censor_dat(get(expr_values), censor_val)),
                                      bins = 500, fun = sample_one) +
      viridis::scale_fill_viridis(name = expr_values, option = "viridis")
  }

  return(p)
}

#' plot_boxes
#' Plot boxplots
#'
#' @param dat A `data.frame`.
#' @param x x axis variable.
#' @param y y axis variable.
#' @param fill_by Variable to fill by.
#' @param color_by (optional) Variable to color by.
#' @param facet_by (optional) Variable to facet by.
#' @param ncol number of column in facetted plot.
#' @param title Graph title.
#'
#' @return a `ggplot` object
#' @export

plot_boxes <- function(dat, x, y, fill_by, color_by = NULL,
                       facet_by = NULL, ncol = NULL, title = NULL, scales = "fixed") {

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = get(x), y = get(y))) +
    viridis::scale_fill_viridis() +
    mytheme$standard() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.background = ggplot2::element_rect(colour = "grey20")
    )

  if (is.null(color_by))
    p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = get(fill_by)),
                                   outlier.size = 0.1, lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by)
  else
    p <- p + geom_boxplot(ggplot2::aes(fill = get(fill_by), color = get(color_by)),
                          outlier.size = 0.1, lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by, color = color_by)

  if (!is.null(facet_by))
    p <- p + ggplot2::facet_wrap(~get(facet_by), scale = scales, ncol = ncol)

  if (is.null(title))
    p <- p + ggplot2::ggtitle(paste(x, y, sep = " - "))
  else
    p <- p + ggtitle(title)

  return(p)
}

#' plot_violins
#' Plot violin plots
#'
#' @param dat A `data.frame`.
#' @param x x axis variable.
#' @param y y axis variable.
#' @param fill_by Variable to fill by.
#' @param color_by (optional) Variable to color by.
#' @param facet_by (optional) Variable to facet by.
#' @param ncol number of column in facetted plot.
#' @param title Graph title.
#'
#' @return a `ggplot` object
#' @export

plot_violin <- function(dat, x, y, fill_by, color_by = NULL,
                        facet_by = NULL, ncol = NULL, title = NULL,
                        scales = "fixed") {

  p <- ggplot(dat, aes(x = get(x), y = get(y))) +
    viridis::scale_fill_viridis() +
    mytheme$standard() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(colour = "grey20")
    )

  if (is.null(color_by))
    p <- p + geom_violin(ggplot2::aes(fill = get(fill_by)),
                         draw_quantiles = 0.5, scale = "width", lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by)
  else
    p <- p + geom_violin(ggplot2::aes(fill = get(fill_by), color = get(color_by)),
                         draw_quantiles = 0.5, scale = "width", lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by, color = color_by)

  if (!is.null(facet_by))
    p <- p + facet_wrap(~get(facet_by), ncol = ncol, scales = scales)

  if (is.null(title))
    p <- p + ggtitle(paste(x, y, sep = " - "))
  else
    p <- p + ggtitle(title)

  return(p)
}

## Expression of markers in violin plots, with significant brackets.
plot_expression_violin <- function(exprs_dat, channels_oi, cell_types, facet_by, y_label, immune_adjustment = FALSE, tip.length = 0.05) {
  exprs_dat2 <- exprs_dat |> 
    filter(cell_type %in% cell_types) |> 
    tibble::as_tibble() |>
    tidyr::drop_na(cell_type) |> 
    dplyr::filter(donor_type != "LongDuration") |> 
    dplyr::filter(channel %in% channels_oi) |> 
    dplyr::mutate(channel = factor(as.character(channel), levels = channels_oi)) |>
    dplyr::mutate(donor_type = ordered(donor_type, levels = c("NoDiabetes", "sAAb+", "mAAb+", "RecentOnset")))
  
  exprs_dat_plot <- exprs_dat2 |> 
    dplyr::group_by(case_id, donor_type, channel, cell_type) |>
    dplyr::summarise(exprs = mean(exprs), .groups = "drop")
  
  ## Calcualte LMM between ND and other stages.
  p_vals <- exprs_dat2 |> 
    tidyr::nest(data = -c(channel, cell_type)) |>
    dplyr::mutate(test = purrr::map(data, \(.x){
      # LMM
      test <- lmerTest::lmer(data = .x, exprs ~ donor_type + (1 | case_id))
      # Linear test:
      m.emm <- emmeans::emmeans(test, ~ donor_type, pbkrtest.limit = 10000)
      # Pairwise comparisons (except for sAAb+ vs RecentOnset -> not required)
      # FIXME: might need to include more if LD in there.
      pairwise <- m.emm |> 
        emmeans::contrast("pairwise", adjust = "none") |> 
        as_tibble() |> 
        filter(contrast != "(sAAb+) - RecentOnset")
      #p_vals_t <- contrast(m.emm, 'poly') |> 
      #  as_tibble() |> 
      #  filter(contrast == "linear")
      #bind_rows(pairwise, p_vals_t)
    })) |> 
    tidyr::unnest(test) |> 
    mutate(p.adj = p.adjust(p.value, method = "fdr"))

  exprs_dat_plot <- exprs_dat2 |> 
    group_by(case_id, donor_type, channel, cell_type) |>
    summarise(exprs = mean(exprs), .groups = "drop")

  ## This is for immune adjustment (age + ICU time)
  if (immune_adjustment) {
    p_vals <- exprs_dat2 |> 
      nest(data = -c(channel, cell_type)) |>
      dplyr::mutate(test = purrr::map(data, \(.x){
        # LMM
        test <- lmerTest::lmer(data = .x, exprs ~ donor_type + age_group.y + ICU_time_days.y + (1 | case_id))
        # Linear test:
        m.emm <- emmeans::emmeans(test, ~ donor_type, pbkrtest.limit = 30000)
        # Pairwise comparisons (except for sAAb+ vs RecentOnset -> not required)
        pairwise <- m.emm |> 
          emmeans::contrast("pairwise", adjust = "none") |> 
          as_tibble() |> 
          filter(contrast != "(sAAb+) - RecentOnset")
        #p_vals_t <- contrast(m.emm, 'poly') |> 
        #  as_tibble() |> 
        #  filter(contrast == "linear")
        #bind_rows(pairwise, p_vals_t)
      })) |> 
      tidyr::unnest(test) |> 
      mutate(p.adj = p.adjust(p.value, method = "fdr"))
  }

  p_vals <- p_vals |> 
    mutate(group1 = contrast, group2 = contrast) |>  
    mutate(group1 = stringr::str_replace(group1, "\\s.*", "")) |> 
    mutate(group2 = stringr::str_replace(group2, "^(?:\\S+\\s+\\S+\\s+)", "")) |> 
    mutate(group1 = stringr::str_remove(group1, "\\(")) |>
    mutate(group2 = stringr::str_remove(group2, "\\(")) |> 
    mutate(group1 = stringr::str_remove(group1, "\\)")) |> 
    mutate(group2 = stringr::str_remove(group2, "\\)"))

  # Use this for actual P-Values.
  p_vals <- p_vals |> 
    mutate(.y. = "exprs",
           p = p.value,
           statistic = t.ratio,
           p.adj = p.adjust(p, method = "fdr"),
           p.signif = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "•", ""))
           ) |> 
    select(channel, .y., group1, group2,
           statistic, p, p.adj, p.signif, cell_type)
  
  # get max-y-value PER channel.
  p_vals <- exprs_dat_plot |> 
      group_by(channel, cell_type) |> 
      summarise(max = max(exprs), .groups = "drop") |>
      mutate(pos_pvals = max*1.05) |> 
      dplyr::right_join(p_vals, by = c("channel", "cell_type"))

  ## Linear
  # p_vals_linear <- p_vals  |> filter(group1 == "linear")
  #p_vals_other <- p_vals |> 
  #  filter((group1 != "linear" & p.signif != "" )| group1 == "linear")
  
  ## Calculate p-values with ggpubr.
  # Use this for right formatting of DF
  stat.test <- exprs_dat_plot |> 
    group_by(channel, cell_type) |> 
    tukey_hsd(exprs ~ donor_type)  |> 
    filter(!(group1 == "sAAb+" & group2 == "RecentOnset")) |> 
    select(channel, group1, group2, cell_type) |> 
    add_xy_position(x = "donor_type", step.increase = 0.07, dodge = 0.8) |> 
    distinct()
  
  stat.test.test <- p_vals |> 
   left_join(stat.test, by = c("channel", "cell_type", "group1", "group2")) |> 
   # mutate(xmin = ifelse(group1 == "linear", 1, xmin),
   #       xmax = ifelse(group2 == "linear", 4, xmax),
   #       y.position = ifelse(group1 == "linear", max, y.position),
   #       groups = ifelse(group1 == "linear", 
   #             list(c("NoDiabetes", "RecentOnset")), groups)) |>
    mutate(y.position = max * 1.3)  |> 
    filter(p.adj < 0.1) |> 
    group_by(channel, cell_type) %>%
    mutate(y.position = dplyr::first(y.position) * (1 + (0.07) * (row_number() - 1))) %>%
    ungroup()
  #stat.test.test.2 <- stat.test.test |> 
  #  filter(group1 == "linear") |> 
  #  mutate(p.adj = ifelse(p.adj >= 1e-2, 
  #       formatC(p.adj, format = "f", digits = 2), # normal decimal
  #       formatC(p.adj, format = "e", digits = 1)  # scientific
  #))
  # 

  # Combine formatting + P-value dataframe.
  #stat.test.test <- p_vals |> 
  #  left_join(stat.test, by = c("channel", "cell_type", ".y.", "group1", "group2")) |> 
  #  mutate(y.position = max * 1.3) |> 
  #  mutate(y.position = ifelse(group2 == "mAAb+", y.position*1.05, y.position)) |>
  #  mutate(y.position = ifelse(group2 == "RecentOnset", y.position*1.1, y.position))

  # Plot.
  p <- ggviolin(exprs_dat_plot, x = "donor_type", y = "exprs", fill = "donor_type", 
          facet.by = facet_by, palette = palettes$stages, 
          draw_quantiles = c(0.5), add = "jitter", 
          add.params = list(width = 0.1, alpha = 1, size = 3)) + 
    stat_pvalue_manual(stat.test.test, label = "p.signif", tip.length = tip.length,
                      size = 10, label.size = 10, vjust = 0.4) + 
    #stat_pvalue_manual(stat.test.test.2, label = "p.adj", tip.length = 0,
    #                  size = 10, label.size = 10) + 
    # stat_pvalue_manual(p_vals_linear, label = "p.signif", tip.length = 0.01)
    facet_wrap(~get(facet_by), scales = "free_y") +
    labs(x = "Stage", y = y_label) + 
    mytheme$large_new() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+ 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  return(p)
}

## Expression of markers between Insulitic vs Non-Insulitic OR between ICI and IDI.
plot_ici_insulitic_expression <- function(df, features_oi, cell_types, cur_var, x_var, group1, group2, facet_by, y_label, cur_palette) {

  exprs_macs <- df |> 
    filter(cell_type %in% cell_types) |>
    tibble::as_tibble() |>
    tidyr::drop_na(cell_type) |>
    filter(donor_type %in% c("RecentOnset")) |>
    pivot_longer(cols = all_of(features_oi), names_to = "channel", values_to = "MeanExprs") |> 
    mutate(channel = factor(channel, levels = features_oi))

  if (cur_var == "insulitic") {
    exprs_macs_islet <- exprs_macs |> 
      group_by(case_id, donor_type, image_id, insulitic, cell_type, channel) |>
      summarise(MeanExprs = mean(MeanExprs), .groups = "drop")
  } else {
    exprs_macs_islet <- exprs_macs |> 
      group_by(case_id, donor_type, image_id, islet_type, cell_type, channel) |>
      summarise(MeanExprs = mean(MeanExprs), .groups = "drop")
  }
  
  cols_islet_plot <- exprs_macs |> 
    dplyr::rename(cur_var = !!cur_var) |>
    group_by(case_id, donor_type, cell_type, channel, cur_var) |>
    summarise(MeanExprs = mean(MeanExprs), .groups = "drop")

  ## Calculate LMM between ND and other stages.
  cur_formula <- formula(paste("MeanExprs ~", cur_var, "+ (1 | case_id)"))
  # cur_formula_2 <- formula(paste("MeanExprs ~", cur_var, "+ (", cur_var, "| case_id)"))

  # lmer_intercept <- lmerTest::lmer(data = exprs_macs |> filter(channel == "XBP1"), cur_formula)
  # lmer_slope <- lmerTest::lmer(data = exprs_macs |> filter(channel == "XBP1"), cur_formula_2)
  # anova(lmer_intercept, lmer_slope)

  p_vals <- exprs_macs_islet  |> 
    nest(data = -c(channel, cell_type)) |> 
    dplyr::mutate(test = purrr::map(data, ~ broom.mixed::tidy(lmerTest::lmer(data = ., cur_formula)))) |>
    tidyr::unnest(test) |>
    ungroup() |> 
    filter(effect == "fixed", term != "(Intercept)") |> 
    mutate(islet_type = stringr::str_remove(term, "islet_type")) |> 
    mutate(insulitic = stringr::str_remove(term, "insulitic")) |>
    select(-c(data, term))

  # p_vals <- bind_rows(p_vals1, p_vals2)

  # Use this for actual P-Values.
  p_vals <- p_vals |> 
    mutate(.y. = "MeanExprs",
          group1 = group1, group2 = group2,
          p = p.value,
          p.adj = p.adjust(p, method = "fdr"),
           p.signif = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "•", ""))
          ) |> 
    select(channel, .y., group1, group2,
          statistic, p, p.adj, p.signif, cell_type, estimate)

  # get max-y-value PER channel.
  p_vals <- cols_islet_plot  |> 
    group_by(channel, cell_type) |> 
    summarise(max = max(MeanExprs), .groups = "drop") |>
    mutate(pos_pvals = max*1.05) |> 
    dplyr::right_join(p_vals, by = c("channel", "cell_type"))

  ## Calculate p-values with ggpubr.
  #  Use this for right formatting of DF
  cur_formula2 <- formula(paste("MeanExprs ~", cur_var))
  stat.test <- exprs_macs |> 
    group_by(channel, cell_type) |> 
    t_test(cur_formula2, ref.group = group1) |> 
    adjust_pvalue()  |> 
    add_significance("p.adj") |> 
    add_xy_position(x = x_var, step.increase = 0.07, dodge = 0.8) |> 
    select(-c(statistic, p, p.adj))

  # Combine formatting + P-value dataframe.
  stat.test.test <- p_vals |> 
    left_join(stat.test, by = c("channel", "cell_type", ".y.", "group1", "group2")) |> 
    mutate(y.position = max * 1.3)
    #mutate(y.position = ifelse(group2 == "mAAb+", y.position*1.05, y.position)) |>
    #mutate(y.position = ifelse(group2 == "RecentOnset", y.position*1.1, y.position))

  # Plot: 
  if (x_var == cur_var) {
    x_var <- "cur_var"
  }
  p <- ggviolin(cols_islet_plot, x = x_var, y = "MeanExprs", fill = "cur_var",
          facet.by = facet_by, palette = cur_palette, 
          draw_quantiles = c(0.25, 0.5, 0.75), add = "point", 
          add.params = list(alpha = 1, size = 1.5)) + 
    stat_pvalue_manual(stat.test.test, label = "p.signif", tip.length = 0.1,
                      size = 10, label.size = 10, vjust = 0.5) + 
    facet_wrap(~get(facet_by), scales = "free_y") +
    labs(x = "", y = y_label) + 
    mytheme$standard_new() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  return(p)
}

plot_clusters_if_score <- function(df, cur_clusters, ct_type, ncol = NULL, nrow = NULL) {
  ## Filter to results of interest.
  ## df <- score_clust[[1]]; ct_type <- "cluster"; cur_clusters <- c("T_CD8_1", "T_CD8_2", "T_CD8_3", "T_CD8_4", "T-naive", "T_CD8_6", "T-EMRA", "T_CD8_7")
  if (ct_type == "cell_type") {
    data_islet <- df |> 
      filter(CellType %in% cur_clusters) |> 
      dplyr::rename("immune_cluster" = "CellType") |> 
      inner_join(case_stages, by = c("case_id", "donor_type")) |>
      as_tibble()
  } else if (ct_type == "cluster") {
    data_islet <- df |> 
      filter(Cluster %in% cur_clusters) |> 
      dplyr::rename("immune_cluster" = "Cluster") |> 
      inner_join(case_stages, by = c("case_id", "donor_type")) |>
      as_tibble()
  } else {
    stop("ct_type must be either 'cell_type' or 'cluster'")
  }

  ## log-scale the score.
  data_islet <- data_islet |> 
    mutate(NormScoreArea = log1p(NormScoreArea))

  ## Calculate p-values for the boxplots.
  p_vals <- data_islet  |> 
    nest(data = -c(immune_cluster)) |> 
    mutate(test = purrr::map(data, \(.x) {
      lm_model <- lmerTest::lmer(NormScoreArea ~ donor_type + age_group + ICU_time_days + (1|case_id), data = .x)
      emmeans::emmeans(lm_model, ~ donor_type, pbkrtest.limit = 5000) |> 
        emmeans::contrast("pairwise", adjust = "none") |> 
        broom.mixed::tidy() |> filter(contrast != "(sAAb+) - RecentOnset")
    })) |>
    tidyr::unnest(test)
             
  p_vals <- p_vals |> 
    mutate(group1 = contrast, group2 = contrast) |>  
    mutate(group1 = stringr::str_replace(group1, "\\s.*", "")) |> 
    mutate(group2 = stringr::str_replace(group2, "^(?:\\S+\\s+\\S+\\s+)", "")) |> 
    mutate(group1 = stringr::str_remove(group1, "\\(")) |>
    mutate(group2 = stringr::str_remove(group2, "\\(")) |> 
    mutate(group1 = stringr::str_remove(group1, "\\)")) |> 
    mutate(group2 = stringr::str_remove(group2, "\\)"))

  p_vals <- p_vals |>
    mutate(donor_type = stringr::str_remove(term, "donor_type")) |> 
    group_by(immune_cluster) |>
    mutate(.y. = "NormScoreArea",
           p = p.value,
           p.adj = p.adjust(p, method = "fdr"),
           p.signif = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "•", ""))
           ) |> 
    select(.y., group1, group2, p, p.adj, p.signif, immune_cluster) |> 
    filter(p.adj < 0.05)
  
  ## Data to plot (by case).
  data_tt <- data_islet |> 
      group_by(immune_cluster, donor_type, case_id) |> 
      summarize(NormScoreArea = mean(NormScoreArea), .groups = "drop")
  
  # get max-y-value PER channel.
  p_vals <- data_tt |> 
      group_by(immune_cluster) |> 
      summarise(max = max(NormScoreArea), .groups = "drop") |>
      mutate(pos_pvals = max*1.05) |> 
      dplyr::right_join(p_vals, by = c("immune_cluster"))

  # ggpubr (for helping with formatting)
  stat.test <- data_tt  |> 
    group_by(immune_cluster) |> 
    tukey_hsd(NormScoreArea ~ donor_type) |> 
    filter(!(group1 == "sAAb+" & group2 == "RecentOnset")) |> 
    select(immune_cluster, group1, group2) |> 
    add_xy_position(x = "donor_type", step.increase = 0.07, dodge = 0.8)

  stat.test.test <- p_vals |> 
      left_join(stat.test, by = c("immune_cluster", "group1", "group2")) |> 
      select(-p) |> 
      mutate(y.position = max * 1.1)  |> 
      group_by(immune_cluster) %>%
      mutate(y.position = first(y.position) * (1 + (0.07) * (row_number() - 1))) %>%
      ungroup()
  
  # Boxplots:
  bxp <- ggboxplot(data_tt, x = "donor_type", y = "NormScoreArea", facet.by = "immune_cluster", scales = "free_y",
                  fill = "donor_type", add = "jitter",
                  add.params = list(size = 3, width = 0.2),
                  palette = palettes$stages, alpha = 1) + 
        labs(x = "Stage", y = "Infiltration Score")

  p <- bxp + 
    stat_pvalue_manual(stat.test.test, label = "p.signif", tip.length = 0.02,
                      size = 10, label.size = 10, 
                      step.increase = 0) + 
    mytheme$large_new() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) + 
    scale_fill_manual(values = palettes$stages) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

 return(p)
}

plot_abundance_stages <- function(stats, dat, cur_celltypes, 
                                  x, y, facet_by, fill_by, 
                                  ct_type = "cell_type", 
                                  nrow = NULL, parse_labels = FALSE,
                                  title = NULL, label_map = NULL,
                                  step.increase = 0.07) {
  require(ggpubr)
  if (ct_type == "cell_type") {
    df <- dat |> 
      filter(cell_type %in% cur_celltypes) |> 
      inner_join(case_stages2, by = c("case_id", "donor_type")) |>
      as_tibble()
  } else if (ct_type == "cluster") {
    df <- dat |> 
      filter(immune_cluster %in% cur_celltypes) |> 
      as_tibble() |> 
      dplyr::mutate(cell_type = immune_cluster)
  } else {
    stop("ct_type must be either 'cell_type' or 'cluster'")
  }

  ## Log-scale.
  df <- df |> 
    mutate(mean_density = log1p(mean_density))

  p_vals <- stats |> 
    filter(cell_type %in% cur_celltypes) |> 
    select(cell_type, group1, group2, PValue, FDR)

  # Use this for actual P-Values.
  # FDRs are as computed by edgeR.
  p_vals <- p_vals |> 
    mutate(.y. = "mean_density",
           p = PValue,
           statistic = F,
           p.adj = FDR,
           p.signif = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", "•", ""))
           ) |> 
    select(.y., group1, group2,
           statistic, p, p.adj, p.signif, cell_type)
  
  # get max-y-value PER channel.
  p_vals <- df |> 
      group_by(cell_type) |> 
      summarise(max = max(mean_density), .groups = "drop") |>
      mutate(pos_pvals = max*1.05) |> 
      dplyr::right_join(p_vals, by = c("cell_type"))

  ## Calculate p-values with ggpubr.
  # Use this for right formatting of DF
  stat.test <- df |> 
    mutate(donor_type = factor(donor_type, levels = c("NoDiabetes", "sAAb+", "mAAb+", "RecentOnset", "LongDuration"))) |>
    group_by(cell_type) |> 
    tukey_hsd(mean_density ~ donor_type)  |> 
    filter(!(group1 == "sAAb+" & group2 == "RecentOnset")) |> 
    filter(!(group1 == "sAAb+" & group2 == "LongDuration")) |>
    filter(!(group1 == "mAAb+" & group2 == "LongDuration")) |>
    filter(!(group2 == "sAAb+" & group1 == "RecentOnset")) |> 
    filter(!(group2 == "sAAb+" & group1 == "LongDuration")) |>
    filter(!(group2 == "mAAb+" & group1 == "LongDuration")) |>
    select(group1, group2, cell_type) |> 
    add_xy_position(x = "donor_type", step.increase = step.increase, dodge = 0.8)
  
  stat.test.test <- p_vals |> 
    left_join(stat.test, by = c("cell_type", "group1", "group2")) |> 
    select(-c(statistic, p)) |> 
    mutate(y.position = max * 1.1)  |> 
    filter(p.adj < 0.1) |> 
    group_by(cell_type) %>%
    mutate(y.position = first(y.position) * (1 + (step.increase) * (row_number() - 1))) %>%
    ungroup()

  # Plot.
  # Relevel + relabel the facet column using your map
  if (parse_labels & !is.null(label_map)){
    df[[facet_by]] <- factor(df[[facet_by]],
                          levels = label_map)
    stat.test.test[[facet_by]] <- factor(stat.test.test[[facet_by]],
                                        levels = label_map)
  }

  p <- ggboxplot(df, x = x, y = y, 
           fill = fill_by, palette = NULL, outliers = FALSE, 
           add = "jitter", add.params = list(size = 3, width = 0.2)) +
    labs(x = "Stage", y = "log1p: Cell Density (cells/mm^2)") +
    # mytheme$large_new() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = palettes$stages) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    stat_pvalue_manual(stat.test.test, label = "p.signif", tip.length = 0.02,
                       size = 10, label.size = 10, inherit.aes = FALSE) +
    ggtitle(title)
  
  if (parse_labels){
    message("parsing labels ...")
    p <- p + facet_wrap(
        vars(.data[[facet_by]]), 
        scales = "free_y", 
        axes = "all", 
        nrow=nrow,
        #labeller = if (isTRUE(parse_labels)) lab else label_value)
          labeller = as_labeller(label_map, default = label_parsed))
  } else {
    p <- p + facet_wrap(~get(facet_by), scales = "free_y", axes = "all", nrow=nrow)
  }
   
  return(p)
}

#' plot_density
#' Plot Density plots
#'
#' @param dat A `data.frame`.
#' @param x x axis variable.
#' @param y y axis variable.
#' @param fill_by Variable to fill by.
#' @param color_by (optional) Variable to color by.
#' @param facet_by (optional) Variable to facet by.
#' @param ncol number of column in facetted plot.
#' @param title Graph title.
#'
#' @return a `ggplot` object
#' @export

plot_density <- function(dat, x, y, fill_by, color_by = NULL,
                         facet_by = NULL, ncol = NULL, title = NULL) {
  require(ggridges)

  p <- ggplot2::ggplot(dat, aes(x = get(x), y = get(y))) +
    scale_fill_viridis() +
    mytheme$standard() +
    theme(
      panel.background = element_rect(colour = "grey20")
    )

  if (is.null(color_by))
    p <- p + geom_density_ridges(ggplot2::aes(fill = get(fill_by)),
                                 alpha = 0.5, lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by)
  else
    p <- p + geom_density_ridges(ggplot2::aes(fill = get(fill_by),
                                              color = get(color_by)),
                                 alpha = 0.5, lwd = 0.3) +
      ggplot2::labs(x = x, y = y, fill = fill_by, color = color_by)

  if (!is.null(facet_by))
    p <- p + facet_wrap(~get(facet_by), ncol = ncol)

  if (is.null(title))
    p <- p + ggtitle(paste(x, y, sep = " - "))
  else
    p <- p + ggtitle(title)

  return(p)
}

#' calc_compo
#' Calculate cell type composition and cell density (number of cells per mm^2).
#' This function can be used with one cluster and one or two grouping variables.
#' It returns the following results:
#' For fractions: out of all the cells in the defined Cluster (`clust`),
#' what is the fraction of cells that belong to each group (`group1`, one grouping variable)
#' or subgroup (`group1` X `group2`, two grouping variables) ?
#' For density: in the area defined by `group1` (one grouping variable) or
#' `group1` X `group2` (two grouping variables), what is the cell density in
#' the defined Cluster (`clust`) ?
#'
#' @param x A `SingleCellExperiment` or `SpatialExperiment` object.
#' @param clust the first cell type or cluster column to use.
#' @param orderclust a vector that indicates how variables in `clust` should be ordered.
#' @param group1 first variable to group by.
#' @param order1 a vector that indicates how variables in `group1` should be ordered.
#' @param group2 (optional) second variable to group by.
#' @param order2 (optional) a vector that indicates how variables in `group2` should be ordered.
#'
#' @return a `tibble` data frame.
#' @export

calc_compo <- function(x, clust, orderclust, group1, order1,
                       group2 = NULL, order2 = NULL) {
  name_clust <- "CellsPerCluster"

  # ONE GROUPING VARIABLE
  if (is.null(group2)) {
    # Calculate total cell number per cluster
    total_cells <- as_tibble(SummarizedExperiment::colData(x)) |>
      dplyr::add_count(get(clust), name = name_clust) |>
      select(c(Cluster = all_of(clust), all_of(name_clust))) |>
      distinct()

    # Calculate number of cells of each cluster per group
    compo <- as_tibble(SummarizedExperiment::colData(x)) |>
      group_by(Group1 = get(group1), Cluster = get(clust)) |>
      dplyr::count(name = "CellNb") |>
      ungroup() |>
      tidyr::complete(Group1, Cluster, fill = list(CellNb = 0))

    # Calculate cluster fraction per group
    compo <- compo |>
      inner_join(total_cells, by = "Cluster") |>
      filter(get(name_clust) > 0) |>
      mutate(fraction = CellNb / get(name_clust))

    # Calculate image area as the sum of cell areas
    area <- as_tibble(SummarizedExperiment::colData(x)) |>
      group_by(Group1 = get(group1)) |>
      summarise(Area = sum(cell_area) / 10^6)

    # Calculate cell density (cell number per mm^2)
    compo <- compo |>
      inner_join(area, by = "Group1") |>
      mutate(density = CellNb / Area)

    # Re-order the dataset
    compo <- compo |>
      mutate(Group1 = factor(Group1, levels = order1),
             Cluster = factor(Cluster, levels = rev(orderclust))) |>
      arrange(Group1, Cluster)

    return(compo |>
             filter(!is.na(Cluster) & !is.na(Group1)))
  } else { # TWO GROUPING VARIABLES
    if (is.null(order2))
      stop("Provide an ordering vector for the second grouping variable")

    # Calculate total cell number per cluster and grouping variable 2
    total_cells <- as_tibble(SummarizedExperiment::colData(x)) |>
      dplyr::add_count(get(clust), get(group2), name = name_clust) |>
      select(c(Cluster = all_of(clust), Group2 = all_of(group2),
               all_of(name_clust))) |>
      distinct()

    # Calculate number of cells of each cluster per group
    compo <- as_tibble(SummarizedExperiment::colData(x)) |>
      group_by(Group1 = get(group1),
               Group2 = get(group2),
               Cluster = get(clust)) |>
      dplyr::count(name = "CellNb") |>
      ungroup() |>
      tidyr::complete(Group1, Group2, Cluster, fill = list(CellNb = 0))

    # Calculate cluster fraction per group
    compo <- compo |>
      inner_join(total_cells, by = c("Cluster", "Group2")) |>
      filter(get(name_clust) > 0) |>
      mutate(fraction = CellNb / get(name_clust))

    # Calculate image area
    area <- as_tibble(SummarizedExperiment::colData(x)) |>
      group_by(Group1 = get(group1), Group2 = get(group2)) |>
      summarise(Area = sum(cell_area) / 10^6, .groups = "keep")

    # Calculate cell density (cell number per mm^2)
    compo <- compo |>
      inner_join(area, by = c("Group1", "Group2")) |>
      mutate(density = CellNb / Area)

    # Re-order the dataset
    compo <- compo |>
      mutate(Group1 = factor(Group1, levels = order1),
             Group2 = factor(Group2, levels = order2),
             Cluster = factor(Cluster, levels = rev(orderclust))) |>
      arrange(Group1, Group2, Cluster)

    compo <- compo |>
      filter(!is.na(Cluster) & !is.na(Group1) & !is.na(Group2))
  }
  return(compo)
}


#' YOLOader
#' Wrapper to load images and masks as `CytoImageList` objects for plotting with `cytomapper`.
#' No checks whatsoever so make sure the arguments and outputs are right.
#'
#' @param x a `SingleCellExperiment` or `SpatialExperiment` object.
#' @param image_dir directory containing the images to load
#' @param image_names names of the images to load.
#' @param suffix_rem suffix to remove from image names.
#' @param suffix_add suffix to add to image names.
#' @param bit.depth image bit depth.
#' @param type either "stacks" or "masks", depending on the object to return.
#' @param ... additional arguments to pass to the cytomapper::loadImages() function.
#' @return a `CytoImageList` object, containin either image stacks or masks.

yoloader <- function(x, image_dir, image_names,
                     suffix_rem = "", suffix_add = "",
                     bit_depth = 16, type, ...) {
  require(cytomapper)

  image_list <- file.path(image_dir, image_names)

  # Test if the image list exist
  test_exist <- which(!file.exists(image_list))
  if (length(test_exist) > 0) {
    stop(c("The following images were not found:\n",
           paste(image_list[test_exist], collapse = "\n")))
  } else {
    # Load and scale the images
    images <- loadImages(image_list, ...)
    # images <- scaleImages(images, (2 ^ bit.depth) - 1)

    # Add image names to metadata
    mcols(images)$ImageName <- gsub(suffix_rem, "", names(images))
    mcols(images)$ImageName <- paste0(mcols(images)$ImageName, suffix_add)

    # Add channel names
    if (type == "stacks") {
      channelNames(images) <- rownames(x)
    }

    return(images)
  }
}


#' Plotting themes for ggplot2
mytheme <- list(
  standard = function(base_size = 16, base_family = "Arial") {
    theme(
      plot.title = element_text(size = 24),

      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),

      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.line  = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),

      strip.text = element_text(size = 12),
      strip.background = element_blank(),

      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0.1, "line")
    )
  },
  standard_new = function(base_size = 22, base_family = "Arial") {
    theme(
      plot.title = element_text(size = 28),

      legend.title = element_text(size = 25),
      legend.text = element_text(size = 22),

      axis.title = element_text(size = 22),
      axis.text = element_text(size = 22),
      axis.line  = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),

      strip.text = element_text(size = 25),
      strip.background = element_blank(),

      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0.3, "line")
    )
  },
  large = function(base_size = 24, base_family = "Arial") {
    theme(
      plot.title = element_text(size = 32),

      legend.title = element_text(size = 24),
      legend.text = element_text(size = 24),

      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      axis.line  = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),

      strip.text = element_text(size = 20),
      strip.background = element_blank(),

      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0.5, "line")
    )
  },
  large_new = function(base_size = 32, base_family = "Arial") {
    theme(
      plot.title = element_text(size = 37),

      legend.title = element_text(size = 35),
      legend.text = element_text(size = 32),

      axis.title = element_text(size = 35),
      axis.text = element_text(size = 32),
      axis.line  = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),

      strip.text = element_text(size = 32),
      strip.background = element_blank(),

      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0.5, "line")
    )
  }
)


#' Color palettes for plotting
palettes <- list(
  colors = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
             "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
             "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
             "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
             "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
             "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"),

  colors50 = c("#C4625D", "#DE7390", "#B35A77", "#87638F", "#B65E46",
               "#3A7BB4", "#CDA12C", "#D46A42", "#93539D", "#56A354",
               "#F67A0D", "#C78DAB", "#D06C78", "#A9572E", "#B06A29",
               "#4CAD4E", "#419584", "#BF862B", "#735B81", "#449D72",
               "#7A7380", "#8F4A68", "#FFC81D", "#566B9B", "#48A460",
               "#999999", "#FFDD25", "#EB7AA9", "#E585B8", "#F9F432",
               "#AB3A4E", "#3A85A8", "#E41A1C", "#3D8D96", "#E57227",
               "#B791A5", "#629363", "#C72A35", "#FFF12D", "#F581BE",
               "#6E8371", "#D689B1", "#FF9E0C", "#C3655F", "#EBD930",
               "#DCBD2E", "#A25392", "#FF8904", "#FFB314", "#A8959F"),
  
  colorsfull = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = TRUE)],

  stages = c(
    "NoDiabetes" = "#3182BD", "sAAb+" = "#FED976", "mAAb+" = "#FD8D3C",
    "RecentOnset" = "#E31A1C", "LongDuration" = "#800026",
    "Pancreatitis" = "#252525"
  ),

  cases = c(
    "6034" = "#ff996c", "6048" = "#2040d2", "6055" = "#00a824",
    "6126" = "#ae25c6", "6178" = "#00c751", "6229" = "#8a37d2",
    "6234" = "#0a8c00", "6289" = "#c356f2", "6420" = "#6ece53",
    "6428" = "#7e49e3", "6447" = "#8eae00", "6488" = "#5d59f3",
    "6509" = "#b5b200", "6531" = "#0061ef", "6547" = "#639900",
    "6044" = "#0150d7", "6090" = "#01c25f", "6123" = "#f42ac3",
    "6147" = "#007611", "6151" = "#ff60ed", "6171" = "#5acf7b",
    "6181" = "#a800a5", "6301" = "#8cc869", "6303" = "#732cae",
    "6310" = "#d3b733", "6314" = "#017bfe", "6347" = "#ff8612",
    "6397" = "#0159cc", "6400" = "#e58d00", "6421" = "#8077ff",
    "6433" = "#f4aa2f", "6437" = "#543cb7", "6483" = "#a88f00",
    "6494" = "#b87cff", "6496" = "#667c00", "6517" = "#e07bff",
    "6518" = "#019c59", "6532" = "#e300a6", "6538" = "#458335",
    "6553" = "#ff5bcb", "6558" = "#adc259", "6562" = "#234bb0",
    "8011" = "#ff7125", "6197" = "#017eed", "6388" = "#ff4c21",
    "6424" = "#38a6ff", "6429" = "#f72e28", "6450" = "#6688ff",
    "6505" = "#d15700", "6512" = "#a88fff", "6521" = "#e9ae41",
    "6549" = "#5841a7", "8002" = "#cc4000", "6209" = "#9e9eff",
    "6224" = "#bc2100", "6228" = "#7e74cb", "6247" = "#b40f00",
    "6324" = "#725bae", "6362" = "#a06a00", "6367" = "#ff89f9",
    "6380" = "#ad4900", "6396" = "#5e4395", "6405" = "#ff8f4b",
    "6414" = "#703994", "6449" = "#fea35d", "6456" = "#89268a",
    "6469" = "#f4a861", "6520" = "#a30083", "6526" = "#b92e00",
    "6533" = "#fa94f5", "6534" = "#953002", "6550" = "#dc95e7",
    "6551" = "#d30026", "6563" = "#a16cbe", "8005" = "#ff6a4d",
    "8008" = "#83327a", "6063" = "#ff234b", "6135" = "#d26aa8",
    "6145" = "#a01f1b", "6208" = "#ff6ec7", "6264" = "#ff876b",
    "6285" = "#a00367", "6321" = "#ff625b", "6328" = "#9d326a",
    "6418" = "#cb0036", "6422" = "#ff7eb4", "6458" = "#9d3432",
    "6510" = "#f20082", "6514" = "#d56959", "6519" = "#d5007d",
    "6036" = "#b74d44", "6043" = "#ff599b", "6061" = "#a70431",
    "6150" = "#ff5989", "6225" = "#b6436c", "6506" = "#eb0054",
    "6522" = "#c2016a"
  ),

  casestages = c(
    "6034" = "#3182BD", "6048" = "#3182BD", "6055" = "#3182BD",
    "6126" = "#3182BD", "6178" = "#3182BD", "6229" = "#3182BD",
    "6234" = "#3182BD", "6289" = "#3182BD", "6420" = "#3182BD",
    "6428" = "#3182BD", "6447" = "#3182BD", "6488" = "#3182BD",
    "6509" = "#3182BD", "6531" = "#3182BD", "6547" = "#3182BD",
    "6044" = "#FED976", "6090" = "#FED976", "6123" = "#FED976",
    "6147" = "#FED976", "6151" = "#FED976", "6171" = "#FED976",
    "6181" = "#FED976", "6301" = "#FED976", "6303" = "#FED976",
    "6310" = "#FED976", "6314" = "#FED976", "6347" = "#FED976",
    "6397" = "#FED976", "6400" = "#FED976", "6421" = "#FED976",
    "6433" = "#FED976", "6437" = "#FED976", "6483" = "#FED976",
    "6494" = "#FED976", "6496" = "#FED976", "6517" = "#FED976",
    "6518" = "#FED976", "6532" = "#FED976", "6538" = "#FED976",
    "6553" = "#FED976", "6558" = "#FED976", "6562" = "#FED976",
    "8011" = "#FED976", "6197" = "#FD8D3C", "6388" = "#FD8D3C",
    "6424" = "#FD8D3C", "6429" = "#FD8D3C", "6450" = "#FD8D3C",
    "6505" = "#FD8D3C", "6512" = "#FD8D3C", "6521" = "#FD8D3C",
    "6549" = "#FD8D3C", "8002" = "#FD8D3C", "6209" = "#E31A1C",
    "6224" = "#E31A1C", "6228" = "#E31A1C", "6247" = "#E31A1C",
    "6324" = "#E31A1C", "6362" = "#E31A1C", "6367" = "#E31A1C",
    "6380" = "#E31A1C", "6396" = "#E31A1C", "6405" = "#E31A1C",
    "6414" = "#E31A1C", "6449" = "#E31A1C", "6456" = "#E31A1C",
    "6469" = "#E31A1C", "6520" = "#E31A1C", "6526" = "#E31A1C",
    "6533" = "#E31A1C", "6534" = "#E31A1C", "6550" = "#E31A1C",
    "6551" = "#E31A1C", "6563" = "#E31A1C", "8005" = "#E31A1C",
    "8008" = "#E31A1C", "6063" = "#800026", "6135" = "#800026",
    "6145" = "#800026", "6208" = "#800026", "6264" = "#800026",
    "6285" = "#800026", "6321" = "#800026", "6328" = "#800026",
    "6418" = "#800026", "6422" = "#800026", "6458" = "#800026",
    "6510" = "#800026", "6514" = "#800026", "6519" = "#800026",
    "6036" = "#252525", "6043" = "#252525", "6061" = "#252525",
    "6150" = "#252525", "6225" = "#252525", "6506" = "#252525",
    "6522" = "#252525"
  ),

  batch = c(
    "#FB8072", "#7BAFDE", "#B17BA6", "#FDB462"
  )
)

#' Parameters to save plots
plotsave_param <- list(
  device = "png", units = "mm",
  width = 300, height = 200, dpi = 300
)

plotsave_param_large <- list(
  device = "png", units = "mm",
  width = 600, height = 400, dpi = 300
)


# Myeloid.
myeloid_labels <- c(
  "Myeloid_1" = "Exocrine low act.",
  "Myeloid_2" = "M1/M2-like low act.",
  "Myeloid_3" = "M1/M2-like act.",
  "Myeloid_4" = "cDC",
  "Myeloid_5" = "CD11c+",
  "Myeloid_6" = "Exocrine act.",
  "Myeloid_7" = "Exocrine MPO+",
  "Myeloid_8" = "CD54+",
  "Myeloid_Other" = "Ambiguous"
)

# T_CD4.
tcd4_labels <- c(
  "T_CD4_1" = "CD4+ Other",
  "T_CD4_2" = "T-EM-like low act.",
  "T_CD4_3" = "PD1+ low act.",
  "T_CD4_4" = "T-CM-like act.",
  "T_CD4_5" = "PD1+ act.",
  "T_CD4_6" = "T-EM-like act.",
  "T_CD4_7" = "PD1low act.", # expression("PD1"^"low"~"act."),
  "T_CD4_8" = "PD1low low act.", # expression("PD1"^"low"~"low act."),
  "T_CD4_9" = "CD73+",
  "T_CD4_10" = "T-reg",
  "T_CD4_Other" = "Ambiguous"
)

tcd4_labels3 <- c(
  "T_CD4_1" = "'T-CD4 Other'",
  "T_CD4_2" = "'T-EM-like low act.'",
  "T_CD4_3" = "PD1*'+'~low~act.", 
  "T_CD4_4" = "'T-CM-like act.'",
  "T_CD4_5" = "PD1*'+'~act.",
  "T_CD4_6" = "'T-EM-like act.'",
  "T_CD4_7" = "PD1^low~act.", # expression("PD1"^"low"~"act."),
  "T_CD4_8" = "PD1^low~low~act.", # expression("PD1"^"low"~"low act."),
  "T_CD4_9" = "'CD73+'",
  "T_CD4_10" = "'T-reg'",
  "T_CD4_Other" = "Ambiguous"
)

# T-CD8.
tcd8_labels <- c(
  "T_CD8_1" = "T-CD8 Other",
  "T_CD8_2" = "T-RM low act.",
  "T_CD8_3" = "GranB+",
  "T_CD8_4" = "T-EM/CM-like",
  "T-naive" = "T-naive",
  "T_CD8_6" = "T-RM act.",
  "T_CD8_7" = "T-exeff",
  "T-EMRA" = "T-EMRA"
)

tcd8_labels3 <- c(
  "T_CD8_1" = "'T-CD8 Other'",
  "T_CD8_2" = "'T-RM low act.'",
  "T_CD8_3" = "'GranB+'", 
  "T_CD8_4" = "'T-EM/CM-like'",
  "T-naive" = "'T-naive'",
  "T_CD8_6" = "'T-RM act.'",
  "T_CD8_7" = "'T-ex'^eff",
  "T-EMRA" = "'T-EMRA'"
)

## Labels for the cell types.
all_labels <- c(
  "Myeloid_1" = "'M: Exocrine low act.'",
  "Myeloid_2" = "'M: M1/M2-like low act.'",
  "Myeloid_3" = "'M: M1/M2-like act.'",
  "Myeloid_4" = "'cDC'",
  "Myeloid_5" = "'M: CD11c+'",
  "Myeloid_6" = "'M: Exocrine act.'",
  "Myeloid_7" = "'M: Exocrine MPO+'",
  "Myeloid_8" = "'M: CD54+'",
  "T_CD4_1" = "'T-CD4: T-CD4 Other'",
  "T_CD4_2" = "'T-CD4: T-EM-like low act.'",
  "T_CD4_3" = "'T-CD4: PD1+ low act.'",
  "T_CD4_4" = "'T-CD4: T-CM-like act.'",
  "T_CD4_5" = "'T-CD4: PD1+ act.'",
  "T_CD4_6" = "'T-CD4: T-EM-like act.'",
  "T_CD4_7" = "'T-CD4: PD1'^low~act.", # expression("PD1"^"low"~"act."),
  "T_CD4_8" = "'T-CD4: PD1'^low~low~act.", # expression("PD1"^"low"~"low act."),
  "T_CD4_9" = "'T-CD4: CD73+'",
  "T_CD4_10" = "'T-reg'",
  "T_CD8_1" = "'T-CD8: T-CD8 Other'",
  "T_CD8_2" = "'T-CD8: T-RM low act.'",
  "T_CD8_3" = "'T-CD8: GranB+'",
  "T_CD8_4" = "'T-CD8: T-EM/CM-like'",
  "T-naive" = "'T-CD8: T-naive'",
  "T_CD8_6" = "'T-CD8: T-RM act.'",
  "T_CD8_7" = "'T-CD8: T-ex'^eff",
  "T-EMRA" = "'T-CD8: T-EMRA'",
  "B" = "'B-cell'"
)


## Labels for the cell types.
all_labels2 <- c(
  "B" = "'B-cell'",
  "T-naive" = "'T-CD8: T-naive'",
  "Myeloid_3" = "'M: M1/M2-like act.'",
  "T_CD8_7" = "'T-CD8: T-ex'^eff",
  "Myeloid_5" = "'M: CD11c+'",
  "T_CD8_2" = "'T-CD8: T-RM low act.'",
  "Myeloid_4" = "'cDC'",
  "T_CD4_5" = "'T-CD4: PD1+ act.'"
)

all_labels2 <- factor(all_labels2, levels = all_labels2)

CN_labels <- c("B-cell" = "B-cell rich/TLS-like",
               "T-cell 1" = "T-CD8 > T-CD4",
                "T-cell 2" = "Exocrine T-CD8",
                "Exocrine Inf" = "Exocrine inflammation",
                "Mixed Inf" = "Innate inflammation",
                "Mixed Inf mAAB+" = "Innate inflammation mAAb+",
                "Exocrine CTRL" = "Exocrine Control",
                "Myeloid Inf" = "Myeloid inflammation",
                "Myeloid Inf mAAB+" = "Myeloid inflammation mAAb+",
                "Myeloid Inf Onset" = "Myeloid inflammation Onset",
                "T-cell 3" = "T-CD4 > T-CD8",
                "Beta" = "Beta",
                "Alpha" = "Alpha",
                "Delta" = "Delta",
                "Islet-Edge CTRL" = "Islet-edge Control",
                "Islet-Edge mAAB+" = "Islet-edge mAAb+",
                "Islet-Edge Onset" = "Islet-edge Onset",
                "Islet-Endothelial" = "Islet-endothelial"
)   

CN_labels_full <- c(CN_labels,
                "Neutrophil" = "Neutrophil",
                "Nerves" = "Nerves",
                "Stroma" = "Stroma",
                "CD303_VIM" = "Stroma CD303+/VIM+",
                "Stroma CTRL" = "Exocrine/Stroma",
                "Ductal" = "Ductal",
                "Acinar" = "Acinar",
                "Endothelial" = "Endothelial",
                "Smooth-Muscle" = "Smooth-Muscle",
                "Mixed-Islet" = "Mixed-islet",
                "Islet-Other" = "Islet-other"
                ) 

major_labels <- c("T_CD4" = "T-CD4", 
                  "T_CD8" = "T-CD8", 
                  "Myeloid" = "Myeloid", 
                  "Neutrophil" = "Neutrophil")

immune_labels <- c("T_CD4" = "T-CD4", 
                  "T_CD8" = "T-CD8", 
                  "B" = "B-cell", 
                  "NK" = "NK",
                  "Neutrophil" = "Neutrophil", 
                  "Myeloid" = "Myeloid", 
                  "T_DN" = "T-DN", 
                  "CD303_VIM" = "CD303+/VIM+")

full_immune_labels <- c(immune_labels, 
  "Nerves" = "Nerves",
  "Stroma" = "Stroma",
  "Ductal" = "Ductal",
  "Acinar" = "Acinar",
  "Endothelial" = "Endothelial",
  "SmoothMuscle" = "Smooth-Muscle",
  "Mesenchymal" = "Mesenchymal",
  "Beta" = "Beta",
  "Alpha" = "Alpha",
  "Delta" = "Delta",
  "IsletOther" = "Islet-other")

breaks <- c(0, 0.5, 0.7, 1)  # Define custom breaks for the color scale
colors <- viridis::viridis(4) 
# Create a color function with more emphasis on the lower values
col_fun <- circlize::colorRamp2(breaks, colors)

# INsulitic ROIs:
palettes$insulitic <- c("Insulitic" = "#E69F00", "Non-Insulitic" = "#009E73")
