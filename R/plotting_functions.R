#' Plot of posterior fit
#'
#' @description
#' `plot_fit()` plots the posterior fit for the desired observations, overlain on the true data.
#'
#' @param t Time points for the x axis. Either a vector that is common for all observations, or a matrix for which
#' each row corresponds to an observation.
#' @param y Matrix of data. Each row corresponds to one observation.
#' @param fit Matrix of fit. Each row corresponds to one observation.
#' @param indices Vector indicating observations for which to plot the fit. Default is `sample(1:nrow(y), 5)`.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"y"`.
#' @param title Title for the figure. Default is `""`.
#' @param line_colors Vectors of colors for the lines of fit. Default assumes 5 lines are plotted, and is
#' `c("red", "blue", "green", "purple", "orange")`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_fit <- function(t, y, fit, indices = sample(1:nrow(y), 5),
                     xlab = "x", ylab =  "y", title = "",
                     line_colors = c("red", "blue", "green", "purple", "orange")) {

  N <- nrow(y)
  n <- ncol(y)

  num_indices <- length(indices)

  if (length(line_colors) != num_indices) {
    stop("Provide line colors of the same length as the number of fitted lines to be plotted.")
  }

  if (is.matrix(t)){
    if (nrow(t) == 1 | ncol(t) == 1) {
      t <- as.numeric(t)
      if (length(t) != n) {
        stop("The dimesions of 't' and 'y' don't match.")
      }
      x <- rep(t, times = N)
      x_ind <- rep(t, times = num_indices)
    } else {
      if (nrow(t) != N | ncol(t) != n) {
        stop("The dimesions of 't' and 'y' don't match.")
      }
      x <- c(t(t))
      x_ind <- matrix(t[indices, ], nrow = n * num_indices, byrow = TRUE)
      x_ind <- as.numeric(x_ind)
    }
  } else {
    if (length(t) != n) {
      stop("The dimesions of 't' and 'y' don't match.")
    }
    x <- rep(t, times = N)
    x_ind <- rep(t, times = num_indices)
  }

  background_df <- data.frame(x = x,
                              y = as.vector(y),
                              group = "Background")

  index_df <- data.frame(x = x_ind,
                         y = as.vector(t(y[indices,])),
                         group = factor(rep(indices, each = n)))

  lines_df <- data.frame(x = x_ind,
                         y = as.vector(t(fit[indices, ])),
                         group = factor(rep(indices, each = n)))

  p <- ggplot2::ggplot() +
       ggplot2::geom_point(data = background_df, ggplot2::aes(x = .data$x, y = .data$y),
                           shape = 4, size = 2, color = "grey80") +
       ggplot2::geom_point(data = index_df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group),
                           shape = 1, size = 2) +
       ggplot2::geom_line(data = lines_df,
                          ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, color =.data$group),
                          size = 0.7) +
       ggplot2::scale_color_manual(values = line_colors) +
       ggplot2::labs(x = xlab, y = ylab, title = title ) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")

  return(p)
}


#' Plot of posterior feature shape
#'
#' @description
#' `plot_feature()` plots the posterior shape function for the desired feature. The posterior quantiles of interest
#' are also represented.
#'
#' @param eval_t Time points for the horizontal axis.
#' @param shape Posterior central summary of the shape function.
#' @param quantile_low Posterior lower quantile of the shape function.
#' @param quantile_high Posterior upper quantile of the shape function.
#' @param color Line color for the shape function. Default is `"black"`.
#' @param alpha Opacity of the color for the quantile ribbon. Default is `0.2`.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"f"'`.
#' @param title Title for the figure. Default is `""`.
#' @param background_ind Indices of observations to plot. Default is `NULL`, in which case no observations are
#' displayed with the posterior shape.
#' @param t Time points at which data is evaluated. Default is `NULL`, for the case when no observations
#' are displayed with the posterior shape.
#' @param y Matrix of data. Each row corresponds to one observation. Default is `NULL`, for the case when no observations
#' are displayed with the posterior shape.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_feature <- function(eval_t, shape, quantile_low, quantile_high, color = "black", alpha = 0.2,
                         xlab = "x", ylab =  "f", title = "",
                         background_ind = NULL, t = NULL, y = NULL) {

  n <- length(eval_t)
  t1 <- eval_t[1]
  tn <- eval_t[n]

  plot_df <- data.frame(x = eval_t,
                        f = shape,
                        f_low = quantile_low,
                        f_high = quantile_high)

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data$x, y = .data$f)) +
       ggplot2::geom_line(color = color) +
       ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$f_low, ymax = .data$f_high), fill = color, alpha = alpha) +
       ggplot2::labs(x = xlab, y = ylab, title = title) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0)),
                      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(r = 0))) +
       ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

  if (!is.null(background_ind)) {
    for (i in 1:length(background_ind)){
      line_data <- data.frame(x = t, y = y[background_ind[i], ])
      p <- p + ggplot2::geom_line(data = line_data,
                                  ggplot2::aes(x = .data$x, y = .data$y),
                                  color = "lightgrey",
                                  linetype = "dashed")
    }
  }

  return(p)
}


#' Plot of posterior time transformation functions
#'
#' @description
#' `plot_tt()` plots time-transformation functions. Options allow to plot all together or split into two groups, as
#' well as color based on covariate value.
#'
#'
#' @param eval_t Time points for the horizontal axis.
#' @param tt_functions Matrix of time transformation functions, where each row corresponds to one function.
#' @param group Vector of 0 and 1s indicating group belonging for each function. Default is `NULL`, for when the group
#' is not provided and all time transformation functions are plotted together.
#' @param covariate Vector of covariate values associated with each function. Default is `NULL`, for when covariate
#' value is not provided.
#' @param group_colors Line colors for each group. Default is `c("gold1", "cornflowerblue")`.
#' @param group_labels Labels to use to denote each group. Default is `c("0", "1")`.
#' @param covariate_gradient Two colors determining the gradient for the covariate values. Default is
#' `c("yellow", "blue")`.
#' @param xlab Label for the horizontal axis. Default is `"Physical Time"`.
#' @param ylab Label for the vertical axis. Default is `"Stochastic Time"`.
#' @param title Title for the figure. Default is `""`.
#' @param covlab Title for the legend for covariate values. Default is `"Age"`.
#'
#' @return
#' A ggplot figure.
#'
#' @export
#'
#' @importFrom rlang .data
plot_tt <- function(eval_t, tt_functions, group = NULL, covariate = NULL,
                    group_colors = c("gold1", "cornflowerblue"),
                    group_labels = c("0", "1"),
                    covariate_gradient = c("yellow", "blue"),
                    xlab = "Physical Time", ylab = "Stochastic Time", title = "", covlab = "Age") {

  N <- nrow(tt_functions)
  n <- length(eval_t)
  t1 <- eval_t[1]
  tn <- eval_t[n]

  if (is.null(group)){
    plot_df <- data.frame(x = eval_t,
                          y = t(tt_functions))
    long_df <- tidyr::pivot_longer(plot_df, cols = -.data$x, names_to = "Row_Index", values_to = "Values")

    p <- ggplot2::ggplot(data = long_df, ggplot2::aes(x = .data$x,
                                                        y = .data$Values,
                                                        group = .data$Row_Index,
                                                        color = factor(.data$Row_Index))) +
         ggplot2::geom_line(linewidth = 0.2) +
         ggplot2::labs(x = xlab, y = ylab, title = title) +
         ggplot2::scale_color_viridis_d(option = "A", direction = 1) +
         ggplot2::theme_classic() +
         ggplot2::theme(legend.position = "none",
                        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                        axis.line.x = ggplot2::element_line(color = "black"),
                        axis.line.y = ggplot2::element_line(color = "black"),
                        panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                        panel.background = ggplot2::element_blank(),
                        plot.margin = ggplot2::margin(10, 10, 10, 10)) +
         ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
         ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(t1, tn))
  } else {
    if (is.null(covariate)) {
      plot_df <- data.frame(x = rep(eval_t, N),
                            y = c(t(tt_functions)),
                            group = rep(1:N, each = length(eval_t)),
                            vec = rep(group, each = length(eval_t)))

      p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data$x,
                                                        y = .data$y,
                                                        group = .data$group,
                                                        color = factor(.data$vec))) +
           ggplot2::geom_line() +
           ggplot2::labs(x = xlab, y = ylab, title = title, color = "Group") +
           ggplot2::scale_color_manual(labels = stats::setNames(group_labels, as.character(unique(plot_df$vec))),
                                       values = stats::setNames(group_colors, as.character(unique(plot_df$vec)))) +
           ggplot2::theme_classic() +
           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                          axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                          axis.line.x = ggplot2::element_line(color = "black"),
                          axis.line.y = ggplot2::element_line(color = "black"),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank(),
                          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                          panel.background = ggplot2::element_blank(),
                          plot.margin = ggplot2::margin(10, 10, 10, 10)) +
           ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
           ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(t1, tn))
    } else {
      plot_df <- data.frame(x = rep(eval_t, N),
                            y = c(t(tt_functions)),
                            group = rep(1:N, each = length(eval_t)),
                            vec = rep(group, each = length(eval_t)),
                            cov = rep(covariate, each = length(eval_t)))

      df_0 <- plot_df[plot_df$vec == 0, ]
      df_1 <- plot_df[plot_df$vec == 1, ]

      p0 <- ggplot2::ggplot(data = df_0, ggplot2::aes(x = .data$x,
                                                      y = .data$y,
                                                      group = .data$group,
                                                      color = .data$cov)) +
            ggplot2::geom_line() +
            ggplot2::labs(x = xlab, y = ylab, color = covlab) +
            ggplot2::scale_color_gradient(low = covariate_gradient[1], high = covariate_gradient[2]) +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(group_labels[1]) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                           axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                           axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                           axis.line.x = ggplot2::element_line(color = "black"),
                           axis.line.y = ggplot2::element_line(color = "black"),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                           panel.background = ggplot2::element_blank(),
                           plot.margin = ggplot2::margin(10, 10, 10, 10)) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
            ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(t1, tn))

      p1 <- ggplot2::ggplot(data = df_1, ggplot2::aes(x = .data$x,
                                                      y = .data$y,
                                                      group = .data$group,
                                                      color = .data$cov)) +
        ggplot2::geom_line() +
        ggplot2::labs(x = xlab, y = ylab, color = covlab) +
        ggplot2::scale_color_gradient(low = covariate_gradient[1], high = covariate_gradient[2]) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(group_labels[2]) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                       axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                       axis.line.x = ggplot2::element_line(color = "black"),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                       panel.background = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(10, 10, 10, 10)) +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(t1, tn))

      p <- gridExtra::grid.arrange(p0, p1, ncol = 2,
                                   top = grid::textGrob(title,
                                                             gp = grid::gpar(fontsize =15, font = 2)))
    }
  }
  return(p)
}

#' Plot of the data
#'
#' @description
#' `plot_data()` plots the rows of the matrix of observations.
#'
#' @param t  Time points for the x axis. Either a vector that is common for all observations, or a matrix for which
#' each row corresponds to an observation.
#' @param y Matrix of data. Each row corresponds to one observation.
#' @param indices Vector indicating which observations to plot. Default is `1:nrow(y)`, indicating that all obersvations
#' should be plotted.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"y"`.
#' @param title Title for the figure. Default is `""`.
#'
#' @return
#' A gglpot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_data <- function(t, y, indices = 1:nrow(y), xlab = "x", ylab = "y", title = "") {
  n <- ncol(y)
  N <- nrow(y)
  t1 <- min(t)
  tn <- max(t)

  if (is.matrix(t)){
    if (nrow(t) == 1 | ncol(t) == 1) {
      t <- as.numeric(t)
      if (length(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      data <- data.frame(x = rep(t, each = N), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
    } else {
      if (nrow(t) != N | ncol(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      data <- data.frame(x = c(t), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
    }
  } else {
    if (length(t) != n) {
      stop("The dimensions of 't' and 'y' don't match.")
    }
    data <- data.frame(x = rep(t, each = N), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
  }
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes(x = .data$x, y = .data$y, group = .data$Row_Index, color = .data$Row_Index)) +
       ggplot2::geom_line() +
       ggplot2::labs(x = xlab, y = ylab, title = title) +
       ggplot2::scale_color_viridis_d(option = "A", direction = 1) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none",
                      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                      axis.line.x = ggplot2::element_line(color = "black"),
                      axis.line.y = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                      panel.background = ggplot2::element_blank()) +
       ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

  return(p)
}

#' Plot a histogram and density estimate for a sample
#'
#' @description
#' `plot_posterior_hist()` returns the histogram and density estimate for a sample. In particular,
#' it allows to illustrate the posterior distribution of parameters or quantities of interest.
#'
#' @param sample Sample of values for which the distribution is of interest.
#' @param colors Vector of two colors. The first one is for the fill color of the histogram bins, the
#' second one is the line color for the density. Default is `c("gold", "darkblue")`.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"Density"`.
#' @param title Title for the figure. Default is `""`.
#' @param binwidth Width of histogram bins. Default is `0.03`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_posterior_hist <- function(sample, colors = c("gold", "darkblue"), xlab = "x", ylab = "Density",
                                title = "", binwidth = 0.03) {
  data <- data.frame("sample" = sample)

  temp_hist <- ggplot2::ggplot(data, ggplot2::aes(x = .data$sample)) +
               ggplot2::geom_histogram(binwidth = binwidth,
                                       ggplot2::aes(y = ggplot2::after_stat(.data$density)))

  hist_data <- ggplot2::ggplot_build(temp_hist)$data[[1]]
  max_density <- max(hist_data$density)

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$sample)) +
       ggplot2::geom_histogram(binwidth = binwidth,
                               ggplot2::aes(y = ggplot2::after_stat(.data$density)),
                               color = "black",
                               fill = colors[1]) +
       ggplot2::geom_density(lwd = 2, color = colors[2]) +
       ggplot2::labs(x = xlab, y = ylab, title = title) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none",
                      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                      axis.line.x = ggplot2::element_line(color = "black"),
                      axis.line.y = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                      panel.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(data), max(data))) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max_density + 0.2))

  return(p)
}

#' Plot a boxplot of individual mixed membership components
#'
#' @description
#' `plot_member_boxplot()` plots a boxplot for the membership to the provided feature, separated
#' by group. In the case study, the group corresponds to the clinical designation.
#'
#' @param pi Values of the mixed membership componenet for each individual.
#' @param group Vector of 0 and 1s indicating group belonging for each individual.
#' @param colors Vector of two colors, where each corresponds to the color of the boxplot for one group.
#' Default is `c("cornflowerblue", "gold")`.
#' @param group_labels Labels to use to denote each group. Default is `c("0", "1")`.
#' @param xlab Label for the horizontal axis. Default is `"Designation"`.
#' @param ylab Label for the vertical axis. Default is `"Value"`.
#' @param title Title for the figure. Default is `""`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
plot_member_boxplot <- function(pi, group, colors = c("cornflowerblue", "gold1"),
                            group_labels = c("0", "1"), xlab = "Designation",
                            ylab = "Value", title = "") {
  data_pi <- data.frame("pi" = pi,
                        "group" = group)
  data_0 <- data_pi[data_pi$group == 0, ]
  data_1 <- data_pi[data_pi$group == 1, ]

  p <- ggplot2::ggplot() +
       ggplot2::geom_boxplot(data = data_0,
                       mapping = ggplot2::aes(x = group_labels[1], y = .data$pi),
                       fill = ggplot2::alpha(colors[1], 0.5)) +
       ggplot2::geom_boxplot(data = data_1,
                        mapping = ggplot2::aes(x = group_labels[2], y = .data$pi),
                        fill = ggplot2::alpha(colors[2], 0.5)) +
       ggplot2::labs(x = xlab, y = ylab, title = title) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      legend.position = "none",
                      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                      axis.line.x = ggplot2::element_line(color = "black"),
                      axis.line.y = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                      panel.background = ggplot2::element_blank()) +
       ggplot2::geom_jitter(data = data_0, ggplot2::aes(x = group_labels[1], y = .data$pi),
                            color = "black", width = 0.1) +
       ggplot2::geom_jitter(data = data_1, ggplot2::aes(x = group_labels[2], y = .data$pi),
                            color = "black", width = 0.1)

  return(p)

}

#' Plot of multiple posterior feature shapes
#'
#' @description
#' `plot_multiple_feature()` plots the posterior shape function for the desired feature obtained through
#' different methods or simulations. The posterior quantiles of interest are also represented.
#'
#' @param eval_t Time points for the horizontal axis.
#' @param shape List whose elements are the posterior central summary of the shape function.
#' @param quantile_low List whose elements are the posterior lower quantile of the shape function.
#' @param quantile_high List whose elements are the posterior upper quantile of the shape function.
#' @param colors Vector of colors for the line and ribbon for each different shape.
#' @param alpha Opacity of the color for the quantile ribbon. Default is `0.2`.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"f"`.
#' @param title Title for the figure. Default is `""`.
#' @param background_ind Indices of observations to plot. Default is `NULL`, in which case no observations
#' are displayed with the posterior shape.
#' @param t Time points at which data is evaluated. Default is `NULL`, for the case when no observations
#' are displayed with the posterior shape.
#' @param y Matrix of data. Each row corresponds to one observation. Default is `NULL`, for the case when
#' no observations are displayed with the posterior shape.
#' @param legend_position Vector of coordinates for the location of the legend. Default is `c(0.8, 0.8)`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_multiple_feature <- function(eval_t, shape, quantile_low, quantile_high, colors,
                                  alpha = 0.2, xlab = "x", ylab = "f", title = "",
                                  background_ind = NULL, t = NULL, y = NULL,
                                  legend_position = c(0.8, 0.8)) {
  n <- length(eval_t)
  t1 <- eval_t[1]
  tn <- eval_t[n]

  list_names <- names(shape)
  m <- length(list_names)

  plot_dfs <- list()

  for (i in 1:m) {
    df <- data.frame(x = eval_t,
                     f = shape[i],
                     f_low = quantile_low[i],
                     f_high = quantile_high[i])
    df$warping <- list_names[i]
    plot_dfs[[paste0("df", i)]] <- df
  }

  plot_df <- do.call(rbind, plot_dfs)
  color_mapping <- stats::setNames(colors, list_names)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$f)) +
       ggplot2::geom_line(ggplot2::aes(color = .data$warping)) +
       ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$f_low,
                                         ymax = .data$f_high,
                                         fill = .data$warping),
                            alpha = 0.2) +
       ggplot2::labs(x = xlab, y = ylab, title = title) +
       ggplot2::scale_color_manual(values = color_mapping) +
       ggplot2::scale_fill_manual(values = color_mapping) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust =0.5),
                      legend.position = legend_position,
                      legend.title = ggplot2::element_blank(),
                      legend.spacing.x = grid::unit(0, "cm"),
                      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                      axis.line.x = ggplot2::element_line(color = "black"),
                      axis.line.y = ggplot2::element_line(color = "black"),
                      axis.ticks.y = ggplot2::element_line(color = "black"),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                      panel.background = ggplot2::element_blank()) +
       ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

  return(p)
}


