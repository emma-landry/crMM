# Case-study-specific helper functions used in the EEG application.
# These functions rely on model components or assumptions that do not
# generalize to arbitrary K-feature mixed membership models.

#' Posterior estimates for the time-transformation regression parameter
#'
#' @description
#' `posterior_regression()` obtains the posterior summary statistics for the time-transformation regression
#' parameter. It calculates the desired summary statistics / moments of the posterior MCMC samples for the
#' regression parameter. If the sample provided considers the analysis of the EEG case study, the
#' mean time-transformation functions associated with inputted ages, and stratified by clinical designation,
#' are also fitted.
#'
#'
#' @param crMM_samples An object of class `crMM_Obj`.
#' @param moment Posterior moment of interest. Options are: `"mean"`, `"median"`, `"sd"` or values
#' between 0 and 1 for quantiles. The default is `"mean"`.
#' @param h The number of inner knots used to define the B-spline for the time-transformation functions.
#' @param degree Degree of the piecewise cubic polynomial of B-splines for the time-transformation
#' functions. The default is `3` for cubic splines.
#' @param intercept If `TRUE`, an intercept is included in the B-spline basis for the time-transformation
#' functions. The default is `FALSE`.
#' @param t_boundary A vector with the first and last observation time point. Default is `c(0,1)`.
#' @param eval_points When `CaseStudy = TRUE`, number of points at which to evaluate the time-transformation
#' functions using B-spines. Default is `100`.
#' @param CaseStudy Logical that indicates whether the sample consists of the analysis of the EEG case study.
#' The default is `FALSE`.
#' @param ages When `CaseStudy = TRUE`, the age values at which to evaluate the regression mean for the
#' time-transformation functions. Default is `NULL`, for the general case.
#'
#' @return
#' When `CaseStudy = FALSE`, returns the posterior sample summary for the regression coefficient.
#'
#' When `CaseStudy = TRUE`, returns a list containing:
#'
#' * `B`: the posterior sample summary of the regression coefficient
#' * `means0`: a matrix with the regression mean function evaluated at different ages, for TD individuals
#' * `means1`: a matrix with the regression mean function evaluated at different ages, for ASD individuals
#'
#' @export
#'
posterior_regression <- function(crMM_samples,  moment = "mean", h, degree = 3, intercept = FALSE,
                                 CaseStudy = FALSE, ages = NULL, t_boundary = c(0, 1), eval_points = 100) {

  if (!inherits(crMM_samples, "crMM_Obj")) {
    stop("'crMM_samples' must be an object of class 'crMM_Obj'.")
  }

  if (intercept == FALSE) {
    df <- h + degree
  } else {
    df <- h + degree + 1
  }

  if (CaseStudy == FALSE) {
    l <- ncol(crMM_samples$B) / (df - 2)

    if (round(l) != l) {
      stop("The provided B-spline parameters do not match the dimensions implied by the regression
           coefficient.")
    }

    if (typeof(moment) == "character") {
      if (moment == "mean") {
        B <- apply(crMM_samples$B, 2, mean)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else if (moment == "median") {
        B <- apply(crMM_samples$B, 2, stats::median)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else if (moment == "sd") {
        B <- apply(crMM_samples$B, 2, stats::sd)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    } else {
      if (moment >= 0 & moment <= 1) {
        B <- apply(crMM_samples$B, 2, stats::quantile, probs = moment)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    }
    return(B)
  } else {
    l <- 3

    if (ncol(crMM_samples$B) != l * (df - 2)) {
      stop("The provided B-spline parameters do not match the dimensions implied by the regression
           coefficient.")
    }

    if (is.null(ages)) {
      stop("'ages' needs to be provided for case study analysis.")
    }

    if (typeof(moment) == "character") {
      if (moment == "mean") {
        B <- apply(crMM_samples$B, 2, mean)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else if (moment == "median") {
        B <- apply(crMM_samples$B, 2, stats::median)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else if (moment == "sd") {
        B <- apply(crMM_samples$B, 2, stats::sd)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    } else {
      if (moment >= 0 & moment <= 1) {
        B <- apply(crMM_samples$B, 2, stats::quantile, probs = moment)
        B <- matrix(B, nrow = l, ncol = df - 2, byrow = F)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    }
    n_ages <- length(ages)
    X0 <- cbind(ages, rep(0, n_ages), rep(0, n_ages))
    X1 <- cbind(ages, rep(1, n_ages), ages)

    t1 <- min(t_boundary)
    tn <- max(t_boundary)
    knots <- seq(t1, tn, length.out = h + 2)[2:(h + 1)]

    phi_mean <- identityTT(Boundary.knots = t_boundary, knots = knots,
                           degree = degree, intercept = intercept)
    jupp_mean <- jupp(phi_mean)[-c(1, df)]
    jupp_mean_mat <- matrix(rep(jupp_mean, n_ages), nrow = n_ages, byrow = T)

    reg_mean0 <- apply(cbind(cbind(rep(0, n_ages), X0 %*% B  + jupp_mean_mat), rep(1, n_ages)), 1, juppinv)
    reg_mean1 <- apply(cbind(cbind(rep(0, n_ages), X1 %*% B  + jupp_mean_mat), rep(1, n_ages)), 1, juppinv)

    eval_t <- seq(t1, tn, length.out = eval_points)
    basis <- splines::bs(x = eval_t, knots = knots, degree = degree, intercept = intercept)

    means0 <- basis %*% reg_mean0
    means1 <- basis %*% reg_mean1

    return(list(B = B,
                means0 = means0,
                means1 = means1))
  }
}

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
  t1 <- min(t)
  tn <- max(t)

  num_indices <- length(indices)

  if (length(line_colors) != num_indices) {
    stop("Provide line colors of the same length as the number of fitted lines to be plotted.")
  }

  if (is.matrix(t)){
    if (nrow(t) == 1 | ncol(t) == 1) {
      t <- as.numeric(t)
      if (length(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      x <- rep(t, times = N)
      x_ind <- rep(t, times = num_indices)
    } else {
      if (nrow(t) != N | ncol(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      x <- c(t(t))
      x_ind <- matrix(t[indices, ], nrow = n * num_indices, byrow = TRUE)
      x_ind <- as.numeric(x_ind)
    }
  } else {
    if (length(t) != n) {
      stop("The dimensions of 't' and 'y' don't match.")
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
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none",
                   axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                   axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                   axis.line.x = ggplot2::element_line(color = "black"),
                   axis.line.y = ggplot2::element_line(color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                   panel.background = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(10, 10, 10, 10)) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(t1, tn))

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
                   axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
                   axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0)),
                   axis.line.x = ggplot2::element_line(color = "black"),
                   axis.line.y = ggplot2::element_line(color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
                   panel.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

  if (!is.null(background_ind)) {
    for (i in 1:length(background_ind)){
      line_data <- data.frame(x = t, y = y[background_ind[i], ])
      p <- p + ggplot2::geom_line(data = line_data,
                                  ggplot2::aes(x = .data$x, y = .data$y),
                                  color = "gray68",
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
#' @details
#' When `group` is provided, functions are split by group membership (e.g., TD vs ASD
#' in the EEG case study). When `covariate` is also provided, lines are colored by the
#' covariate value within each group.
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
#' @param viridis_option String indicating color map to use for `scale_color_viridis_d()`. Default is `"C"`.
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
                    viridis_option = "C",
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
      ggplot2::geom_line(linewidth = 0.4) +
      ggplot2::labs(x = xlab, y = ylab, title = title) +
      ggplot2::scale_color_viridis_d(option = viridis_option, direction = 1) +
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
        ggplot2::scale_color_viridis_c(option = viridis_option, direction = -1) +
        #ggplot2::scale_color_gradient(low = covariate_gradient[1], high = covariate_gradient[2]) +
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
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(tt_functions), max(tt_functions)))

      p1 <- ggplot2::ggplot(data = df_1, ggplot2::aes(x = .data$x,
                                                      y = .data$y,
                                                      group = .data$group,
                                                      color = .data$cov)) +
        ggplot2::geom_line() +
        ggplot2::labs(x = xlab, y = ylab, color = covlab) +
        ggplot2::scale_color_viridis_c(option = viridis_option, direction = -1) +
        #ggplot2::scale_color_gradient(low = covariate_gradient[1], high = covariate_gradient[2]) +
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
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(tt_functions), max(tt_functions)))

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
#' @param viridis_option String indicating color map to use for `scale_color_viridis_d()`. Default is `"C"`.
#' If `NULL`, plots in a single color rather than with the viridis gradient.
#' @param color Color to plot in if `viridis_option` is `NULL`. Default is `"grey"`.

#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_data <- function(t, y, indices = 1:nrow(y), xlab = "x", ylab = "y", title = "", viridis_option = "C",
                      color = "grey") {
  n <- ncol(y)
  N <- nrow(y)
  y <- y[indices, ]
  t1 <- min(t)
  tn <- max(t)

  if (is.matrix(t)){
    if (nrow(t) == 1 | ncol(t) == 1) {
      N <- nrow(y)
      t <- as.numeric(t)
      if (length(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      data <- data.frame(x = rep(t, each = N), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
    } else {
      if (nrow(t) != N | ncol(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      t <- t[indices, ]
      N <- nrow(y)
      data <- data.frame(x = c(t), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
    }
  } else {
    N <- nrow(y)
    if (length(t) != n) {
      stop("The dimensions of 't' and 'y' don't match.")
    }
    data <- data.frame(x = rep(t, each = N), y = c(y), Row_Index = paste0("X", rep(1:N, n)))
  }

  if (!is.null(viridis_option)) {
    p <- ggplot2::ggplot(data = data,
                         ggplot2::aes(x = .data$x, y = .data$y, group = .data$Row_Index,
                                      color = .data$Row_Index)) +
      ggplot2::geom_line() +
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
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

    p <- p + ggplot2::scale_color_viridis_d(option = viridis_option, direction = 1)
  } else {
    p <- ggplot2::ggplot(data = data,
                         ggplot2::aes(x = .data$x, y = .data$y, group = .data$Row_Index,
                                      color = .data$Row_Index, linetype = .data$Row_Index)) +
      ggplot2::geom_line() +
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
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn))

    available_linetypes <- c("solid", "dashed", "dotted", "dotdash",
                             "longdash", "twodash", "1F", "2F", "3F",
                             "4F", "12345678", "longdash", "twodash")
    unique_rows <- unique(data$Row_Index)
    set.seed(15)
    random_linetypes <- sample(available_linetypes, length(unique_rows), replace = TRUE)
    p <- p + ggplot2::scale_color_manual(values = rep(color, length(unique(data$Row_Index)))) +
      ggplot2::scale_linetype_manual(values = random_linetypes)
  }

  return(p)
}

#' Plot of the data, split by group and colored by covariate value
#'
#' @description
#' `plot_data_cd()` plots the rows of the matrix of observations, in separate panels by clinical designation.
#' The lines are colored by age value for each individual.
#'
#' @details
#' This function is intended for the EEG case study, where `group` corresponds
#' to clinical designation (TD vs ASD).
#'
#' @param t  Time points for the x axis. Either a vector that is common for all observations, or a matrix for which
#' each row corresponds to an observation.
#' @param y Matrix of data. Each row corresponds to one observation.
#' @param group Vector of 0 and 1s indicating group belonging for each function. Default is `NULL`, for when the group
#' is not provided and all time transformation functions are plotted together.
#' @param covariate Vector of covariate values associated with each function. Default is `NULL`, for when covariate
#' value is not provided.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"y"`.
#' @param title Title for the figure. Default is `""`.
#' @param covlab Title for the legend for covariate values. Default is `"Age"`.
#' @param group_labels Labels to use to denote each group. Default is `c("0", "1")`.
#' @param viridis_option String indicating color map to use for `scale_color_viridis_d()`. Default is `"C"`.
#' If `NULL`, plots in a single color rather than with the viridis gradient.
#' @param legend_position Vector of coordinates for the location of the legend. Default is `c(0.85, 0.78)`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_data_cd <- function(t, y, group, covariate,
                         xlab = "x", ylab = "y", title = "", covlab = "Age",
                         group_labels = c("0", "1"), viridis_option = "C",
                         legend_position = c(0.85, 0.78)){

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
      plot_df <- data.frame(x = rep(t, N),
                            y = c(t(y)),
                            group = rep(1:N, each = n),
                            vec = rep(group, each = n),
                            cov = rep(covariate, each = n))
    } else {
      if (nrow(t) != N | ncol(t) != n) {
        stop("The dimensions of 't' and 'y' don't match.")
      }
      plot_df <- data.frame(x = c(t(t)),
                            y = c(t(y)),
                            group = rep(1:N, each = n),
                            vec = rep(group, each = n),
                            cov = rep(covariate, each = n))
    }
  } else {
    if (length(t) != n) {
      stop("The dimensions of 't' and 'y' don't match.")
    }
    plot_df <- data.frame(x = rep(t, N),
                          y = c(t(y)),
                          group = rep(1:N, each = n),
                          vec = rep(group, each = n),
                          cov = rep(covariate, each = n))
  }

  df_0 <- plot_df[plot_df$vec == 0, ]
  df_1 <- plot_df[plot_df$vec == 1, ]

  p0 <- ggplot2::ggplot(data = df_0, ggplot2::aes(x = .data$x,
                                                  y = .data$y,
                                                  group = .data$group,
                                                  color = .data$cov)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = xlab, y = ylab, color = covlab) +
    ggplot2::scale_color_viridis_c(option = viridis_option, direction = -1) +
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
                   plot.margin = ggplot2::margin(10, 10, 10, 10),
                   legend.position = legend_position) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(y) - 0.01, max(y) + 0.01))

  p1 <- ggplot2::ggplot(data = df_1, ggplot2::aes(x = .data$x,
                                                  y = .data$y,
                                                  group = .data$group,
                                                  color = .data$cov)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = xlab, y = ylab, color = covlab) +
    ggplot2::scale_color_viridis_c(option = viridis_option, direction = -1) +
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
                   plot.margin = ggplot2::margin(10, 10, 10, 10),
                   legend.position = legend_position) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(t1, tn)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(y) - 0.01, max(y) + 0.01))

  p <- gridExtra::grid.arrange(p0, p1, ncol = 2,
                               top = grid::textGrob(title,
                                                    gp = grid::gpar(fontsize = 15, font = 2)))
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
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max_density + 0.3))

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
                          fill = ggplot2::alpha(colors[1], 1)) +
    ggplot2::geom_boxplot(data = data_1,
                          mapping = ggplot2::aes(x = group_labels[2], y = .data$pi),
                          fill = ggplot2::alpha(colors[2], 1)) +
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
                     f = shape[[i]],
                     f_low = quantile_low[[i]],
                     f_high = quantile_high[[i]])
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
                         alpha = alpha) +
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

  if (!is.null(background_ind)) {
    for (i in 1:length(background_ind)){
      line_data <- data.frame(x = t, y = y[background_ind[i], ])
      p <- p + ggplot2::geom_line(data = line_data,
                                  ggplot2::aes(x = .data$x, y = .data$y),
                                  color = "gray68",
                                  linetype = "dashed")
    }
  }

  return(p)
}

#' Plot shape functions from different samples
#'
#' @description
#' `plot_identifiability()` plots the posterior shape function for the desired feature obtained through
#' different methods or simulations. It focuses on the identifiability, for plots with confidence bands
#' and legends see `plot_multiple_feature()`.
#'
#' @param eval_t Time points for the horizontal axis.
#' @param shape List whose elements are the values of the shape function for different simulations.
#' @param colors Vectors of colors for the line of each shape. Default is
#' `("red", "blue", "green", "orange", "purple")`, assuming 5 lines are plotted.
#' @param xlab Label for the horizontal axis. Default is `"x"`.
#' @param ylab Label for the vertical axis. Default is `"f"`.
#' @param title Title for the figure. Default is `""`.
#'
#' @return
#' A ggplot figure
#'
#' @export
#'
#' @importFrom rlang .data
#'
plot_identifiability <- function(eval_t, shape,
                                 colors = c("red", "blue", "green", "orange", "purple"),
                                 xlab = "x", ylab = "f", title = "") {
  n <- length(eval_t)
  t1 <- eval_t[1]
  tn <- eval_t[n]

  list_names <- names(shape)
  m <- length(list_names)

  if (length(colors) != m) {
    stop("Please provide as many colors as there are samples to represent.")
  }

  plot_df <- data.frame(x = eval_t)

  for (i in 1:m) {
    plot_df[[paste0("sample", i)]] <- shape[[i]]
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x))

  for (i in 1:m) {
    col_name <- paste0("sample", i)
    p <- p + ggplot2::geom_line(ggplot2::aes(y = !!rlang::sym(col_name)), color = colors[i])
  }

  p <- p + ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust =0.5),
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


