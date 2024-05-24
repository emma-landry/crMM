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
      x <- matrix(t, nrow = n * N, byrow = TRUE)
      x <- as.numeric(x)
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
                              y = y,
                              group = "Background")

  group <- factor(rep(indices, each = n))

  index_df <- data.frame(x = x_ind,
                         y = as.vector(t(y[indices,])),
                         group = group)

  lines_df <- data.frame(x = x_ind,
                         y = as.vector(t(fit[indices, ])),
                         group = group)

  p <- ggplot2::ggplot() +
       ggplot2::geom_point(data = background_df, ggplot2::aes(x = x, y = y), shape = 4, size = 2, color = "grey80") +
       ggplot2::geom_point(data = index_df, ggplot2::aes(x = x, y = y, color = group), shape = 1, size = 2) +
       ggplot2::geom_line(data = lines_df, ggplot2::aes(x = x, y = y, group = group, color = group), size = 0.7) +
       ggplot2::scale_color_manual(values = line_colors) +
       ggplot2::labs(x = xlab, y = ylab, title = title ) +
       ggplot2::theme_classic() +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")

  return(p)
}


