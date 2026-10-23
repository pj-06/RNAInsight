############################################################
# RNAInsight Dashboard
# 01_theme.R : colour tokens, ggplot2 theme, plotly styling
############################################################
#
# One dark visual system is shared by every chart in the app
# (ggplot2, plotly and the 3D/2D structure viewer), so
# the dashboard reads as one product. Colours are defined once
# here and referred to by role everywhere else.
############################################################

#' Colour tokens used across the dashboard.
#'
#' Series colours follow a validated colour-blind-safe order
#' (blue, orange, aqua, yellow, magenta) stepped for a dark
#' surface. RNA families (11 classes) are never given 11 hues:
#' family identity is carried by axis position or labels, and
#' colour is used to highlight one family against the rest.
RI_COLORS <- list(
  page        = "#0b0f0e",
  surface     = "#121816",
  card        = "#161d1b",
  card_alt    = "#1b2421",
  border      = "#26302c",
  grid        = "#232c29",
  axis        = "#3a4540",
  text        = "#eef1ef",
  text_2      = "#b9c2bd",
  muted       = "#84908a",
  accent      = "#3ecf8e",
  accent_dim  = "#1f6b4c",
  neutral     = "#4a5752",
  series      = c("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181"),
  seq_low     = "#16222b",
  seq_mid     = "#2a78d6",
  seq_high    = "#cde2fb",
  div_neg     = "#3987e5",
  div_mid     = "#383835",
  div_pos     = "#e66767"
)

#' Five models in a fixed colour order (colour follows the model,
#' never its rank, so filtering never repaints a model).
MODEL_ORDER <- c(
  "Decision Tree (C5.0)",
  "Random Forest",
  "Gradient Boosting Machine",
  "K-Nearest Neighbors",
  "Multinomial Logistic Regression"
)

MODEL_SHORT <- c(
  "Decision Tree (C5.0)"            = "C5.0",
  "Random Forest"                   = "Random Forest",
  "Gradient Boosting Machine"       = "GBM",
  "K-Nearest Neighbors"             = "KNN",
  "Multinomial Logistic Regression" = "Multinomial LR"
)

MODEL_COLORS <- stats::setNames(RI_COLORS$series, MODEL_ORDER)

#' ggplot2 theme matching the dashboard surface.
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
theme_rnainsight <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = RI_COLORS$card, colour = NA),
      panel.background  = ggplot2::element_rect(fill = RI_COLORS$card, colour = NA),
      panel.grid.major  = ggplot2::element_line(colour = RI_COLORS$grid, linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(colour = RI_COLORS$axis, linewidth = 0.3),
      axis.text         = ggplot2::element_text(colour = RI_COLORS$muted),
      axis.title        = ggplot2::element_text(colour = RI_COLORS$text_2),
      plot.title        = ggplot2::element_text(colour = RI_COLORS$text, face = "bold"),
      plot.subtitle     = ggplot2::element_text(colour = RI_COLORS$text_2),
      plot.caption      = ggplot2::element_text(colour = RI_COLORS$muted),
      legend.background = ggplot2::element_rect(fill = RI_COLORS$card, colour = NA),
      legend.key        = ggplot2::element_rect(fill = RI_COLORS$card, colour = NA),
      legend.text       = ggplot2::element_text(colour = RI_COLORS$text_2),
      legend.title      = ggplot2::element_text(colour = RI_COLORS$text_2),
      strip.text        = ggplot2::element_text(colour = RI_COLORS$text, face = "bold"),
      strip.background  = ggplot2::element_rect(fill = RI_COLORS$card_alt, colour = NA)
    )
}

#' Apply the shared dark styling to a plotly object.
#'
#' Works for native plot_ly() charts and for ggplotly()
#' conversions, so both kinds look identical in the app.
#'
#' @param p A plotly object.
#' @param show_legend Show the legend?
#' @param ... Extra arguments passed to plotly::layout().
#' @return The styled plotly object.
ri_plotly <- function(p, show_legend = TRUE, ...) {
  axis_style <- list(
    gridcolor = RI_COLORS$grid,
    zerolinecolor = RI_COLORS$axis,
    linecolor = RI_COLORS$axis,
    tickfont = list(color = RI_COLORS$muted),
    title = list(font = list(color = RI_COLORS$text_2))
  )
  extra <- list(...)
  # Merge caller axis settings on top of the shared style (a plain
  # string title becomes list(text = ...) so the font is kept).
  merge_axis <- function(user) {
    if (is.null(user)) return(axis_style)
    if (!is.null(user$title) && is.character(user$title)) {
      user$title <- list(text = user$title, font = axis_style$title$font)
    }
    utils::modifyList(axis_style, user)
  }
  base <- list(
    paper_bgcolor = "rgba(0,0,0,0)",
    plot_bgcolor = "rgba(0,0,0,0)",
    font = list(
      family = "Inter, system-ui, -apple-system, 'Segoe UI', sans-serif",
      color = RI_COLORS$text_2, size = 12
    ),
    showlegend = show_legend,
    legend = list(bgcolor = "rgba(0,0,0,0)", font = list(color = RI_COLORS$text_2)),
    hoverlabel = list(
      bgcolor = RI_COLORS$card_alt,
      bordercolor = RI_COLORS$border,
      font = list(color = RI_COLORS$text)
    ),
    margin = list(l = 10, r = 10, t = 30, b = 10)
  )
  base$xaxis <- merge_axis(extra$xaxis)
  base$yaxis <- merge_axis(extra$yaxis)
  extra$xaxis <- NULL
  extra$yaxis <- NULL
  # second panel of two-panel subplots
  base$xaxis2 <- utils::modifyList(axis_style, list(title = NULL))
  base$yaxis2 <- utils::modifyList(axis_style, list(title = NULL))
  if (!is.null(extra$xaxis2)) base$xaxis2 <- utils::modifyList(base$xaxis2, extra$xaxis2)
  extra$xaxis2 <- NULL
  if (!is.null(extra$legend)) {
    base$legend <- utils::modifyList(base$legend, extra$legend)
    extra$legend <- NULL
  }
  args <- utils::modifyList(base, extra)
  do.call(plotly::layout, c(list(p), args)) |>
    plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
}

#' Convert a ggplot to plotly with the shared styling.
#'
#' @param g A ggplot object.
#' @param tooltip Aesthetics to show in the tooltip.
#' @param show_legend Show the legend?
#' @return A styled plotly object.
ri_ggplotly <- function(g, tooltip = "text", show_legend = TRUE) {
  ri_plotly(plotly::ggplotly(g, tooltip = tooltip), show_legend = show_legend)
}

#' bslib theme for the whole dashboard.
#'
#' @return A bslib::bs_theme object.
ri_bs_theme <- function() {
  bslib::bs_theme(
    version = 5,
    bg = RI_COLORS$page,
    fg = RI_COLORS$text,
    primary = RI_COLORS$accent,
    secondary = RI_COLORS$neutral,
    base_font = bslib::font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "sans-serif"),
    heading_font = bslib::font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "sans-serif"),
    code_font = bslib::font_collection("JetBrains Mono", "Consolas", "monospace"),
    "border-color" = RI_COLORS$border,
    "card-bg" = RI_COLORS$card,
    "card-border-color" = RI_COLORS$border
  )
}
