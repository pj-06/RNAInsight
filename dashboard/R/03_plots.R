############################################################
# RNAInsight Dashboard
# 03_plots.R : chart builders for the Results section
############################################################
#
# Every function takes the pre-computed dashboard data (see
# 02_data.R) plus the user's choices and returns a plotly or
# ggplot object. Keeping them free of Shiny code means each
# chart can also be rebuilt from the console for the report:
#
#   d <- load_dashboard_data(".", "dashboard")
#   plot_confusion(d, "Random Forest")
#
# Colour use (see 01_theme.R):
#   * RNA families are never given 11 hues - one family can be
#     highlighted (blue) against the rest (grey).
#   * Models keep a fixed colour each.
#   * Magnitudes use one-hue sequential scales, signed values a
#     blue-grey-red diverging scale.
############################################################

HIGHLIGHT_NONE <- "(none)"

#' Colour scale for "one family vs the rest" highlighting.
#'
#' @param highlight Family name or HIGHLIGHT_NONE.
#' @return Named vector for scale_*_manual.
highlight_scale <- function(highlight) {
  if (identical(highlight, HIGHLIGHT_NONE)) {
    c("All records" = RI_COLORS$series[1])
  } else {
    stats::setNames(c(RI_COLORS$series[1], RI_COLORS$neutral), c(highlight, "Other families"))
  }
}

#' Group label for highlighting a single family.
#'
#' @param family Factor/character of families.
#' @param highlight Family name or HIGHLIGHT_NONE.
#' @return Factor with the highlight group of each record.
highlight_group <- function(family, highlight) {
  if (identical(highlight, HIGHLIGHT_NONE)) {
    return(factor(rep("All records", length(family))))
  }
  factor(ifelse(family == highlight, highlight, "Other families"), levels = c(highlight, "Other families"))
}

#' Sequential plotly colourscale (dark surface -> blue -> pale blue).
seq_colorscale <- function() {
  list(c(0, RI_COLORS$seq_low), c(0.55, RI_COLORS$seq_mid), c(1, RI_COLORS$seq_high))
}

#' Diverging plotly colourscale (blue - grey - red).
div_colorscale <- function() {
  list(c(0, RI_COLORS$div_neg), c(0.5, RI_COLORS$div_mid), c(1, RI_COLORS$div_pos))
}

#' Format a proportion as a percentage string.
pct <- function(x, digits = 1) {
  ifelse(is.na(x), "n/a", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

# ==========================================================
# Overview
# ==========================================================

#' Records vs unique ligands per RNA family (grouped bars).
#'
#' @param d Dashboard data.
#' @return plotly object.
plot_family_composition <- function(d) {
  records <- as.data.frame(d$family_counts, stringsAsFactors = FALSE)
  names(records) <- c("Family", "Records")
  records$Ligands <- vapply(records$Family, function(f) {
    length(unique(d$cleaned$SMILES[d$cleaned$Target_Family == f]))
  }, numeric(1))
  records$Family <- factor(records$Family, levels = rev(records$Family))

  plotly::plot_ly(records, y = ~Family) |>
    plotly::add_bars(
      x = ~Records, name = "Records (rows)", marker = list(color = RI_COLORS$series[1]),
      text = ~Records, textposition = "outside", textfont = list(color = RI_COLORS$text_2),
      hovertemplate = "%{y}<br>%{x} records<extra></extra>"
    ) |>
    plotly::add_bars(
      x = ~Ligands, name = "Unique ligands", marker = list(color = RI_COLORS$series[3]),
      text = ~Ligands, textposition = "outside", textfont = list(color = RI_COLORS$text_2),
      hovertemplate = "%{y}<br>%{x} unique ligands<extra></extra>"
    ) |>
    ri_plotly(barmode = "group", bargap = 0.25,
              xaxis = list(title = "Count", gridcolor = RI_COLORS$grid),
              yaxis = list(title = ""),
              legend = list(orientation = "h", x = 0, y = 1.08, font = list(color = RI_COLORS$text_2)))
}

#' Heatmap: how many ligands two RNA families share.
#'
#' @param d Dashboard data.
#' @return plotly object.
plot_family_overlap <- function(d) {
  ov <- d$overlap
  fams <- d$families
  m <- matrix(0, length(fams), length(fams), dimnames = list(fams, fams))
  for (i in seq_len(nrow(ov))) m[ov$Family_A[i], ov$Family_B[i]] <- ov$Shared[i]
  text <- matrix(as.character(m), nrow(m))
  text[m == 0] <- ""

  plotly::plot_ly(
    x = fams, y = fams, z = m, type = "heatmap",
    colorscale = seq_colorscale(), xgap = 2, ygap = 2,
    text = text, texttemplate = "%{text}",
    hovertemplate = "%{y} & %{x}<br>%{z} shared ligands<extra></extra>",
    colorbar = list(title = list(text = "Ligands", font = list(color = RI_COLORS$text_2)),
                    tickfont = list(color = RI_COLORS$muted), outlinewidth = 0)
  ) |>
    ri_plotly(show_legend = FALSE,
              xaxis = list(title = "", tickangle = -40),
              yaxis = list(title = "", autorange = "reversed"))
}

# ==========================================================
# Data exploration
# ==========================================================

#' Histogram of one descriptor with an optional highlighted family.
#'
#' @param d Dashboard data.
#' @param feature Column name.
#' @param highlight Family or HIGHLIGHT_NONE.
#' @return plotly object (built with ggplot2).
plot_distribution <- function(d, feature, highlight) {
  df <- d$cleaned
  df$Group <- highlight_group(df$Target_Family, highlight)
  df$Value <- df[[feature]]

  # Grey "rest" layer first, highlighted family drawn on top of it.
  breaks <- pretty(range(df$Value), n = 30)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = Value, fill = Group)) +
    ggplot2::geom_histogram(data = df[df$Group != highlight, ], breaks = breaks,
                            colour = RI_COLORS$card, linewidth = 0.4, alpha = 0.95) +
    ggplot2::geom_histogram(data = df[df$Group == highlight, ], breaks = breaks,
                            colour = RI_COLORS$card, linewidth = 0.4, alpha = 0.95) +
    ggplot2::scale_fill_manual(values = highlight_scale(highlight), name = NULL) +
    ggplot2::labs(x = feature_label(feature), y = "Records") +
    theme_rnainsight()

  ri_ggplotly(g, tooltip = c("x", "y", "fill")) |>
    plotly::layout(legend = list(orientation = "h", x = 0, y = 1.1))
}

#' Interactive scatter plot of two descriptors.
#'
#' @param d Dashboard data.
#' @param x,y Column names.
#' @param highlight Family or HIGHLIGHT_NONE.
#' @param size_by Column name for marker size, or "" for fixed size.
#' @return plotly object.
plot_scatter <- function(d, x, y, highlight, size_by = "") {
  df <- d$ligands
  fam_flag <- vapply(strsplit(df$Families, ", "), function(f) highlight %in% f, logical(1))
  df$Group <- if (identical(highlight, HIGHLIGHT_NONE)) "All ligands" else
    ifelse(fam_flag, highlight, "Other families")
  df$tip <- paste0("<b>", df$Name, "</b><br>", df$Families,
                   "<br>", x, ": ", signif(df[[x]], 4), "<br>", y, ": ", signif(df[[y]], 4))
  df <- df[order(df$Group != "Other families"), ]

  colors <- if (identical(highlight, HIGHLIGHT_NONE)) c("All ligands" = RI_COLORS$series[1]) else
    stats::setNames(c(RI_COLORS$series[1], RI_COLORS$neutral), c(highlight, "Other families"))
  size <- if (nzchar(size_by)) df[[size_by]] else rep(1, nrow(df))

  plotly::plot_ly(
    df, x = df[[x]], y = df[[y]], type = "scatter", mode = "markers",
    color = ~Group, colors = colors, size = size, sizes = c(40, 400),
    marker = list(opacity = 0.9, sizemode = "area"),
    text = ~tip, hoverinfo = "text", key = ~key, source = "scatter"
  ) |>
    ri_plotly(xaxis = list(title = feature_label(x)), yaxis = list(title = feature_label(y)),
              legend = list(orientation = "h", x = 0, y = 1.08))
}

# ==========================================================
# Feature engineering
# ==========================================================

#' Why each of the 20 candidate features was kept or dropped.
#'
#' Bar length = the feature's strongest absolute correlation with
#' any other candidate; the dashed line is the 0.90 cut-off used
#' by caret::findCorrelation in 04_feature_engineering.R.
#'
#' @param d Dashboard data.
#' @return plotly object.
plot_feature_selection <- function(d) {
  sel <- d$selection
  m <- abs(stats::cor(d$cleaned[, sel$Feature]))
  diag(m) <- 0
  sel$MaxCor <- apply(m, 1, max)
  sel$Partner <- colnames(m)[apply(m, 1, which.max)]
  sel$Status <- ifelse(sel$Kept, "Kept", ifelse(sel$Dropped_NZV, "Dropped: near-zero variance", "Dropped: |r| > 0.90"))
  sel$Type <- ifelse(sel$Feature %in% ENGINEERED_FEATURES, "engineered", "RDKit")
  sel <- sel[order(sel$MaxCor), ]
  sel$Feature <- factor(sel$Feature, levels = sel$Feature)
  colors <- c("Kept" = RI_COLORS$series[3], "Dropped: |r| > 0.90" = RI_COLORS$series[2],
              "Dropped: near-zero variance" = RI_COLORS$neutral)

  plotly::plot_ly(sel, y = ~Feature, x = ~MaxCor, type = "bar", orientation = "h",
                  color = ~Status, colors = colors,
                  text = ~paste0(Feature, " (", Type, ")<br>max |r| = ", round(MaxCor, 3),
                                 " with ", Partner, "<br>", Status),
                  hoverinfo = "text", textposition = "none") |>
    plotly::layout(shapes = list(list(type = "line", x0 = 0.9, x1 = 0.9, y0 = -0.5, y1 = nrow(sel) - 0.5,
                                      line = list(color = RI_COLORS$text_2, dash = "dash", width = 1)))) |>
    ri_plotly(xaxis = list(title = "Strongest |correlation| with another feature", range = c(0, 1.02)),
              yaxis = list(title = ""),
              legend = list(orientation = "h", x = 0, y = 1.06))
}

#' Model-based importance from the tree / ensemble models.
#'
#' @param d Dashboard data.
#' @return plotly object.
plot_model_importance <- function(d) {
  imp <- d$importance
  if (is.null(imp) || !nrow(imp)) {
    return(empty_plot("Model-based importance needs the randomForest / gbm / C50 packages."))
  }
  order <- stats::aggregate(Importance ~ Feature, imp, mean)
  imp$Feature <- factor(imp$Feature, levels = order$Feature[order(order$Importance)])
  imp$Model <- factor(imp$Model, levels = intersect(MODEL_ORDER, unique(imp$Model)))

  plotly::plot_ly(imp, y = ~Feature, x = ~Importance, color = ~Model,
                  colors = MODEL_COLORS[levels(imp$Model)], type = "bar", orientation = "h",
                  hovertemplate = "%{y}<br>%{x:.1f}<extra>%{fullData.name}</extra>") |>
    ri_plotly(barmode = "group", xaxis = list(title = "Scaled importance (0-100)"), yaxis = list(title = ""),
              legend = list(orientation = "h", x = 0, y = 1.06))
}

# ==========================================================
# Model comparison
# ==========================================================


#' Test-set accuracy next to cross-validated accuracy for every model.
#'
#' @param d Dashboard data.
#' @return plotly object.
plot_accuracy_comparison <- function(d) {
  lb <- d$leaderboard
  lb <- lb[order(-lb$Accuracy), ]
  lb$Short <- factor(MODEL_SHORT[lb$Model], levels = MODEL_SHORT[lb$Model])
  plotly::plot_ly(lb, x = ~Short) |>
    plotly::add_bars(y = ~CV_Accuracy, name = "Cross-validated accuracy (mean of 15 folds)",
                     marker = list(color = RI_COLORS$series[1]),
                     error_y = list(type = "data", array = ~CV_Accuracy_SD, color = RI_COLORS$text_2, thickness = 1),
                     text = ~sprintf("%.3f", CV_Accuracy), textposition = "inside",
                     insidetextanchor = "end", textfont = list(color = "#ffffff"),
                     hovertemplate = "%{x}<br>CV accuracy %{y:.3f}<extra></extra>") |>
    plotly::add_bars(y = ~Accuracy, name = "Test accuracy (145 held-out records)",
                     marker = list(color = RI_COLORS$series[2]),
                     text = ~sprintf("%.3f", Accuracy), textposition = "inside",
                     insidetextanchor = "end", textfont = list(color = "#ffffff"),
                     hovertemplate = "%{x}<br>Test accuracy %{y:.3f}<extra></extra>") |>
    ri_plotly(barmode = "group", bargap = 0.3,
              yaxis = list(title = "Accuracy", range = c(0, 1)), xaxis = list(title = ""),
              legend = list(orientation = "h", x = 0, y = 1.12))
}

# ==========================================================
# Class-level diagnostics
# ==========================================================

#' Confusion-matrix heatmap.
#'
#' @param d Dashboard data.
#' @param model Model name.
#' @param normalise If TRUE, show each actual class as row percentages.
#' @return plotly object.
plot_confusion <- function(d, model, normalise = FALSE) {
  ev <- d$evaluation[[model]]
  if (is.null(ev)) return(empty_plot("This model could not be loaded (is its R package installed?)."))
  cm <- ev$confusion
  fams <- d$families
  m <- matrix(0, length(fams), length(fams), dimnames = list(fams, fams)) # rows = actual
  for (i in seq_len(nrow(cm))) m[cm$Actual[i], cm$Predicted[i]] <- cm$Count[i]
  counts <- m
  if (normalise) {
    m <- sweep(m, 1, pmax(rowSums(m), 1), "/")
  }
  label <- if (normalise) matrix(ifelse(counts == 0, "", sprintf("%.0f%%", 100 * m)), nrow(m)) else
    matrix(ifelse(counts == 0, "", as.character(counts)), nrow(m))

  plotly::plot_ly(
    x = fams, y = fams, z = m, type = "heatmap", colorscale = seq_colorscale(), xgap = 2, ygap = 2,
    customdata = counts, text = label, texttemplate = "%{text}",
    hovertemplate = "Actual %{y}<br>Predicted %{x}<br>%{customdata} records<extra></extra>",
    colorbar = list(title = list(text = if (normalise) "Row %" else "Records", font = list(color = RI_COLORS$text_2)),
                    tickfont = list(color = RI_COLORS$muted), outlinewidth = 0,
                    tickformat = if (normalise) ".0%" else "")
  ) |>
    ri_plotly(show_legend = FALSE,
              xaxis = list(title = "Predicted family", tickangle = -40),
              yaxis = list(title = "Actual family", autorange = "reversed"))
}

# ==========================================================
# Helpers
# ==========================================================

#' Placeholder plot carrying a message.
#'
#' @param message Text to show.
#' @return plotly object.
empty_plot <- function(message) {
  plotly::plot_ly() |>
    plotly::layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(list(text = message, showarrow = FALSE, x = 0.5, y = 0.5,
                              xref = "paper", yref = "paper", font = list(color = RI_COLORS$muted)))
    ) |>
    ri_plotly(show_legend = FALSE)
}
