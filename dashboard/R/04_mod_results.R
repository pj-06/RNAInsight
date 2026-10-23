############################################################
# RNAInsight Dashboard
# 04_mod_results.R : "Results & Visualizations" section
############################################################
#
# Shiny module with five tabs:
#   Overview | Data Exploration | Feature Engineering |
#   Model Comparison | Class Diagnostics
#
# Every chart's controls sit in that chart's own card header,
# so each dropdown is next to the plot it changes. All charts
# come from 03_plots.R; this file only lays out the cards and
# wires inputs to outputs.
############################################################

# ----------------------------------------------------------
# Layout helpers
# ----------------------------------------------------------

#' Small pill naming the package that draws a chart.
#'
#' @param ... Package names, e.g. "ggplot2", "plotly".
#' @return Tag list of badges.
tool_badges <- function(...) {
  htmltools::tags$span(
    class = "ri-badges",
    lapply(c(...), function(t) htmltools::tags$span(class = paste("ri-badge", paste0("ri-badge-", t)), t))
  )
}

#' Dashboard card with a title, tool badges, its own controls and a note.
#'
#' @param title Card title.
#' @param body Card content (usually an output).
#' @param tools Character vector of packages used.
#' @param note Optional one-line "how to read this" text.
#' @param controls Optional inputs shown in the card header, next to the chart.
#' @param height Optional card body height (CSS).
#' @param full_screen Allow the bslib full-screen toggle.
#' @param fill Let the card body stretch to fill its container (FALSE for
#'   tables, so the card is only as tall as its rows).
#' @return bslib card.
ri_card <- function(title, body, tools = NULL, note = NULL, controls = NULL,
                    height = NULL, full_screen = TRUE, fill = TRUE) {
  bslib::card(
    full_screen = full_screen,
    fill = fill,
    class = "ri-card",
    bslib::card_header(
      htmltools::div(
        class = "ri-card-head",
        htmltools::div(class = "ri-card-title", title, if (length(tools)) do.call(tool_badges, as.list(tools))),
        if (!is.null(controls)) htmltools::div(class = "ri-card-controls", controls)
      )
    ),
    bslib::card_body(body, height = height, fill = fill, fillable = fill, class = "ri-card-body"),
    if (!is.null(note)) bslib::card_footer(class = "ri-card-note", note)
  )
}

#' KPI tile.
#'
#' @param value Headline value.
#' @param label What the value is.
#' @param detail Optional smaller detail line.
#' @return div tag.
kpi_tile <- function(value, label, detail = NULL) {
  htmltools::div(
    class = "ri-kpi",
    htmltools::div(class = "ri-kpi-label", label),
    htmltools::div(class = "ri-kpi-value", value),
    if (!is.null(detail)) htmltools::div(class = "ri-kpi-detail", detail)
  )
}

#' Compact labelled select for a card header.
#'
#' The label sits to the left of the dropdown, so several controls
#' fit on one line beside the card title.
#'
#' @param inputId Input id.
#' @param label Short label (or NULL).
#' @param choices Choices.
#' @param selected Selected value.
#' @param width CSS width of the dropdown.
#' @return Tag.
mini_select <- function(inputId, label, choices, selected = NULL, width = "170px") {
  htmltools::div(
    class = "ri-mini-input",
    if (!is.null(label)) htmltools::tags$label(class = "ri-mini-label", `for` = inputId, label),
    shiny::selectInput(inputId, NULL, choices = choices, selected = selected, width = width)
  )
}

#' UI of the Results & Visualizations section.
#'
#' @param id Module id.
#' @param d Dashboard data (for choices and static text).
#' @return Tag list.
results_ui <- function(id, d) {
  ns <- shiny::NS(id)
  feature_choices <- stats::setNames(c(RAW_DESCRIPTORS, ENGINEERED_FEATURES),
                                     feature_label(c(RAW_DESCRIPTORS, ENGINEERED_FEATURES)))
  family_choices <- c("None" = HIGHLIGHT_NONE, d$families)
  models_available <- intersect(MODEL_ORDER, d$leaderboard$Model)
  model_choices <- stats::setNames(models_available, MODEL_SHORT[models_available])
  best <- d$leaderboard[which.max(d$leaderboard$Accuracy), ]

  bslib::navset_underline(
    id = ns("tabs"),

    # ---------------- Overview ----------------
    bslib::nav_panel(
      "Overview",
      htmltools::div(
        class = "ri-hero",
        htmltools::div(class = "ri-hero-kicker", "RNA target family prediction from ligand descriptors"),
        htmltools::h2("Can a ligand's molecular descriptors predict which RNA family it binds?"),
        htmltools::p("RNAInsight classifies RNA-binding small molecules from the R-BIND database into 11 RNA target ",
                     "families using 14 RDKit descriptors. This section collects the results and visualizations ",
                     "of the pipeline; the Ligand Structure Explorer shows how each descriptor is computed on the molecule.")
      ),
      htmltools::div(
        class = "ri-kpi-row",
        kpi_tile(format(nrow(d$cleaned)), "Molecule-Target Records", "After cleaning and deduplication"),
        kpi_tile(format(nrow(d$ligands)), "Unique Ligands", "Each tested against 1 to 11 RNAs"),
        kpi_tile(format(length(d$families)), "RNA Target Families", "8 to 288 records per family"),
        kpi_tile("14 + 6", "Molecular Descriptors", "RDKit-computed + engineered"),
        kpi_tile(format(nrow(d$leaderboard)), "Classification Models", "Identical split and CV folds"),
        kpi_tile(sprintf("%.1f%%", 100 * best$Accuracy), "Best Test Accuracy", MODEL_SHORT[[best$Model]])
      ),
      htmltools::div(
        class = "ri-pipeline",
        lapply(list(
          c("1", "Data Source", "R-BIND dataset"),
          c("2", "Descriptor Generation", "Python + RDKit"),
          c("3", "Data Cleaning", "750 records, 11 families"),
          c("4", "Feature Engineering", "6 derived, 15 selected"),
          c("5", "Database Storage", "SQLite"),
          c("6", "Model Training", paste(nrow(d$leaderboard), "classifiers with caret")),
          c("7", "Dashboard", "R Shiny")
        ), function(s) htmltools::div(class = "ri-step",
                                      htmltools::span(class = "ri-step-n", s[1]),
                                      htmltools::div(htmltools::tags$b(s[2]), htmltools::br(), htmltools::span(s[3]))))
      ),
      bslib::layout_columns(
        col_widths = c(6, 6),
        ri_card("Records vs unique ligands per family",
                plotly::plotlyOutput(ns("composition"), height = "420px"),
                tools = "plotly",
                note = "750 records come from only 152 distinct structures: most ligands were assayed against more than one RNA."),
        ri_card("Ligands shared between RNA families",
                plotly::plotlyOutput(ns("overlap"), height = "420px"),
                tools = "plotly",
                note = "Each cell counts ligands labelled with both families. The diagonal is the family's own ligand count. 'Other' shares ligands with every family, which is why it is the hardest class to predict.")
      )
    ),

    # ---------------- Data exploration ----------------
    bslib::nav_panel(
      "Data Exploration",
      bslib::layout_columns(
        col_widths = c(5, 7),
        ri_card("Distribution", plotly::plotlyOutput(ns("eda_hist"), height = "420px"),
                tools = c("ggplot2", "plotly"),
                controls = htmltools::tagList(
                  mini_select(ns("hist_feature"), "Descriptor", feature_choices, "LogP", "190px"),
                  mini_select(ns("hist_highlight"), "Highlight", family_choices, "Riboswitch", "150px")
                ),
                note = "Highlighted family in blue, drawn over all other records in grey (overlaid, not stacked)."),
        ri_card("Ligand scatter", plotly::plotlyOutput(ns("eda_scatter"), height = "420px"),
                tools = "plotly",
                controls = htmltools::tagList(
                  mini_select(ns("sc_x"), "x", feature_choices, "LogP", "150px"),
                  mini_select(ns("sc_y"), "y", feature_choices, "TPSA", "150px"),
                  mini_select(ns("sc_size"), "Size", c("Fixed" = "", "Records" = "N_Records", feature_choices),
                              "N_Records", "140px"),
                  mini_select(ns("sc_highlight"), "Highlight", family_choices, "Riboswitch", "140px")
                ),
                note = "One point per unique ligand (152). Click a point to open that ligand in the Structure Explorer.")
      )
    ),

    # ---------------- Feature engineering ----------------
    bslib::nav_panel(
      "Feature Engineering",
      bslib::layout_columns(
        col_widths = c(6, 6),
        ri_card("Feature selection: why 15 of 20 candidates were kept",
                plotly::plotlyOutput(ns("fe_selection"), height = "520px"),
                tools = "plotly",
                note = "Dashed line = the |r| > 0.90 cut-off in 04_feature_engineering.R. From each correlated pair, caret::findCorrelation() drops one feature and keeps the other (e.g. drops MW and keeps MolMR). No feature had near-zero variance."),
        ri_card("Model-based feature importance",
                plotly::plotlyOutput(ns("fe_importance_model"), height = "520px"),
                tools = "plotly", note = "caret::varImp() of the tree ensembles, scaled 0-100 within each model.")
      ),
      ri_card(
        "The six engineered features", full_screen = FALSE,
        htmltools::div(
          class = "ri-formula-grid",
          formula_tile("HBondTotal", "HBA + HBD", "total hydrogen-bonding capacity"),
          formula_tile("LipinskiViolations", "[MW>500] + [LogP>5] + [HBA>10] + [HBD>5]", "drug-likeness"),
          formula_tile("AromaticRingFraction", "NumAromaticRings / RingCount", "0 when there are no rings"),
          formula_tile("FlexibilityIndex", "RotatableBonds / HeavyAtomCount", "flexibility per atom"),
          formula_tile("HeteroatomFraction", "NumHeteroAtoms / HeavyAtomCount", "heteroatom density"),
          formula_tile("TPSA_per_MW", "TPSA / MolecularWeight", "polarity normalised by size")
        )
      )
    ),

    # ---------------- Model comparison ----------------
    bslib::nav_panel(
      "Model Comparison",
      htmltools::div(
        class = "ri-section-intro compact",
        htmltools::p(nrow(d$leaderboard), " algorithms trained with caret on the same stratified 80/20 split. ",
                     "Cross-validated accuracy is the mean over the same 5-fold \u00d7 3 cross-validation folds of the ",
                     "training data; test accuracy is measured once on the ", d$test_size, " held-out test records.")
      ),
      ri_card("Cross-validated vs test accuracy", plotly::plotlyOutput(ns("mc_profile"), height = "440px"),
              tools = "plotly",
              note = "Models sorted by test accuracy. Error bars show the standard deviation across the 15 CV folds; similar CV and test accuracy means the model is not overfitting."),
      ri_card("Leaderboard", DT::DTOutput(ns("mc_table"), height = "auto"), tools = "DT",
              full_screen = FALSE, fill = FALSE)
    ),

    # ---------------- Class diagnostics ----------------
    bslib::nav_panel(
      "Class Diagnostics",
      ri_card("Confusion matrix", plotly::plotlyOutput(ns("cd_confusion"), height = "600px"),
              tools = "plotly",
              controls = htmltools::tagList(
                mini_select(ns("cm_model"), "Model", model_choices, models_available[1], "190px"),
                mini_select(ns("cm_normalise"), "Show", c("Counts" = "counts", "Row %" = "percent"), "counts", "120px")
              ),
              note = "Rows = actual family, columns = prediction; the diagonal holds the correct predictions. Most errors are records labelled 'Other' whose ligand also has a specific family.")
    )
  )
}

#' One engineered-feature formula tile.
#'
#' @param name Feature name.
#' @param formula Formula text.
#' @param note Short explanation.
#' @return Tag.
formula_tile <- function(name, formula, note) {
  htmltools::div(class = "ri-formula", htmltools::div(class = "ri-formula-name", name),
                 htmltools::code(formula), htmltools::div(class = "ri-formula-note", note))
}

# ----------------------------------------------------------
# Server
# ----------------------------------------------------------

#' Server of the Results & Visualizations section.
#'
#' @param id Module id.
#' @param d Dashboard data.
#' @param open_ligand Callback(key) that opens a ligand in the explorer.
results_server <- function(id, d, open_ligand) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Overview ----
    output$composition <- plotly::renderPlotly(plot_family_composition(d))
    output$overlap <- plotly::renderPlotly(plot_family_overlap(d))

    # ---- Data exploration ----
    output$eda_hist <- plotly::renderPlotly(plot_distribution(d, input$hist_feature, input$hist_highlight))
    output$eda_scatter <- plotly::renderPlotly({
      plot_scatter(d, input$sc_x, input$sc_y, input$sc_highlight, input$sc_size) |>
        plotly::event_register("plotly_click")
    })

    # suppressWarnings: the event is queried before the scatter's tab has
    # been opened (and the plot registered) the first time.
    scatter_click <- shiny::reactive(suppressWarnings(plotly::event_data("plotly_click", source = "scatter")))
    shiny::observeEvent(scatter_click(), {
      click <- scatter_click()
      if (!is.null(click$key)) open_ligand(as.integer(click$key[[1]]))
    })

    # ---- Feature engineering ----
    output$fe_selection <- plotly::renderPlotly(plot_feature_selection(d))
    output$fe_importance_model <- plotly::renderPlotly(plot_model_importance(d))

    # ---- Model comparison ----
    output$mc_profile <- plotly::renderPlotly(plot_accuracy_comparison(d))
    output$mc_table <- DT::renderDT({
      lb <- d$leaderboard[order(-d$leaderboard$Accuracy), ]
      show <- data.frame(
        Model = lb$Model,
        `Test accuracy` = lb$Accuracy,
        `CV accuracy (mean)` = lb$CV_Accuracy,
        `CV accuracy (SD)` = lb$CV_Accuracy_SD,
        check.names = FALSE
      )
      DT::datatable(show, rownames = FALSE, class = "compact ri-dt",
                    options = list(dom = "t", pageLength = 20, ordering = TRUE)) |>
        DT::formatRound(2:4, 3)
    })

    # ---- Class diagnostics ----
    output$cd_confusion <- plotly::renderPlotly(
      plot_confusion(d, input$cm_model, identical(input$cm_normalise, "percent"))
    )
  })
}
