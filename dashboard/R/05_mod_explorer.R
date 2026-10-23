############################################################
# RNAInsight Dashboard
# 05_mod_explorer.R : Ligand Structure Explorer section
############################################################
#
# Written entirely in R:
#   * 3D structure      - r3dmol (R htmlwidget)
#   * 2D structure      - ggplot2 (plot_structure_2d())
#   * descriptor list   - shiny::radioButtons
#   * calculation panel - shiny::renderUI from feature_calculation()
#   * ligand list       - DT table with row selection
#   * "where this ligand sits" histogram - plotly
# The chemistry behind every highlight lives in 07_structure.R.
############################################################

#' Choices for the ligand picker: "Name - families" -> key.
#'
#' @param ligs Ligand table.
#' @return Named vector.
ligand_choices <- function(ligs) {
  ligs <- ligs[order(ligs$Name), ]
  stats::setNames(ligs$key, paste0(ligs$Name, "  ·  ", ligs$Families))
}

#' Radio-button choices for the descriptor list of one ligand.
#'
#' Each choice shows the ligand's value next to the descriptor name.
#'
#' @param mol Prepared molecule.
#' @return list(names = list of tags, values = character vector)
descriptor_choices <- function(mol) {
  d <- mol$descriptors
  row <- function(value, label, cls = "") {
    htmltools::span(class = paste("qx-desc-row", cls),
                    htmltools::span(class = "qx-desc-val", value),
                    htmltools::span(class = "qx-desc-name", label))
  }
  names_raw <- list(row("—", "Structure only (element colours)", "qx-desc-none"))
  values <- "none"
  for (f in names(FEATURE_INFO)) {
    info <- FEATURE_INFO[[f]]
    names_raw[[length(names_raw) + 1]] <- row(format_value(d[[f]], min(info$digits, 2)), info$label)
    values <- c(values, f)
  }
  for (f in names(ENGINEERED_INFO)) {
    e <- ENGINEERED_INFO[[f]]
    names_raw[[length(names_raw) + 1]] <- row(format_value(e$calc(d), 3), e$label, "qx-desc-eng")
    values <- c(values, f)
  }
  list(names = names_raw, values = values)
}

# ----------------------------------------------------------
# UI
# ----------------------------------------------------------

#' UI of the Ligand Structure Explorer.
#'
#' @param id Module id.
#' @param d Dashboard data.
#' @param default_key Ligand shown first.
#' @return Tag list.
explorer_ui <- function(id, d, default_key) {
  ns <- shiny::NS(id)
  first <- descriptor_choices(d$molecules[[as.character(default_key)]])
  has_r3dmol <- requireNamespace("r3dmol", quietly = TRUE)

  similar_card <- function(k) {
    htmltools::div(
      class = "qx-sim-card",
      shiny::plotOutput(ns(paste0("sim_plot_", k)), height = "140px"),
      shiny::uiOutput(ns(paste0("sim_info_", k))),
      shiny::actionButton(ns(paste0("sim_open_", k)), "Open", class = "qx-btn qx-btn-small")
    )
  }

  htmltools::div(
    class = "qx-root",

    # ---- top bar ---------------------------------------------------------
    htmltools::div(
      class = "qx-topbar",
      htmltools::div(
        class = "qx-title",
        htmltools::div(
          htmltools::div(class = "qx-kicker", "Ligand structure explorer"),
          htmltools::div(class = "qx-name", shiny::textOutput(ns("name"), inline = TRUE))
        ),
      ),
      htmltools::div(
        class = "qx-controls",
        shiny::selectInput(ns("family"), NULL, c("All RNA families" = "", d$families), width = "190px"),
        shiny::selectizeInput(ns("ligand"), NULL, choices = ligand_choices(d$ligands), selected = default_key,
                              width = "380px", options = list(placeholder = "Search 152 ligands..."))
      )
    ),

    # ---- stage: toolbar | viewer | panels -------------------------------
    htmltools::div(
      class = "qx-stage",
      htmltools::div(
        class = "qx-toolbar",
        htmltools::div(class = "qx-tool-head", "View"),
        shiny::radioButtons(ns("view"), NULL,
                            choices = if (has_r3dmol) c("3D" = "3d", "2D" = "2d") else c("2D" = "2d"),
                            selected = if (has_r3dmol) "3d" else "2d", inline = TRUE),
        shiny::uiOutput(ns("view_tools")),
        htmltools::div(class = "qx-tool-sep"),
        shiny::actionButton(ns("clear"), "Clear highlight", class = "qx-btn")
      ),
      htmltools::div(
        class = "qx-viewer",
        htmltools::div(
          class = "qx-ident",
          htmltools::div(class = "qx-ident-label", "Ligand SMILES"),
          shiny::textOutput(ns("smiles"), container = htmltools::div),
          htmltools::div(class = "qx-ident-label", "Formula"),
          shiny::textOutput(ns("formula"), container = htmltools::div)
        ),
        shiny::uiOutput(ns("viewer"), class = "qx-canvas"),
        shiny::uiOutput(ns("legend"))
      ),
      htmltools::div(
        class = "qx-side",
        htmltools::div(
          class = "qx-panel",
          htmltools::div(class = "qx-panel-head", shiny::textOutput(ns("ranking_title"), inline = TRUE)),
          DT::DTOutput(ns("ranking"))
        ),
        htmltools::div(
          class = "qx-panel qx-dock",
          htmltools::div(class = "qx-panel-head", "Molecular descriptors",
                         htmltools::span(class = "qx-hint", "pick one to see how it is calculated · grey = engineered in R")),
          shiny::radioButtons(ns("feature"), NULL, choiceNames = first$names, choiceValues = first$values,
                              selected = "none")
        )
      )
    ),

    # ---- analysis strip ----------------------------------------------------
    htmltools::div(
      class = "qx-strip",
      htmltools::div(class = "qx-card qx-calc-card",
                     htmltools::div(class = "qx-card-head", "How it is calculated"),
                     shiny::uiOutput(ns("calc"))),
      htmltools::div(class = "qx-card",
                     htmltools::div(class = "qx-card-head", shiny::textOutput(ns("context_title"), inline = TRUE)),
                     plotly::plotlyOutput(ns("context"), height = "250px"),
                     htmltools::div(class = "qx-card-foot", shiny::textOutput(ns("context_text"))))
    ),
    htmltools::div(
      class = "qx-strip qx-strip-2",
      htmltools::div(class = "qx-card",
                     htmltools::div(class = "qx-card-head", "RNA targets this ligand was tested on (R-BIND records)"),
                     shiny::uiOutput(ns("targets"))),
      htmltools::div(class = "qx-card",
                     htmltools::div(class = "qx-card-head", "Most similar ligands (14-descriptor distance)"),
                     htmltools::div(class = "qx-similar", similar_card(1), similar_card(2), similar_card(3)))
    )
  )
}

# ----------------------------------------------------------
# Server
# ----------------------------------------------------------

#' Server of the Ligand Structure Explorer.
#'
#' @param id Module id.
#' @param d Dashboard data.
#' @param selected reactiveVal holding the selected ligand key
#'        (shared with the Results section).
explorer_server <- function(id, d, selected) {
  shiny::moduleServer(id, function(input, output, session) {
    ligs <- d$ligands
    scaled <- scale(ligs[, RAW_DESCRIPTORS])
    scaled[is.na(scaled)] <- 0
    has_r3dmol <- requireNamespace("r3dmol", quietly = TRUE)

    current <- shiny::reactive({
      key <- selected()
      shiny::req(key)
      ligs[ligs$key == key, ]
    })
    mol <- shiny::reactive(d$molecules[[as.character(selected())]])
    feature <- shiny::reactive(if (is.null(input$feature)) "none" else input$feature)
    highlight <- shiny::reactive(feature_highlight(mol(), feature()))

    # ---- selection plumbing ----
    shiny::observeEvent(input$ligand, {
      key <- suppressWarnings(as.integer(input$ligand))
      if (!is.na(key) && !identical(key, selected())) selected(key)
    }, ignoreInit = TRUE)

    shiny::observeEvent(selected(), {
      key <- selected()
      if (!identical(as.character(key), input$ligand)) {
        shiny::updateSelectizeInput(session, "ligand", selected = key)
      }
      ch <- descriptor_choices(d$molecules[[as.character(key)]])
      shiny::updateRadioButtons(session, "feature", choiceNames = ch$names, choiceValues = ch$values,
                                selected = feature())
    })

    shiny::observeEvent(input$family, {
      sub <- if (nzchar(input$family)) ligs[grepl(input$family, ligs$Families, fixed = TRUE), ] else ligs
      keep <- if (selected() %in% sub$key) selected() else sub$key[which.max(sub$N_Records)]
      shiny::updateSelectizeInput(session, "ligand", choices = ligand_choices(sub), selected = keep)
      selected(keep)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$clear, shiny::updateRadioButtons(session, "feature", selected = "none"))

    # ---- header ----
    output$name <- shiny::renderText(paste0(current()$Name, "  ·  ", current()$Primary_Family))
    output$smiles <- shiny::renderText(mol()$smiles)
    output$formula <- shiny::renderText(sprintf("%s    %d heavy atoms · %d atoms with H",
                                                mol()$formula, mol()$n_heavy, length(mol()$el)))

    # ---- view-specific controls ----
    # 3D view: style, hydrogens, rotation. 2D view: atom names, save image.
    choice_or <- function(value, default) if (is.null(value)) default else value
    output$view_tools <- shiny::renderUI({
      ns <- session$ns
      if (!identical(input$view, "3d")) {
        return(htmltools::tagList(
          htmltools::div(class = "qx-tool-head", "Atom names"),
          shiny::radioButtons(ns("names"), NULL, c("Off" = "off", "On" = "on"),
                              selected = shiny::isolate(choice_or(input$names, "off")), inline = TRUE),
          shiny::downloadButton(ns("download"), "Save 2D image", class = "qx-btn")
        ))
      }
      htmltools::tagList(
        htmltools::div(class = "qx-tool-head", "3D style"),
        shiny::selectInput(ns("style"), NULL, c("Ball & stick" = "ballstick", "Sticks" = "stick",
                                                "Space filling" = "sphere"),
                           selected = shiny::isolate(choice_or(input$style, "ballstick")), width = "100%"),
        htmltools::div(class = "qx-tool-head", "Hydrogens"),
        shiny::radioButtons(ns("show_h"), NULL, c("Show" = "show", "Hide" = "hide"),
                            selected = shiny::isolate(choice_or(input$show_h, "show")), inline = TRUE),
        htmltools::div(class = "qx-tool-head", "Rotation"),
        shiny::radioButtons(ns("spin"), NULL, c("Off" = "off", "On" = "on"),
                            selected = shiny::isolate(choice_or(input$spin, "off")), inline = TRUE)
      )
    })

    # ---- structure viewer ----
    output$viewer <- shiny::renderUI({
      if (identical(input$view, "3d") && has_r3dmol) {
        r3dmol::r3dmolOutput(session$ns("v3d"), height = "600px")
      } else {
        shiny::plotOutput(session$ns("v2d"), height = "600px")
      }
    })
    if (has_r3dmol) {
      output$v3d <- r3dmol::renderR3dmol({
        # depending on the view makes the viewer redraw each time 3D is shown again
        shiny::req(identical(input$view, "3d"))
        structure_3d(mol(), highlight(), style = choice_or(input$style, "ballstick"),
                     show_h = !identical(input$show_h, "hide"),
                     show_labels = FALSE,
                     spin = identical(input$spin, "on"))
      })
    }
    output$v2d <- shiny::renderPlot({
      shiny::req(identical(input$view, "2d"))
      plot_structure_2d(mol(), highlight(), show_names = identical(input$names, "on"))
    }, bg = "transparent", res = 110)

    output$download <- shiny::downloadHandler(
      filename = function() {
        paste0(gsub("[^A-Za-z0-9_-]+", "_", mol()$name), "_", feature(), ".png")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot_structure_2d(mol(), highlight(), show_names = identical(input$names, "on")),
                        width = 10, height = 6.5, dpi = 200, bg = RI_COLORS$page)
      }
    )

    output$legend <- shiny::renderUI({
      h <- highlight()
      item <- function(col, label) {
        htmltools::span(class = "qx-legend-item", htmltools::tags$i(style = paste0("background:", col)), label)
      }
      if (feature() == "none") {
        els <- c("C", "N", "O", "S", "F", "H")
        return(htmltools::div(class = "qx-legend", htmltools::div(class = "qx-legend-title", "Element colours"),
                              lapply(els, function(e) item(element_color(e), e))))
      }
      title <- if (feature() %in% names(FEATURE_INFO)) FEATURE_INFO[[feature()]]$label else
        ENGINEERED_INFO[[feature()]]$label
      htmltools::div(
        class = "qx-legend",
        htmltools::div(class = "qx-legend-title", paste("Showing:", title)),
        lapply(seq_len(nrow(h$legend)), function(k) item(h$legend$color[k], h$legend$label[k])),
        item(HL$dim, "not counted")
      )
    })

    # ---- calculation panel ----
    output$calc <- shiny::renderUI({
      if (feature() == "none") {
        return(htmltools::div(
          htmltools::div(class = "qx-calc-name", "Pick a descriptor"),
          htmltools::p(class = "qx-calc-what",
                       "Pick any of the 14 RDKit descriptors (or the 6 engineered features) in the panel above. ",
                       "The atoms, bonds or rings it counts light up on the structure, and this panel shows the ",
                       "calculation step by step, checked against the value stored in the dataset.")
        ))
      }
      calc <- feature_calculation(mol(), feature())
      value_txt <- format_value(calc$value, calc$digits)
      table_tag <- NULL
      if (!is.null(calc$table)) {
        tab <- calc$table
        table_tag <- htmltools::div(
          class = "qx-table-wrap",
          htmltools::tags$table(
            class = "qx-table",
            htmltools::tags$thead(htmltools::tags$tr(lapply(names(tab), htmltools::tags$th))),
            htmltools::tags$tbody(
              if (nrow(tab)) lapply(seq_len(nrow(tab)), function(r) {
                htmltools::tags$tr(lapply(tab[r, ], function(x) htmltools::tags$td(as.character(x))))
              }) else htmltools::tags$tr(htmltools::tags$td(colspan = ncol(tab), class = "qx-none", "none in this ligand"))
            )
          )
        )
      }
      htmltools::tagList(
        htmltools::div(
          class = "qx-calc-top",
          htmltools::div(class = "qx-calc-value", value_txt,
                         if (nzchar(calc$unit)) htmltools::span(class = "qx-unit", calc$unit)),
          htmltools::div(htmltools::div(class = "qx-calc-name", calc$label),
                         htmltools::tags$code(class = "qx-rdkit", calc$rdkit))
        ),
        htmltools::p(class = "qx-calc-what", calc$what),
        htmltools::div(class = "qx-formula", calc$formula),
        if (!is.null(calc$worked)) htmltools::div(class = "qx-formula sub",
                                                  paste(calc$worked, "=", format_value(calc$value, 3))),
        table_tag,
        if (!calc$engineered) htmltools::div(
          class = paste("qx-check", if (calc$ok) "ok" else "warn"),
          sprintf("%s Breakdown total %s = stored dataset value %s", if (calc$ok) "✓" else "!",
                  format_value(calc$total, calc$digits), value_txt)
        )
      )
    })

    # ---- ligand list (DT, click a row to open it) ----
    ranking_family <- shiny::reactive({
      if (!is.null(input$family) && nzchar(input$family)) input$family else current()$Primary_Family
    })
    ranking_table <- shiny::reactive({
      fam <- ranking_family()
      sub <- ligs[grepl(fam, ligs$Families, fixed = TRUE), ]
      sub[order(-sub$N_Records, sub$Name), ]
    })
    output$ranking_title <- shiny::renderText(paste0("Ligands binding ", ranking_family()))
    output$ranking <- DT::renderDT({
      sub <- ranking_table()
      show <- data.frame(`#` = seq_len(nrow(sub)), Ligand = sub$Name, Records = sub$N_Records, check.names = FALSE)
      DT::datatable(
        show, rownames = FALSE, class = "compact ri-dt qx-rank",
        selection = list(mode = "single", selected = which(sub$key == shiny::isolate(selected()))),
        options = list(dom = "t", paging = FALSE, scrollY = "190px", ordering = FALSE,
                       columnDefs = list(list(width = "28px", targets = 0)))
      )
    })
    shiny::observeEvent(input$ranking_rows_selected, {
      sub <- ranking_table()
      key <- sub$key[input$ranking_rows_selected]
      if (length(key) && !identical(key, selected())) selected(key)
    })

    # ---- feature context ----
    context_feature <- shiny::reactive({
      f <- feature()
      if (!(f %in% names(ligs))) "MolecularWeight" else f
    })
    output$context_title <- shiny::renderText(paste0("Where this ligand sits: ", feature_label(context_feature())))
    output$context <- plotly::renderPlotly({
      f <- context_feature()
      lig <- current()
      v <- lig[[f]]
      plotly::plot_ly(x = ligs[[f]], type = "histogram", nbinsx = 30,
                      marker = list(color = RI_COLORS$neutral, line = list(color = RI_COLORS$card, width = 1)),
                      hovertemplate = "%{x}: %{y} ligands<extra></extra>", name = "All ligands") |>
        plotly::layout(shapes = list(list(type = "line", x0 = v, x1 = v, y0 = 0, y1 = 1, yref = "paper",
                                          line = list(color = RI_COLORS$accent, width = 3))),
                       annotations = list(list(x = v, y = 1, yref = "paper", text = lig$Name, showarrow = FALSE,
                                               yanchor = "bottom", font = list(color = RI_COLORS$accent, size = 11)))) |>
        ri_plotly(show_legend = FALSE, xaxis = list(title = feature_label(f)), yaxis = list(title = "Ligands"),
                  margin = list(l = 10, r = 10, t = 24, b = 10)) |>
        plotly::config(displayModeBar = FALSE)
    })
    output$context_text <- shiny::renderText({
      f <- context_feature()
      lig <- current()
      sprintf("%s = %s, higher than %.0f%% of the 152 ligands.", f, format(signif(lig[[f]], 5)),
              100 * mean(ligs[[f]] < lig[[f]]))
    })

    # ---- RNA targets table ----
    output$targets <- shiny::renderUI({
      lig <- current()
      rows <- d$cleaned[unlist(lig$Rows), c("X.", "Target_Family", "Assay")]
      rows <- rows[order(rows$Target_Family == "Other", as.character(rows$Target_Family)), ]
      htmltools::div(
        class = "qx-targets",
        htmltools::tags$table(
          class = "qx-table",
          htmltools::tags$thead(htmltools::tags$tr(htmltools::tags$th("RNA family"), htmltools::tags$th("Assay target (as reported)"))),
          htmltools::tags$tbody(lapply(seq_len(nrow(rows)), function(k) {
            assay <- rows$Assay[k]
            htmltools::tags$tr(
              htmltools::tags$td(htmltools::span(class = "qx-fam", as.character(rows$Target_Family[k]))),
              htmltools::tags$td(if (is.na(assay) || !nzchar(assay)) htmltools::span(class = "qx-none", "(no specific target reported → 'Other')") else assay)
            )
          }))
        )
      )
    })

    # ---- similar ligands ----
    similar <- shiny::reactive({
      i <- which(ligs$key == selected())
      dist <- sqrt(rowSums((scaled - matrix(scaled[i, ], nrow(scaled), ncol(scaled), byrow = TRUE))^2))
      dist[i] <- Inf
      top <- order(dist)[1:3]
      data.frame(key = ligs$key[top], dist = dist[top], name = ligs$Name[top],
                 families = ligs$Families[top], formula = ligs$Formula[top])
    })
    lapply(1:3, function(k) {
      output[[paste0("sim_plot_", k)]] <- shiny::renderPlot({
        plot_structure_2d(d$molecules[[as.character(similar()$key[k])]], NULL, mini = TRUE)
      }, bg = "transparent", res = 96)
      output[[paste0("sim_info_", k)]] <- shiny::renderUI({
        s <- similar()[k, ]
        htmltools::tagList(
          htmltools::div(class = "qx-sim-name", s$name),
          htmltools::div(class = "qx-sim-meta", s$families),
          htmltools::div(class = "qx-sim-row", htmltools::span(s$formula), htmltools::span(sprintf("distance %.2f", s$dist)))
        )
      })
      shiny::observeEvent(input[[paste0("sim_open_", k)]], selected(similar()$key[k]))
    })
  })
}
