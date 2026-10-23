############################################################
# RNAInsight Dashboard
# 08_styles.R : page styles
############################################################
#
# The dashboard's look (dark theme, card layout, explorer
# panels) as one style sheet held in an R string, added to the
# page with htmltools::tags$style() in app.R, so the project has
# no separate front-end files. Colours mirror RI_COLORS in
# 01_theme.R.
############################################################

#' Style rules for the whole dashboard.
RI_STYLES <- r"---[
/* =========================================================================
 * RNAInsight dashboard - styles
 * One dark visual system (tokens mirror dashboard/R/01_theme.R).
 * ======================================================================= */
:root {
  --page: #0b0f0e;
  --surface: #121816;
  --card: #161d1b;
  --card-alt: #1b2421;
  --border: #26302c;
  --border-soft: rgba(255, 255, 255, 0.07);
  --text: #eef1ef;
  --text-2: #b9c2bd;
  --muted: #84908a;
  --accent: #3ecf8e;
  --accent-dim: #1f6b4c;
  --blue: #3987e5;
  --orange: #d95926;
  --radius: 14px;
  --dock-w: 640px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: Consolas, "Cascadia Mono", monospace;
  color-scheme: dark;
}

html, body { background: var(--page); color: var(--text); font-family: var(--font); }
body { font-size: 14px; }
a { color: var(--accent); }
code, pre { font-family: var(--mono); }

/* ---------------- navbar ---------------- */
.navbar { background: rgba(11, 15, 14, 0.94) !important; border-bottom: 1px solid var(--border);
  backdrop-filter: blur(8px); padding: 0 !important; min-height: 64px; }
.navbar > .container-fluid { max-width: 1680px; margin: 0 auto; padding: 0 28px; gap: 28px; min-height: 64px;
  align-items: center; flex-wrap: nowrap; }
.navbar-brand { margin-right: 8px; padding: 0; display: flex; align-items: center; height: 64px; }
.navbar .navbar-collapse { align-self: stretch; }
.navbar .navbar-nav { gap: 6px; align-items: center; height: 100%; }
.navbar .nav-link { color: var(--text-2) !important; font-weight: 500; font-size: 0.95rem; padding: 21px 14px !important;
  border-radius: 0; }
.navbar .nav-link.active, .navbar .nav-link:hover { color: var(--text) !important; }
.navbar .nav-link.active { box-shadow: inset 0 -2px 0 var(--accent); }
.ri-brand { display: flex; align-items: center; gap: 1px; padding-right: 26px; border-right: 1px solid var(--border); line-height: 1; }
.ri-brand-mark { font-weight: 800; letter-spacing: 0.02em; color: var(--accent); font-size: 1.4rem; }
.ri-brand-name { font-weight: 600; font-size: 1.4rem; color: var(--text); }
.ri-nav-link { color: var(--muted) !important; font-size: 0.85rem; text-decoration: none; padding: 0 4px; }
.ri-nav-link:hover { color: var(--text) !important; }
@media (max-width: 990px) { .ri-brand { border-right: 0; padding-right: 0; } .navbar .nav-link { padding: 10px 6px !important; } }

/* ---------------- overview hero ---------------- */
.ri-hero { padding: 6px 0 14px; }
.ri-hero-kicker { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.12em; color: var(--accent); font-weight: 600; margin-bottom: 6px; }
.ri-hero h2 { font-size: 1.5rem; font-weight: 700; margin: 0 0 6px; letter-spacing: -0.01em; }
.ri-hero p { color: var(--text-2); max-width: 980px; margin: 0; }

.ri-page { padding: 20px 28px 40px; max-width: 1680px; margin: 0 auto; }
@media (max-width: 600px) { .ri-page { padding: 12px 16px 32px; } }

/* ---------------- section tabs ---------------- */
.ri-page > .nav-underline, .ri-page .nav-underline { border-bottom: 1px solid var(--border); margin-bottom: 18px; gap: 6px; flex-wrap: wrap; }
.nav-underline .nav-link { color: var(--muted); font-weight: 500; padding: 10px 6px; }
.nav-underline .nav-link.active { color: var(--text); border-bottom-color: var(--accent); }

.ri-section-intro h2 { font-size: 1.45rem; font-weight: 700; margin: 4px 0 6px; letter-spacing: -0.01em; }
.ri-section-intro p { color: var(--text-2); max-width: 980px; margin-bottom: 14px; }
.ri-section-intro.compact p { margin-bottom: 12px; }

/* ---------------- KPI tiles ---------------- */
.ri-kpi-row { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 12px; margin-bottom: 14px; }
@media (max-width: 1200px) { .ri-kpi-row { grid-template-columns: repeat(3, minmax(0, 1fr)); } }
@media (max-width: 600px) { .ri-kpi-row { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
.ri-kpi { background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); padding: 16px 18px;
  border-top: 2px solid var(--accent-dim); }
.ri-kpi-label { color: var(--muted); font-size: 0.7rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.08em; }
.ri-kpi-value { font-size: 2rem; font-weight: 700; color: var(--text); line-height: 1.15; margin: 6px 0 4px;
  font-variant-numeric: tabular-nums; }
.ri-kpi-detail { color: var(--text-2); font-size: 0.78rem; }

/* ---------------- pipeline ---------------- */
.ri-pipeline { display: grid; grid-template-columns: repeat(7, minmax(0, 1fr)); gap: 8px; margin-bottom: 16px; }
@media (max-width: 1100px) { .ri-pipeline { grid-template-columns: repeat(4, minmax(0, 1fr)); } }
@media (max-width: 600px) { .ri-pipeline { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
.ri-step { display: flex; gap: 10px; align-items: flex-start; background: var(--surface); border: 1px solid var(--border);
  border-radius: 12px; padding: 10px 12px; font-size: 0.8rem; color: var(--muted); }
.ri-step b { color: var(--text); font-weight: 600; font-size: 0.85rem; }
.ri-step-n { flex: none; width: 22px; height: 22px; border-radius: 50%; background: var(--accent-dim); color: var(--accent);
  display: grid; place-items: center; font-weight: 700; font-size: 0.75rem; }

/* ---------------- cards ---------------- */
.ri-card.card { background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); margin-bottom: 16px; }
.ri-card .card-header { background: transparent; border-bottom: 1px solid var(--border-soft); padding: 12px 16px; }
.ri-card .card-footer { background: transparent; border-top: 1px solid var(--border-soft); }
.ri-card-head { display: flex; justify-content: space-between; align-items: center; gap: 10px; flex-wrap: wrap; }
.ri-card-title { font-weight: 600; color: var(--text); display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.ri-card-controls { display: flex; gap: 8px; align-items: center; }
.ri-card-note { color: var(--muted); font-size: 0.8rem; padding: 8px 16px; }
.ri-card-body { padding: 8px 12px; }
.ri-stack { display: flex; flex-direction: column; }
.ri-badges { display: inline-flex; gap: 4px; }
.ri-badge { font-size: 0.68rem; font-weight: 600; padding: 2px 7px; border-radius: 999px; border: 1px solid var(--border);
  color: var(--muted); font-family: var(--mono); }
.ri-badge-ggplot2 { color: #86b6ef; border-color: #1c3a5e; }
.ri-badge-plotly { color: #3ecf8e; border-color: #1f4a39; }
.ri-badge-DT { color: #c8a4ff; border-color: #3b2c59; }
.ri-missing { color: var(--muted); padding: 40px; text-align: center; }

/* ---------------- inputs ---------------- */
.form-control, .form-select, .selectize-input, .selectize-dropdown {
  background: var(--surface) !important; color: var(--text) !important; border-color: var(--border) !important;
}
.selectize-input { box-shadow: none !important; }
.selectize-input > input { color: var(--text) !important; }
.selectize-dropdown .option { color: var(--text-2); }
.selectize-dropdown .active { background: var(--card-alt) !important; color: var(--text) !important; }
.selectize-input .item { color: var(--text); }
.selectize-control.multi .selectize-input > div { background: var(--card-alt); color: var(--text); border-radius: 6px; }
.control-label, .form-label, .shiny-input-container > label { color: var(--muted); font-size: 0.78rem; font-weight: 500; margin-bottom: 4px; }
.ri-card-controls { flex-wrap: wrap; justify-content: flex-end; gap: 10px; }
.ri-mini-input { display: flex; align-items: center; gap: 6px; }
.ri-mini-input .shiny-input-container { margin-bottom: 0; }
.ri-mini-label { color: var(--muted); font-size: 0.75rem; font-weight: 500; margin: 0; white-space: nowrap; }
.ri-mini-check .form-check, .ri-mini-check .checkbox { margin: 0; }
.ri-mini-check label { color: var(--text-2); font-size: 0.8rem; }
.ri-mini-input .selectize-input, .ri-mini-input .form-select { min-height: 30px; padding: 4px 10px; font-size: 0.82rem; }
.form-check-input:checked { background-color: var(--accent); border-color: var(--accent); }

/* ---------------- formulas ---------------- */
.ri-formula-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 10px; }
@media (max-width: 900px) { .ri-formula-grid { grid-template-columns: 1fr; } }
.ri-formula { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; }
.ri-formula-name { font-weight: 600; margin-bottom: 4px; }
.ri-formula code { color: var(--accent); font-size: 0.82rem; }
.ri-formula-note { color: var(--muted); font-size: 0.78rem; margin-top: 4px; }

/* ---------------- DT ---------------- */
table.dataTable.ri-dt, .ri-dt { color: var(--text-2) !important; }
table.dataTable.ri-dt thead th { color: var(--text); border-bottom: 1px solid var(--border) !important; font-weight: 600; }
table.dataTable.ri-dt tbody tr { background: transparent !important; }
table.dataTable.ri-dt tbody td { border-top: 1px solid var(--border-soft) !important; font-variant-numeric: tabular-nums; }
table.dataTable.ri-dt.compact tbody tr:hover { background: var(--card-alt) !important; }
.dataTables_wrapper { color: var(--text-2); }

/* ---------------- ggplot2 studio ---------------- */
.ri-studio-sidebar { background: var(--surface) !important; border-right: 1px solid var(--border); }
.bslib-sidebar-layout { border: 1px solid var(--border) !important; border-radius: var(--radius); background: var(--page) !important; }
.bslib-sidebar-layout > .main { background: var(--page); }
.ri-layer { font-size: 0.7rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; color: var(--accent);
  margin: 10px 0 2px; }
.ri-code { background: var(--surface); color: #cfe9dc; border: 1px solid var(--border); border-radius: 10px; padding: 14px;
  font-size: 0.8rem; white-space: pre; overflow-x: auto; margin: 0; }

/* plotly modebar on dark */
.modebar-btn path { fill: var(--muted) !important; }
.js-plotly-plot .plotly .modebar { background: transparent !important; }

/* tables: card height follows the table instead of filling the row */
.ri-card .html-widget.datatables { height: auto !important; }
.ri-card:has(.datatables) { height: auto !important; }
.ri-nav-link { display: flex; align-items: center; height: 64px; }

/* =========================================================================
 * Ligand Structure Explorer
 * ======================================================================= */
.ri-page-explorer { max-width: 1760px; }
.qx-root { display: flex; flex-direction: column; gap: 14px; }

.qx-topbar { display: flex; justify-content: space-between; align-items: center; gap: 16px; flex-wrap: wrap; }
.qx-title { display: flex; align-items: center; gap: 12px; }
.qx-kicker { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.1em; color: var(--accent); font-weight: 600; }
.qx-name { font-size: 1.35rem; font-weight: 600; color: var(--text); }
.qx-controls { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
.qx-controls .shiny-input-container { margin-bottom: 0; }

/* stage: toolbar | viewer | side panels */
.qx-stage { display: grid; grid-template-columns: 170px minmax(0, 1fr) 640px; gap: 14px; padding: 14px;
  border-radius: 20px; border: 1px solid var(--border);
  background:
    radial-gradient(ellipse 50% 70% at 35% 50%, rgba(31, 107, 76, 0.22) 0%, rgba(11, 15, 14, 0) 70%),
    radial-gradient(circle at 25% 20%, rgba(57, 135, 229, 0.10) 0%, rgba(11, 15, 14, 0) 45%),
    linear-gradient(160deg, #0f1614 0%, #0a0e0d 100%); }
@media (max-width: 1400px) { .qx-stage { grid-template-columns: 170px minmax(0, 1fr); } .qx-side { grid-column: 1 / -1; } }
@media (max-width: 800px) { .qx-stage { grid-template-columns: 1fr; } }

.qx-toolbar { background: rgba(18, 24, 22, 0.82); border: 1px solid var(--border); border-radius: 16px; padding: 14px 12px;
  align-self: start; display: flex; flex-direction: column; gap: 2px; }
.qx-tool-head { font-size: 0.68rem; text-transform: uppercase; letter-spacing: 0.08em; color: var(--accent); font-weight: 600; margin-top: 6px; }
.qx-toolbar .shiny-input-container { margin-bottom: 6px; width: 100% !important; }
/* toolbar choices as segmented buttons (radio inputs hidden) */
.qx-toolbar .shiny-options-group { display: flex; gap: 0; border: 1px solid var(--border); border-radius: 999px;
  overflow: hidden; background: var(--surface); }
.qx-toolbar .radio-inline, .qx-toolbar .form-check-inline { flex: 1; margin: 0 !important; padding: 0 !important; }
.qx-toolbar .radio-inline input, .qx-toolbar .form-check-inline input { position: absolute; opacity: 0; pointer-events: none; }
.qx-toolbar .radio-inline span, .qx-toolbar .form-check-inline .form-check-label { display: block; text-align: center;
  padding: 5px 4px; font-size: 0.8rem; color: var(--text-2); cursor: pointer; margin: 0; }
.qx-toolbar .radio-inline input:checked + span, .qx-toolbar .form-check-inline input:checked + .form-check-label {
  background: var(--accent); color: #06231a; font-weight: 700; }
.qx-tool-sep { height: 1px; background: var(--border); margin: 8px 0; }
.btn.qx-btn { background: var(--card-alt); color: var(--text); border: 1px solid var(--border); font-size: 0.8rem;
  padding: 6px 10px; margin-top: 6px; width: 100%; }
.btn.qx-btn:hover { border-color: var(--accent); color: var(--text); }
.btn.qx-btn-small { width: auto; padding: 3px 12px; margin-top: 6px; }

.qx-viewer { position: relative; min-width: 0; display: flex; flex-direction: column; gap: 8px; }
.qx-canvas { min-height: 600px; }
.qx-ident { max-width: 640px; }
.qx-ident-label { font-size: 0.7rem; color: var(--accent); font-weight: 600; letter-spacing: 0.05em; margin-top: 4px; }
.qx-ident .shiny-text-output { font-family: var(--mono); font-size: 0.78rem; color: #c3cdc8; word-break: break-all; }
.qx-legend { display: flex; flex-wrap: wrap; gap: 6px 14px; align-self: flex-start; background: rgba(18, 24, 22, 0.75);
  border: 1px solid var(--border); padding: 8px 12px; border-radius: 10px; font-size: 0.78rem; color: var(--text-2); }
.qx-legend-title { width: 100%; color: var(--text); font-weight: 600; font-size: 0.8rem; }
.qx-legend-item { display: inline-flex; align-items: center; gap: 6px; }
.qx-legend-item i { width: 11px; height: 11px; border-radius: 50%; display: inline-block; border: 1px solid rgba(255,255,255,0.15); }

.qx-side { display: flex; flex-direction: column; gap: 14px; min-width: 0; }
.qx-panel { background: rgba(18, 24, 22, 0.86); border: 1px solid var(--border); border-radius: 16px; padding: 0 10px 10px; }
.qx-panel-head { font-weight: 600; color: var(--text); font-size: 0.92rem; padding: 12px 4px 6px; display: flex;
  justify-content: space-between; align-items: baseline; gap: 8px; }
.qx-hint { color: var(--muted); font-weight: 400; font-size: 0.72rem; }
.qx-rank table.dataTable tbody tr { cursor: pointer; }
table.dataTable.qx-rank tbody tr.selected > * { box-shadow: inset 0 0 0 9999px rgba(62, 207, 142, 0.16) !important; color: var(--text) !important; }

.qx-dock .shiny-options-group { display: grid; grid-template-columns: 1fr 1fr; gap: 0 14px; }
@media (max-width: 700px) { .qx-dock .shiny-options-group { grid-template-columns: 1fr; } }
.qx-dock .shiny-input-container { width: 100% !important; margin-bottom: 0; }
.qx-dock .radio, .qx-dock .form-check { margin: 0; padding: 0; }
.qx-dock .radio input, .qx-dock .form-check-input { position: absolute; opacity: 0; pointer-events: none; }
.qx-dock .radio label, .qx-dock .form-check-label { display: block; width: 100%; padding: 0; cursor: pointer; }
.qx-desc-row { display: grid; grid-template-columns: 50px 1fr; gap: 8px; padding: 4px 6px; border-radius: 6px;
  border-bottom: 1px solid var(--border-soft); font-size: 0.78rem; line-height: 1.25; color: var(--text-2); }
.qx-desc-row:hover { background: var(--card-alt); color: var(--text); }
.qx-dock input:checked + span .qx-desc-row, .qx-dock input:checked ~ .qx-desc-row,
.qx-dock .form-check-input:checked + .form-check-label .qx-desc-row {
  background: rgba(62, 207, 142, 0.14); color: var(--text); box-shadow: inset 2px 0 0 var(--accent); }
.qx-desc-val { font-family: var(--mono); color: var(--accent); font-variant-numeric: tabular-nums; text-align: right; }
.qx-desc-none { grid-column: 1 / -1; }
.qx-desc-none .qx-desc-val { color: var(--muted); }
.qx-desc-eng .qx-desc-name { color: var(--muted); }

/* analysis strip */
.qx-strip { display: grid; grid-template-columns: 1.4fr 1fr; gap: 14px; }
.qx-strip-2 { grid-template-columns: 1fr 1.5fr; }
@media (max-width: 1250px) { .qx-strip, .qx-strip-2 { grid-template-columns: 1fr; } }
.qx-card { background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); padding: 14px 16px; min-width: 0; }
.qx-card-head { font-weight: 600; color: var(--text); margin-bottom: 8px; }
.qx-card-head-row { display: flex; justify-content: space-between; align-items: center; gap: 8px; }
.qx-card-foot { color: var(--text-2); font-size: 0.8rem; margin-top: 6px; }

.qx-calc-top { display: flex; align-items: center; gap: 18px; margin-bottom: 8px; }
.qx-calc-value { font-size: 2.6rem; font-weight: 600; color: var(--accent); line-height: 1; font-variant-numeric: tabular-nums; }
.qx-unit { font-size: 1rem; color: var(--muted); margin-left: 6px; font-weight: 500; }
.qx-calc-name { font-size: 1.05rem; font-weight: 600; color: var(--text); }
.qx-rdkit { font-size: 0.75rem; color: var(--muted); }
.qx-calc-what { color: var(--text-2); font-size: 0.85rem; margin: 6px 0 8px; }
.qx-formula { font-family: var(--mono); font-size: 0.82rem; color: #cfe9dc; background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 7px 10px; margin-bottom: 8px; }
.qx-formula.sub { color: var(--text); }
.qx-table-wrap { max-height: 230px; overflow-y: auto; border: 1px solid var(--border-soft); border-radius: 8px; }
.qx-table { width: 100%; border-collapse: collapse; font-size: 0.8rem; }
.qx-table th { position: sticky; top: 0; background: var(--card-alt); color: var(--text); font-weight: 600; text-align: left; padding: 6px 10px; }
.qx-table td { padding: 5px 10px; border-top: 1px solid var(--border-soft); color: var(--text-2); font-variant-numeric: tabular-nums; }
.qx-none { color: var(--muted); font-style: italic; }
.qx-check { margin-top: 8px; font-size: 0.8rem; }
.qx-check.ok { color: var(--accent); }
.qx-check.warn { color: #fab219; }
.qx-pill { display: inline-block; font-size: 0.72rem; font-weight: 600; padding: 2px 8px; border-radius: 999px; margin-right: 4px; }
.qx-pill.ok { background: rgba(12, 163, 12, 0.18); color: #5fdc6a; }
.qx-pill.miss { background: rgba(208, 59, 59, 0.18); color: #f08a8a; }
.qx-targets { max-height: 260px; overflow-y: auto; border: 1px solid var(--border-soft); border-radius: 8px; }
.qx-fam { font-weight: 600; color: var(--text); white-space: nowrap; }

.qx-similar { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 10px; }
@media (max-width: 700px) { .qx-similar { grid-template-columns: 1fr; } }
.qx-sim-card { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 8px 10px 10px; min-width: 0; }
.qx-sim-name { color: var(--accent); font-weight: 600; font-size: 0.85rem; margin-top: 4px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.qx-sim-meta { color: var(--muted); font-size: 0.72rem; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.qx-sim-row { display: flex; justify-content: space-between; font-size: 0.72rem; font-family: var(--mono); margin-top: 4px; color: var(--text-2); }

]---"

#' The dashboard style sheet as a <style> tag.
#'
#' @return htmltools tag.
ri_style_tag <- function() {
  htmltools::tags$style(htmltools::HTML(RI_STYLES))
}
