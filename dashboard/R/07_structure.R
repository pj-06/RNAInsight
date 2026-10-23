############################################################
# RNAInsight Dashboard
# 07_structure.R : ligand structures, descriptor breakdowns,
#                  2D (ggplot2) and 3D (r3dmol) drawings
############################################################
#
# Everything the Ligand Structure Explorer draws is computed
# here, in R:
#   * prepare_molecule()   - turns one record of
#                            data/molecules.json into R vectors
#   * feature_highlight()  - which atoms / bonds / rings a
#                            descriptor counts, and their colours
#   * feature_calculation()- the step-by-step calculation table
#   * plot_structure_2d()  - 2D structure drawn with ggplot2
#   * structure_3d()       - interactive 3D viewer via r3dmol
#
# data/molecules.json is written once by
# scripts/06_structure_export.py (RDKit); every breakdown in it
# reproduces the descriptor values stored in
# data/processed/cleaned_dataset.csv.
############################################################

# ----------------------------------------------------------
# Colours
# ----------------------------------------------------------

#' Element colours (dark-background CPK variant).
ELEMENT_COLORS <- c(
  C = "#1fa59b", H = "#e6eeec", N = "#3d8bff", O = "#e0479e", S = "#e7c84a",
  F = "#86e08f", Cl = "#4ccf6e", Br = "#c0603f", I = "#a066c8", P = "#f08c3a",
  B = "#f2b8a0", Se = "#f0a030"
)

#' Highlight colours used by the descriptor views.
HL <- list(
  dim = "#2b3532", dim_h = "#323d3a", hba = "#ff8a3d", hbd = "#4fa3ff",
  rot = "#f5d76e", aromatic = "#b18cff", aliphatic = "#f5a742",
  sp3 = "#3ecf8e", sp2 = "#5c6b66", neg = "#3987e5", mid = "#55605c",
  pos = "#e66767", mass = "#3987e5", tpsa = "#e0479e"
)

#' Colour of an element symbol.
#'
#' @param el Character vector of element symbols.
#' @return Character vector of hex colours.
element_color <- function(el) {
  out <- unname(ELEMENT_COLORS[el])
  ifelse(is.na(out), "#c9a0ff", out)
}

#' Linear mix of two hex colours.
#'
#' @param a,b Hex colours.
#' @param t Mixing weight(s) in [0, 1] (0 = a, 1 = b).
#' @return Hex colour(s).
mix_color <- function(a, b, t) {
  if (!length(a)) return(character(0))
  n <- max(length(a), length(b), length(t))
  ca <- grDevices::col2rgb(rep_len(a, n))
  cb <- grDevices::col2rgb(rep_len(b, n))
  t <- pmin(pmax(rep_len(t, n), 0), 1)
  mixed <- ca + (cb - ca) * rep(t, each = 3)
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ], maxColorValue = 255)
}

#' Sequential ramp: 0 = dark, 0.5 = hue, 1 = near white.
#'
#' @param t Value(s) in [0, 1].
#' @param hue Hex colour of the middle of the ramp.
#' @return Hex colour(s).
seq_ramp <- function(t, hue) {
  t <- pmin(pmax(t, 0), 1)
  ifelse(t < 0.5, mix_color("#24403a", hue, t * 2), mix_color(hue, "#e8fff4", (t - 0.5) * 1.1))
}

#' Diverging ramp: -1 = blue, 0 = grey, +1 = red.
#'
#' @param t Value(s) in [-1, 1].
#' @return Hex colour(s).
div_ramp <- function(t) {
  t <- pmin(pmax(t, -1), 1)
  ifelse(t < 0, mix_color(HL$mid, HL$neg, -t), mix_color(HL$mid, HL$pos, t))
}

# ----------------------------------------------------------
# Molecule records
# ----------------------------------------------------------

#' Convert one parsed JSON molecule into plain R vectors.
#'
#' Atom indices stay 0-based (as in RDKit and the SDF block) so
#' they can be passed straight to the 3D viewer; add 1 to index
#' R vectors.
#'
#' @param m One element of the parsed molecules.json list.
#' @return A list describing the molecule.
prepare_molecule <- function(m) {
  atoms <- m$atoms
  b <- m$breakdown
  num <- function(x) as.numeric(unlist(x))
  int <- function(x) as.integer(unlist(x))
  bonds <- do.call(rbind, lapply(m$bonds, function(x) as.numeric(unlist(x))))
  el <- vapply(atoms, function(a) a$el, character(1))
  # Hydrogens are explicit atoms in the export, so count them from the
  # bonds: each H is bonded to exactly one heavy atom.
  n_heavy <- as.integer(m$n_heavy)
  h_parent <- ifelse(bonds[, 1] >= n_heavy, bonds[, 2], ifelse(bonds[, 2] >= n_heavy, bonds[, 1], NA))
  h_count <- tabulate(h_parent[!is.na(h_parent)] + 1, nbins = length(el))
  h_count[seq_along(el) > n_heavy] <- 0L
  list(
    key = m$key,
    name = m$name,
    smiles = m$smiles,
    formula = m$formula,
    n_heavy = n_heavy,
    el = el,
    h = as.integer(h_count),
    chg = vapply(atoms, function(a) as.integer(a$chg), integer(1)),
    ar = vapply(atoms, function(a) isTRUE(a$ar), logical(1)),
    hyb = vapply(atoms, function(a) a$hyb, character(1)),
    bonds = data.frame(a = as.integer(bonds[, 1]), b = as.integer(bonds[, 2]), order = bonds[, 3]),
    xyz = do.call(rbind, lapply(m$xyz, num)),
    xy = do.call(rbind, lapply(m$xy, num)),
    molblock = m$molblock,
    descriptors = vapply(m$descriptors, as.numeric, numeric(1)),
    breakdown = list(
      mass_avg = num(b$mass_avg),
      mass_mono = num(b$mass_mono),
      charge_corr = as.numeric(b$charge_corr),
      elements = do.call(rbind, lapply(b$elements, function(e) {
        data.frame(el = e$el, n = as.integer(e$n), w = as.numeric(e$w), m = as.numeric(e$m))
      })),
      logp = num(b$logp),
      mr = num(b$mr),
      tpsa = num(b$tpsa),
      hba = int(b$hba),
      hba_why = as.character(unlist(b$hba_why)),
      hbd = int(b$hbd),
      hbd_why = as.character(unlist(b$hbd_why)),
      rot = lapply(b$rot, int),
      rings = lapply(b$rings, function(r) list(atoms = int(r$atoms), size = as.integer(r$size),
                                               aromatic = isTRUE(r$aromatic))),
      hetero = int(b$hetero),
      carbons = int(b$carbons),
      sp3c = int(b$sp3c)
    )
  )
}

#' Atom name used in tables, e.g. "N12" (1-based like chemistry drawings).
#'
#' @param mol Prepared molecule.
#' @param i 0-based atom index (vector).
#' @return Character vector.
atom_name <- function(mol, i) {
  paste0(mol$el[i + 1], i + 1)
}

#' Atom name with attached hydrogens / charge, e.g. "N12 (NH2)".
#'
#' @param mol Prepared molecule.
#' @param i 0-based atom index (vector).
#' @return Character vector.
atom_name_h <- function(mol, i) {
  vapply(i, function(k) {
    h <- mol$h[k + 1]
    chg <- mol$chg[k + 1]
    htxt <- if (h > 0) paste0("H", if (h > 1) h else "") else ""
    ctxt <- if (chg != 0) paste0(if (chg > 0) "+" else "−", if (abs(chg) > 1) abs(chg) else "") else ""
    base <- atom_name(mol, k)
    if (nzchar(htxt) || nzchar(ctxt)) paste0(base, " (", mol$el[k + 1], htxt, ctxt, ")") else base
  }, character(1))
}

#' Parent heavy atom (0-based) of every hydrogen (0-based), by index.
#'
#' @param mol Prepared molecule.
#' @return Integer vector named by hydrogen index.
hydrogen_parents <- function(mol) {
  n <- mol$n_heavy
  b <- mol$bonds
  from_a <- b$a >= n
  from_b <- b$b >= n
  stats::setNames(c(b$b[from_a], b$a[from_b]), c(b$a[from_a], b$b[from_b]))
}

# ----------------------------------------------------------
# Descriptor metadata and highlights
# ----------------------------------------------------------

#' What each of the 14 descriptors means and how RDKit computes it.
FEATURE_INFO <- list(
  MolecularWeight = list(
    label = "Molecular weight", unit = "g/mol", digits = 3, rdkit = "Descriptors.MolWt(mol)",
    what = "Average molecular mass: the standard atomic weight of every atom, implicit hydrogens included, added up.",
    formula = "MW = Σ (number of atoms of element × atomic weight)"),
  ExactMolWt = list(
    label = "Exact mass", unit = "Da", digits = 3, rdkit = "Descriptors.ExactMolWt(mol)",
    what = paste("Monoisotopic mass: uses the mass of each element's most abundant isotope (12C = 12.00000,",
                 "1H = 1.00783), which is what a mass spectrometer measures. It is almost perfectly correlated",
                 "with MW, so feature selection dropped it."),
    formula = "Exact mass = Σ (count × monoisotopic mass) − charge × m(electron)"),
  LogP = list(
    label = "LogP", unit = "", digits = 3, rdkit = "Descriptors.MolLogP(mol)  (Wildman-Crippen)",
    what = paste("Octanol/water partition coefficient. Each atom (hydrogens included) is assigned one of about",
                 "70 Crippen atom types, each with a fixed contribution, and LogP is the sum. Positive means",
                 "lipophilic (red), negative means hydrophilic (blue). Hydrogen contributions are added to",
                 "their parent atom here."),
    formula = "LogP = Σ atom-type contributions"),
  MolMR = list(
    label = "Molar refractivity", unit = "", digits = 3, rdkit = "Descriptors.MolMR(mol)  (Wildman-Crippen)",
    what = paste("Molar refractivity: a measure of molecular volume and polarisability. It uses the same atom",
                 "types as LogP, but each type has a refractivity contribution, so big polarisable atoms",
                 "(S, Br, aromatic C) add the most."),
    formula = "MR = Σ atom-type refractivity contributions"),
  TPSA = list(
    label = "TPSA", unit = "Å²", digits = 2, rdkit = "Descriptors.TPSA(mol)  (Ertl 2000)",
    what = paste("Topological polar surface area: the summed surface of polar atoms (N and O with their",
                 "attached H), read from tabulated fragment values for each N/O environment, so no 3D",
                 "structure is needed. S and P are excluded, as in RDKit's default."),
    formula = "TPSA = Σ fragment contributions of polar N and O atoms"),
  HBA = list(
    label = "H-bond acceptors", unit = "", digits = 0, rdkit = "Descriptors.NumHAcceptors(mol)",
    what = paste("Atoms that can accept a hydrogen bond, matched with RDKit's SMARTS rules. Amide",
                 "nitrogens, and O/N whose lone pair is tied up next to a carbonyl, are not counted."),
    formula = "HBA = number of atoms matching the acceptor pattern"),
  HBD = list(
    label = "H-bond donors", unit = "", digits = 0, rdkit = "Descriptors.NumHDonors(mol)",
    what = paste("Heteroatoms carrying a hydrogen that can be donated in a hydrogen bond (N-H, O-H, S-H).",
                 "Each atom counts once, however many H it has."),
    formula = "HBD = number of atoms matching the donor pattern"),
  RotatableBonds = list(
    label = "Rotatable bonds", unit = "", digits = 0, rdkit = "Descriptors.NumRotatableBonds(mol)  (Strict)",
    what = paste("Single, non-ring bonds between two non-terminal heavy atoms, the bonds the molecule can",
                 "twist around. RDKit's strict rules also skip amide-type C(=O)-N bonds, triple bonds, and",
                 "bonds to CH3, CF3 or t-butyl groups."),
    formula = "RotB = number of bonds matching the rotatable-bond pattern"),
  RingCount = list(
    label = "Rings", unit = "", digits = 0, rdkit = "Descriptors.RingCount(mol)",
    what = "Number of rings in the smallest set of smallest rings (SSSR). Fused systems such as naphthalene count as 2.",
    formula = "Rings = |SSSR|"),
  NumAromaticRings = list(
    label = "Aromatic rings", unit = "", digits = 0, rdkit = "Descriptors.NumAromaticRings(mol)",
    what = paste("SSSR rings in which every bond is aromatic. Flat aromatic rings can stack between RNA",
                 "bases, which matters for intercalators and G-quadruplex binders."),
    formula = "Aromatic rings = rings whose bonds are all aromatic"),
  NumAliphaticRings = list(
    label = "Aliphatic rings", unit = "", digits = 0, rdkit = "Descriptors.NumAliphaticRings(mol)",
    what = "SSSR rings with at least one non-aromatic bond (saturated or partly saturated rings, such as piperazine or a sugar ring).",
    formula = "Aliphatic rings = rings with ≥ 1 non-aromatic bond"),
  HeavyAtomCount = list(
    label = "Heavy atoms", unit = "", digits = 0, rdkit = "Descriptors.HeavyAtomCount(mol)",
    what = "Every atom that is not hydrogen. Used as the size denominator in two engineered features.",
    formula = "Heavy atoms = all non-H atoms"),
  NumHeteroAtoms = list(
    label = "Heteroatoms", unit = "", digits = 0, rdkit = "Descriptors.NumHeteroatoms(mol)",
    what = paste("Atoms other than carbon and hydrogen (N, O, S, halogens, P ...). They carry most of the",
                 "molecule's polarity and its ability to H-bond with RNA."),
    formula = "Heteroatoms = atoms that are not C or H"),
  FractionCSP3 = list(
    label = "Fraction sp³ C", unit = "", digits = 3, rdkit = "Descriptors.FractionCSP3(mol)",
    what = paste("The share of carbons that are sp3 (four single bonds, tetrahedral). High values mean a",
                 "saturated, three-dimensional molecule; low values mean a flat, aromatic one."),
    formula = "Fsp3 = sp3 carbons ÷ all carbons")
)

#' Engineered features (R/04_feature_engineering.R) and the descriptor
#' whose atoms they highlight.
ENGINEERED_INFO <- list(
  HBondTotal = list(label = "H-bond total", formula = "HBA + HBD", show = c("HBA", "HBD"),
                    calc = function(d) d[["HBA"]] + d[["HBD"]]),
  LipinskiViolations = list(label = "Lipinski violations", formula = "[MW > 500] + [LogP > 5] + [HBA > 10] + [HBD > 5]",
                            show = character(0),
                            calc = function(d) (d[["MolecularWeight"]] > 500) + (d[["LogP"]] > 5) +
                              (d[["HBA"]] > 10) + (d[["HBD"]] > 5)),
  AromaticRingFraction = list(label = "Aromatic ring fraction", formula = "NumAromaticRings / RingCount",
                              show = "RingCount",
                              calc = function(d) if (d[["RingCount"]] == 0) 0 else d[["NumAromaticRings"]] / d[["RingCount"]]),
  FlexibilityIndex = list(label = "Flexibility index", formula = "RotatableBonds / HeavyAtomCount",
                          show = "RotatableBonds", calc = function(d) d[["RotatableBonds"]] / d[["HeavyAtomCount"]]),
  HeteroatomFraction = list(label = "Heteroatom fraction", formula = "NumHeteroAtoms / HeavyAtomCount",
                            show = "NumHeteroAtoms", calc = function(d) d[["NumHeteroAtoms"]] / d[["HeavyAtomCount"]]),
  TPSA_per_MW = list(label = "TPSA per MW", formula = "TPSA / MolecularWeight",
                     show = "TPSA", calc = function(d) d[["TPSA"]] / d[["MolecularWeight"]])
)

#' Empty highlight (plain element colours).
empty_highlight <- function() {
  list(atoms = character(0), glow = numeric(0), values = numeric(0),
       bonds = data.frame(a = integer(0), b = integer(0), color = character(0)),
       rings = list(), legend = data.frame(color = character(0), label = character(0)),
       h_color = FALSE, soft = FALSE, keep_element = FALSE)
}

#' Highlight driven by one value per heavy atom (mass, LogP, MR, TPSA).
#'
#' @param values Numeric vector, one value per heavy atom.
#' @param color_fn Function(value, max_abs) returning colours.
#' @param legend Legend data frame (color, label).
#' @param skip_zero Do not highlight atoms with a zero value.
#' @return Highlight list.
value_highlight <- function(values, color_fn, legend, skip_zero = FALSE) {
  h <- empty_highlight()
  idx <- seq_along(values) - 1L
  keep <- if (skip_zero) abs(values) > 1e-9 else rep(TRUE, length(values))
  max_abs <- max(1e-9, abs(values))
  h$atoms <- stats::setNames(color_fn(values[keep], max_abs), idx[keep])
  h$glow <- stats::setNames(0.35 + 0.65 * abs(values[keep]) / max_abs, idx[keep])
  h$values <- stats::setNames(values[keep], idx[keep])
  h$legend <- legend
  h$h_color <- TRUE
  h
}

#' Highlight a set of atoms in one colour.
#'
#' @param idx 0-based atom indices.
#' @param color Hex colour.
#' @param legend_label Legend text.
#' @return Highlight list.
set_highlight <- function(idx, color, legend_label) {
  h <- empty_highlight()
  h$atoms <- stats::setNames(rep(color, length(idx)), idx)
  h$glow <- stats::setNames(rep(1, length(idx)), idx)
  h$legend <- data.frame(color = color, label = legend_label)
  h
}

#' Highlight rings (all, aromatic only or aliphatic only).
#'
#' @param mol Prepared molecule.
#' @param filter NULL, "aromatic" or "aliphatic".
#' @return Highlight list.
ring_highlight <- function(mol, filter = NULL) {
  h <- empty_highlight()
  for (k in seq_along(mol$breakdown$rings)) {
    r <- mol$breakdown$rings[[k]]
    if (identical(filter, "aromatic") && !r$aromatic) next
    if (identical(filter, "aliphatic") && r$aromatic) next
    col <- if (r$aromatic) HL$aromatic else HL$aliphatic
    h$rings[[length(h$rings) + 1]] <- list(atoms = r$atoms, color = col, label = paste0("R", k))
    h$atoms[as.character(r$atoms)] <- col
  }
  h$legend <- data.frame(color = c(HL$aromatic, HL$aliphatic), label = c("aromatic ring", "aliphatic ring"))
  h
}

#' Which atoms, bonds or rings a descriptor counts, and their colours.
#'
#' @param mol Prepared molecule.
#' @param feature Descriptor name, engineered feature name or "none".
#' @return Highlight list (atoms named by 0-based index).
feature_highlight <- function(mol, feature) {
  b <- mol$breakdown
  n <- mol$n_heavy
  if (is.null(feature) || feature == "none") return(empty_highlight())

  if (feature %in% names(ENGINEERED_INFO)) {
    parts <- ENGINEERED_INFO[[feature]]$show
    if (!length(parts)) return(empty_highlight())
    hs <- lapply(parts, function(p) feature_highlight(mol, p))
    h <- hs[[1]]
    if (length(hs) > 1) {
      for (extra in hs[-1]) {
        h$atoms[names(extra$atoms)] <- extra$atoms
        h$glow[names(extra$glow)] <- extra$glow
        h$legend <- rbind(h$legend, extra$legend)
      }
      h$h_color <- TRUE
    }
    return(h)
  }

  switch(
    feature,
    MolecularWeight = value_highlight(b$mass_avg, function(v, m) seq_ramp(v / m, HL$mass),
                                      data.frame(color = c(seq_ramp(0.25, HL$mass), seq_ramp(1, HL$mass)),
                                                 label = c("light atom", "heavy atom (incl. its H)"))),
    ExactMolWt = value_highlight(b$mass_mono, function(v, m) seq_ramp(v / m, HL$mass),
                                 data.frame(color = c(seq_ramp(0.25, HL$mass), seq_ramp(1, HL$mass)),
                                            label = c("light atom", "heavy atom (incl. its H)"))),
    LogP = value_highlight(b$logp, function(v, m) div_ramp(v / m),
                           data.frame(color = c(HL$neg, HL$mid, HL$pos),
                                      label = c("hydrophilic (−)", "≈ 0", "lipophilic (+)"))),
    MolMR = value_highlight(b$mr, function(v, m) seq_ramp(v / m, HL$mass),
                            data.frame(color = c(seq_ramp(0.25, HL$mass), seq_ramp(1, HL$mass)),
                                       label = c("small contribution", "large contribution"))),
    TPSA = value_highlight(b$tpsa, function(v, m) seq_ramp(v / m, HL$tpsa),
                           data.frame(color = c(seq_ramp(0.3, HL$tpsa), seq_ramp(1, HL$tpsa)),
                                      label = c("small polar contribution", "large polar contribution")),
                           skip_zero = TRUE),
    HBA = set_highlight(b$hba, HL$hba, "acceptor atom"),
    HBD = {
      h <- set_highlight(b$hbd, HL$hbd, "donor atom (and its H)")
      h$h_color <- TRUE
      h
    },
    RotatableBonds = {
      h <- empty_highlight()
      if (length(b$rot)) {
        pairs <- do.call(rbind, b$rot)
        h$bonds <- data.frame(a = pairs[, 1], b = pairs[, 2], color = HL$rot)
        h$atoms[as.character(unique(c(pairs)))] <- HL$rot
      }
      h$legend <- data.frame(color = HL$rot, label = "rotatable bond")
      h$soft <- TRUE
      h
    },
    RingCount = ring_highlight(mol),
    NumAromaticRings = ring_highlight(mol, "aromatic"),
    NumAliphaticRings = ring_highlight(mol, "aliphatic"),
    HeavyAtomCount = {
      h <- empty_highlight()
      h$atoms <- stats::setNames(element_color(mol$el[seq_len(n)]), seq_len(n) - 1L)
      h$legend <- data.frame(color = HL$dim_h, label = "hydrogen (not counted)")
      h$keep_element <- TRUE
      h
    },
    NumHeteroAtoms = {
      h <- empty_highlight()
      h$atoms <- stats::setNames(element_color(mol$el[b$hetero + 1]), b$hetero)
      h$glow <- stats::setNames(rep(1, length(b$hetero)), b$hetero)
      h$legend <- data.frame(color = element_color(c("N", "O", "S", "F")), label = c("N", "O", "S", "halogen"))
      h
    },
    FractionCSP3 = {
      h <- empty_highlight()
      h$atoms <- stats::setNames(rep(HL$sp2, length(b$carbons)), b$carbons)
      h$atoms[as.character(b$sp3c)] <- HL$sp3
      h$glow <- stats::setNames(rep(1, length(b$sp3c)), b$sp3c)
      h$legend <- data.frame(color = c(HL$sp3, HL$sp2), label = c("sp³ carbon", "sp² / sp carbon"))
      h
    },
    empty_highlight()
  )
}

#' Final display colour of every atom (heavy atoms and hydrogens).
#'
#' @param mol Prepared molecule.
#' @param h Highlight list.
#' @return Character vector, one colour per atom (R order = 0-based index + 1).
atom_display_colors <- function(mol, h) {
  n_all <- length(mol$el)
  n <- mol$n_heavy
  has <- length(h$atoms) > 0
  cols <- character(n_all)
  for (i in seq_len(n)) {
    key <- as.character(i - 1)
    cols[i] <- if (!has) element_color(mol$el[i]) else if (!is.na(h$atoms[key])) h$atoms[[key]] else
      if (isTRUE(h$soft)) mix_color(element_color(mol$el[i]), HL$dim, 0.65) else HL$dim
  }
  if (n_all > n) {
    parents <- hydrogen_parents(mol)
    for (i in (n + 1):n_all) {
      parent <- parents[as.character(i - 1)]
      pkey <- as.character(parent)
      cols[i] <- if (has && isTRUE(h$h_color) && !is.na(parent) && !is.na(h$atoms[pkey])) {
        mix_color(h$atoms[[pkey]], "#ffffff", 0.35)
      } else if (has) HL$dim_h else ELEMENT_COLORS[["H"]]
    }
  }
  cols
}

# ----------------------------------------------------------
# Step-by-step calculation
# ----------------------------------------------------------

#' Largest per-atom contributions as a table.
contribution_table <- function(mol, values, n_rows, digits, header) {
  idx <- which(abs(values) > 1e-9)
  idx <- idx[order(-abs(values[idx]))]
  idx <- utils::head(idx, n_rows)
  out <- data.frame(atom_name_h(mol, idx - 1), sprintf(paste0("%+.", digits, "f"), values[idx]))
  names(out) <- header
  out
}

#' Rings as a table (optionally aromatic or aliphatic only).
ring_table <- function(mol, filter = NULL) {
  rows <- list()
  for (k in seq_along(mol$breakdown$rings)) {
    r <- mol$breakdown$rings[[k]]
    if (identical(filter, "aromatic") && !r$aromatic) next
    if (identical(filter, "aliphatic") && r$aromatic) next
    tab <- table(mol$el[r$atoms + 1])
    comp <- paste0(names(tab), ifelse(tab > 1, tab, ""), collapse = "")
    rows[[length(rows) + 1]] <- data.frame(Ring = paste0("R", k), Size = paste0(r$size, "-membered"),
                                           Atoms = comp, Type = if (r$aromatic) "aromatic" else "aliphatic")
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(Ring = character(0), Size = character(0), Atoms = character(0), Type = character(0))
}

#' Element table (count x mass = subtotal).
element_table <- function(mol, column, header) {
  e <- mol$breakdown$elements
  out <- data.frame(e$el, e$n, sprintf("%.5f", e[[column]]), sprintf("%.3f", e$n * e[[column]]))
  names(out) <- c("Element", "Count", header, "Subtotal")
  out
}

#' How a descriptor is calculated for one ligand.
#'
#' @param mol Prepared molecule.
#' @param feature Descriptor or engineered-feature name.
#' @return list(label, value, unit, rdkit, what, formula, table, total, ok, engineered)
feature_calculation <- function(mol, feature) {
  b <- mol$breakdown
  d <- mol$descriptors

  if (feature %in% names(ENGINEERED_INFO)) {
    e <- ENGINEERED_INFO[[feature]]
    terms <- e$formula
    for (nm in names(d)) {
      terms <- gsub(paste0("\\b", nm, "\\b"), format(signif(d[[nm]], 5)), terms)
    }
    return(list(label = e$label, value = e$calc(d), unit = "", rdkit = "R/04_feature_engineering.R",
                what = paste("Engineered in R from the RDKit descriptors, so no new chemistry is computed.",
                             "The structure highlights the descriptor(s) it is built from."),
                formula = paste(feature, "=", e$formula), worked = paste("=", terms),
                table = NULL, total = e$calc(d), ok = TRUE, engineered = TRUE, digits = 3))
  }

  info <- FEATURE_INFO[[feature]]
  table <- switch(
    feature,
    MolecularWeight = element_table(mol, "w", "Atomic weight"),
    ExactMolWt = {
      t <- element_table(mol, "m", "Isotope mass")
      if (abs(b$charge_corr) > 0) {
        t <- rbind(t, stats::setNames(data.frame("e⁻", "", "charge correction", sprintf("%.5f", b$charge_corr)), names(t)))
      }
      t
    },
    LogP = {
      t <- contribution_table(mol, b$logp, 8, 3, c("Largest contributions", "Value"))
      rbind(t, stats::setNames(data.frame(c("... all positive atoms", "... all negative atoms"),
                                          sprintf("%+.3f", c(sum(b$logp[b$logp > 0]), sum(b$logp[b$logp < 0])))),
                               names(t)))
    },
    MolMR = contribution_table(mol, b$mr, 10, 3, c("Largest contributions", "Value")),
    TPSA = contribution_table(mol, b$tpsa, 40, 2, c("Polar atom", "Contribution (Å²)")),
    HBA = data.frame(Atom = atom_name_h(mol, b$hba), `Why it counts` = b$hba_why, check.names = FALSE),
    HBD = data.frame(Atom = atom_name_h(mol, b$hbd), `Why it counts` = b$hbd_why, check.names = FALSE),
    RotatableBonds = data.frame(
      `#` = seq_along(b$rot),
      Bond = vapply(b$rot, function(p) paste(atom_name(mol, p[1]), "–", atom_name(mol, p[2])), character(1)),
      check.names = FALSE
    ),
    RingCount = ring_table(mol),
    NumAromaticRings = ring_table(mol, "aromatic"),
    NumAliphaticRings = ring_table(mol, "aliphatic"),
    HeavyAtomCount = {
      e <- b$elements[b$elements$el != "H", ]
      data.frame(Element = e$el, Count = e$n)
    },
    NumHeteroAtoms = {
      e <- b$elements[!(b$elements$el %in% c("C", "H")), ]
      data.frame(Element = e$el, Count = e$n)
    },
    FractionCSP3 = data.frame(
      Quantity = c("sp³ carbons", "all carbons", "ratio"),
      Value = c(length(b$sp3c), length(b$carbons), paste(length(b$sp3c), "÷", length(b$carbons)))
    )
  )

  total <- switch(
    feature,
    MolecularWeight = sum(b$elements$n * b$elements$w),
    ExactMolWt = sum(b$elements$n * b$elements$m) + b$charge_corr,
    LogP = sum(b$logp), MolMR = sum(b$mr), TPSA = sum(b$tpsa),
    HBA = length(b$hba), HBD = length(b$hbd), RotatableBonds = length(b$rot),
    RingCount = length(b$rings),
    NumAromaticRings = sum(vapply(b$rings, `[[`, logical(1), "aromatic")),
    NumAliphaticRings = sum(!vapply(b$rings, `[[`, logical(1), "aromatic")),
    HeavyAtomCount = mol$n_heavy, NumHeteroAtoms = length(b$hetero),
    FractionCSP3 = if (length(b$carbons)) length(b$sp3c) / length(b$carbons) else 0
  )
  tolerance <- if (info$digits == 0) 0 else 0.0015
  list(label = info$label, value = d[[feature]], unit = info$unit, rdkit = info$rdkit, what = info$what,
       formula = info$formula, worked = NULL, table = table, total = total,
       ok = abs(total - d[[feature]]) <= tolerance + 1e-9, engineered = FALSE, digits = info$digits)
}

#' Format a descriptor value with its natural number of decimals.
#'
#' @param x Number.
#' @param digits Decimals (0 = integer).
#' @return Character.
format_value <- function(x, digits) {
  if (digits == 0) format(round(x)) else formatC(x, format = "f", digits = digits)
}

# ----------------------------------------------------------
# 2D drawing (ggplot2)
# ----------------------------------------------------------

#' Bond line segments for a 2D depiction.
#'
#' Single bonds are one line, double bonds two, triple bonds three;
#' aromatic bonds get a dashed inner line. Lines are shortened where
#' they meet a labelled (non-carbon) atom.
#'
#' @param mol Prepared molecule.
#' @param P Matrix of 2D atom positions (heavy atoms).
#' @param bond_len Typical bond length in plot units.
#' @param labelled Logical vector: atom has a text label.
#' @return Data frame of segments (x, y, xend, yend, dashed, a, b).
bond_segments <- function(mol, P, bond_len, labelled) {
  n <- mol$n_heavy
  rings <- mol$breakdown$rings
  ring_centre <- function(a, b) {
    best <- NULL
    for (r in rings) {
      at <- r$atoms
      k <- match(a, at)
      if (is.na(k)) next
      nxt <- at[(k %% length(at)) + 1]
      prv <- at[((k - 2) %% length(at)) + 1]
      if (b == nxt || b == prv) {
        if (is.null(best) || length(at) < best$size) {
          best <- list(centre = colMeans(P[at + 1, , drop = FALSE]), size = length(at))
        }
      }
    }
    best
  }
  rows <- list()
  off <- bond_len * 0.16
  shrink <- bond_len * 0.28
  for (k in seq_len(nrow(mol$bonds))) {
    a <- mol$bonds$a[k]
    b <- mol$bonds$b[k]
    order <- mol$bonds$order[k]
    if (a >= n || b >= n) next
    p1 <- P[a + 1, ]
    p2 <- P[b + 1, ]
    v <- p2 - p1
    len <- sqrt(sum(v^2))
    u <- v / len
    if (labelled[a + 1]) p1 <- p1 + u * shrink
    if (labelled[b + 1]) p2 <- p2 - u * shrink
    nrm <- c(-u[2], u[1])
    ring <- ring_centre(a, b)
    if (!is.null(ring)) {
      mid <- (P[a + 1, ] + P[b + 1, ]) / 2
      if (sum((ring$centre - mid) * nrm) < 0) nrm <- -nrm
    }
    add <- function(q1, q2, dashed = FALSE) {
      rows[[length(rows) + 1]] <<- data.frame(x = q1[1], y = q1[2], xend = q2[1], yend = q2[2],
                                              dashed = dashed, a = a, b = b)
    }
    if (order == 1) {
      add(p1, p2)
    } else if (order == 1.5 || (order == 2 && !is.null(ring))) {
      add(p1, p2)
      t <- 0.14 * len
      add(p1 + nrm * off + u * t, p2 + nrm * off - u * t, dashed = order == 1.5)
    } else if (order == 2) {
      add(p1 + nrm * off / 2, p2 + nrm * off / 2)
      add(p1 - nrm * off / 2, p2 - nrm * off / 2)
    } else {
      add(p1, p2)
      add(p1 + nrm * off, p2 + nrm * off)
      add(p1 - nrm * off, p2 - nrm * off)
    }
  }
  do.call(rbind, rows)
}

#' Draw a ligand in 2D with ggplot2, highlighting what a descriptor counts.
#'
#' @param mol Prepared molecule.
#' @param h Highlight list from feature_highlight() (or NULL).
#' @param show_names Label every heavy atom with its name (e.g. N12).
#' @param mini Compact drawing for small cards (no values, smaller text).
#' @return ggplot object.
plot_structure_2d <- function(mol, h = NULL, show_names = FALSE, mini = FALSE) {
  if (is.null(h)) h <- empty_highlight()
  n <- mol$n_heavy
  P <- mol$xy[seq_len(n), , drop = FALSE]
  bond_len <- 1.5
  labelled <- mol$el[seq_len(n)] != "C" | mol$chg[seq_len(n)] != 0
  has <- length(h$atoms) > 0

  atoms <- data.frame(idx = seq_len(n) - 1L, x = P[, 1], y = P[, 2], el = mol$el[seq_len(n)])
  atoms$label <- vapply(seq_len(n), function(i) {
    hh <- mol$h[i]
    chg <- mol$chg[i]
    paste0(mol$el[i], if (hh > 0) paste0("H", if (hh > 1) hh else "") else "",
           if (chg > 0) "+" else if (chg < 0) "−" else "")
  }, character(1))
  in_hl <- as.character(atoms$idx) %in% names(h$atoms)
  atoms$color <- element_color(atoms$el)
  if (has && !isTRUE(h$keep_element)) atoms$color[!in_hl] <- mix_color(atoms$color[!in_hl], "#0b0f0e", 0.55)

  segs <- bond_segments(mol, P, bond_len, labelled)
  segs$color <- if (!has) "#c7d0cc" else
    ifelse(as.character(segs$a) %in% names(h$atoms) & as.character(segs$b) %in% names(h$atoms), "#dfe6e3", "#55625d")

  size_scale <- if (mini) 0.55 else 1
  g <- ggplot2::ggplot()

  # rings: translucent fill + outline + label
  for (r in h$rings) {
    poly <- data.frame(x = P[r$atoms + 1, 1], y = P[r$atoms + 1, 2])
    g <- g + ggplot2::geom_polygon(data = poly, ggplot2::aes(x, y), fill = r$color, alpha = 0.18,
                                   colour = r$color, linewidth = 1.6 * size_scale)
    if (!mini) {
      g <- g + ggplot2::annotate("text", x = mean(poly$x), y = mean(poly$y), label = r$label,
                                 colour = r$color, fontface = "bold", size = 4)
    }
  }
  # highlighted bonds (rotatable bonds)
  if (nrow(h$bonds)) {
    hb <- data.frame(x = P[h$bonds$a + 1, 1], y = P[h$bonds$a + 1, 2],
                     xend = P[h$bonds$b + 1, 1], yend = P[h$bonds$b + 1, 2])
    g <- g + ggplot2::geom_segment(data = hb, ggplot2::aes(x, y, xend = xend, yend = yend),
                                   colour = HL$rot, alpha = 0.55, linewidth = 6 * size_scale, lineend = "round")
  }
  # glow discs
  if (length(h$glow)) {
    gl <- atoms[match(names(h$glow), as.character(atoms$idx)), ]
    gl$size <- (7 + 9 * h$glow) * size_scale
    gl$fill <- h$atoms[names(h$glow)]
    g <- g + ggplot2::geom_point(data = gl, ggplot2::aes(x, y, size = size, colour = fill), alpha = 0.33)
  }
  # bonds
  g <- g +
    ggplot2::geom_segment(data = segs[!segs$dashed, ], ggplot2::aes(x, y, xend = xend, yend = yend, colour = color),
                          linewidth = 0.7 * size_scale, lineend = "round") +
    ggplot2::geom_segment(data = segs[segs$dashed, ], ggplot2::aes(x, y, xend = xend, yend = yend, colour = color),
                          linewidth = 0.7 * size_scale, linetype = "22")
  # atom labels (heteroatoms, charged atoms)
  lab <- atoms[labelled, ]
  if (nrow(lab)) {
    g <- g + ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = label, colour = color),
                                fontface = "bold", size = if (mini) 2.4 else 4.2)
  }
  # small dots on highlighted carbons
  dots <- atoms[in_hl & !labelled & !(as.character(atoms$idx) %in% names(h$glow)), ]
  if (nrow(dots) && !mini) {
    dots$fill <- h$atoms[as.character(dots$idx)]
    g <- g + ggplot2::geom_point(data = dots, ggplot2::aes(x, y, colour = fill), size = 2.2)
  }
  # contribution values (when few) and atom names
  if (!mini && length(h$values)) {
    keep <- names(h$values)[abs(h$values) > 1e-9]
    if (length(keep) <= 14) {
      vl <- atoms[match(keep, as.character(atoms$idx)), ]
      vl$txt <- sprintf("%.2f", h$values[keep])
      g <- g + ggplot2::geom_text(data = vl, ggplot2::aes(x + 0.42, y + 0.42, label = txt),
                                  colour = RI_COLORS$text, size = 3.1)
    }
  }
  if (show_names && !mini) {
    g <- g + ggplot2::geom_text(data = atoms, ggplot2::aes(x - 0.32, y - 0.36, label = paste0(el, idx + 1)),
                                colour = RI_COLORS$muted, size = 2.6)
  }
  pad <- if (mini) 0.6 else 1
  g +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal(xlim = range(P[, 1]) + c(-pad, pad), ylim = range(P[, 2]) + c(-pad, pad), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   legend.position = "none")
}

# ----------------------------------------------------------
# 3D viewer (r3dmol)
# ----------------------------------------------------------

#' Atom selection by 0-based index, in the format r3dmol / 3Dmol expect.
#'
#' @param idx Integer vector of 0-based atom indices.
#' @return Selection list.
index_sel <- function(idx) {
  sel <- list(index = as.list(as.integer(idx)))
  class(sel) <- "AtomSelectionSpec"
  sel
}

#' Interactive 3D structure with the descriptor highlight (r3dmol).
#'
#' @param mol Prepared molecule.
#' @param h Highlight list.
#' @param style "ballstick", "stick" or "sphere".
#' @param show_h Show hydrogens?
#' @param show_labels Label heteroatoms?
#' @param spin Auto-rotate?
#' @return r3dmol htmlwidget.
structure_3d <- function(mol, h, style = "ballstick", show_h = TRUE, show_labels = FALSE, spin = FALSE) {
  cols <- atom_display_colors(mol, h)
  n_all <- length(mol$el)
  idx_all <- seq_len(n_all) - 1L
  visible <- if (show_h) idx_all else idx_all[seq_len(mol$n_heavy)]

  style_for <- function(col) {
    # colorScheme = NULL: r3dmol's default colour scheme ("default",
    # element colours) would otherwise override the per-atom colour.
    switch(style,
           stick = r3dmol::m_style_stick(radius = 0.17, color = col, colorScheme = NULL),
           sphere = r3dmol::m_style_sphere(scale = 0.9, color = col, colorScheme = NULL),
           c(r3dmol::m_style_stick(radius = 0.13, color = col, colorScheme = NULL),
             r3dmol::m_style_sphere(scale = 0.24, color = col, colorScheme = NULL)))
  }

  v <- r3dmol::r3dmol(viewer_spec = r3dmol::m_viewer_spec(backgroundColor = "#0f1614", antialias = TRUE)) |>
    r3dmol::m_add_model(data = mol$molblock, format = "sdf", keepH = TRUE) |>
    r3dmol::m_set_style(sel = r3dmol::m_sel(), style = r3dmol::m_style_stick(hidden = TRUE))

  # one style call per colour group
  for (col in unique(cols[visible + 1])) {
    group <- visible[cols[visible + 1] == col]
    v <- v |> r3dmol::m_set_style(sel = index_sel(group), style = style_for(col))
  }

  xyz <- function(i) r3dmol::m_vector3(mol$xyz[i + 1, 1], mol$xyz[i + 1, 2], mol$xyz[i + 1, 3])

  if (length(h$glow) && style != "sphere") {
    for (k in names(h$glow)) {
      i <- as.integer(k)
      v <- v |> r3dmol::m_add_sphere(center = xyz(i), radius = 0.45 + 0.55 * h$glow[[k]],
                                     spec = r3dmol::m_shape_spec(color = h$atoms[[k]], opacity = 0.28))
    }
  }
  if (nrow(h$bonds)) {
    v <- v |> r3dmol::m_add_cylinder(
      start = lapply(h$bonds$a, xyz), end = lapply(h$bonds$b, xyz),
      radius = 0.32, color = HL$rot, alpha = 0.45
    )
  }
  for (r in h$rings) {
    at <- r$atoms
    v <- v |> r3dmol::m_add_cylinder(
      start = lapply(at, xyz), end = lapply(c(at[-1], at[1]), xyz),
      radius = 0.24, color = r$color, alpha = 0.4
    )
    centre <- colMeans(mol$xyz[at + 1, , drop = FALSE])
    # sel = NULL: place the label at `position` (a selection would move it
    # to the centre of the selected atoms)
    v <- v |> r3dmol::m_add_label(
      text = r$label, sel = NULL,
      style = r3dmol::m_style_label(fontSize = 11, fontColor = "#0b0f0e", backgroundColor = r$color,
                                    backgroundOpacity = 0.9, alignment = "center",
                                    position = r3dmol::m_vector3(centre[1], centre[2], centre[3]))
    )
  }
  if (show_labels) {
    for (i in which(mol$el[seq_len(mol$n_heavy)] != "C") - 1L) {
      v <- v |> r3dmol::m_add_label(
        text = atom_name(mol, i), sel = NULL,
        style = r3dmol::m_style_label(fontSize = 10, fontColor = "#eef1ef", backgroundColor = "#0b0f0e",
                                      backgroundOpacity = 0.6, alignment = "center",
                                      position = xyz(i))
      )
    }
  }
  # The r3dmol viewer object is reused between redraws in Shiny, so a
  # rotation once started keeps going unless it is explicitly set again:
  # speed 0 stops it when auto-rotate is off.
  v |>
    r3dmol::m_zoom_to() |>
    r3dmol::m_spin(axis = "y", speed = if (spin) 0.6 else 0)
}
