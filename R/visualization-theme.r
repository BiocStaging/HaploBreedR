# ============================================================================
# R/visualization-theme.r
# Colour palettes and ggplot2 themes for HaploBreedR.
# ============================================================================


# ============================================================================
# 1. Palette registry -- 10 named colour sets
# ============================================================================

#' Named HaploBreedR colour palettes
#'
#' A named list of ten colour palettes. The first four are
#' colourblind-safe; the last six are journal-inspired and chosen for
#' print reproduction. Each palette is a character vector of hex colours
#' ordered so the first is the "primary" (superior / highlighted), the
#' second is the secondary, and the rest fill categorical groups.
#'
#' @format A named list of character vectors.
#'
#' @noRd
hbr_palettes <- list(
  haplobreed = c("#D55E00", "#0072B2", "#009E73", "#CC79A7",
                 "#F0E442", "#56B4E9", "#E69F00", "#999999"),
  okabe_ito  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999"),
  npg        = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85"),
  nature     = c("#2166AC", "#B2182B", "#1B7837", "#F1A340",
                 "#762A83", "#666666", "#4393C3", "#D6604D"),
  lancet     = c("#00468B", "#ED0000", "#42B540", "#0099B4",
                 "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6",
                 "#1B1919"),
  jama       = c("#374E55", "#DF8F44", "#00A1D5", "#B24745",
                 "#79AF97", "#6A6599", "#80796B"),
  nejm       = c("#BC3C29", "#0072B5", "#E18727", "#20854E",
                 "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"),
  jco        = c("#0073C2", "#EFC000", "#868686", "#CD534C",
                 "#7AA6DC", "#003C67", "#8F7700", "#3B3B3B",
                 "#A73030"),
  aaas       = c("#3B4992", "#EE0000", "#008B45", "#631879",
                 "#008280", "#BB0021", "#5F559B", "#A20056",
                 "#808180", "#1B1919"),
  frontiers  = c("#348ABD", "#A60628", "#7A68A6", "#467821",
                 "#CF4457", "#188487", "#E24A33")
)

#' Names of the built-in colour palettes
#'
#' @return A character vector of palette names accepted by
#'   [hbr_palette()].
#'
#' @examples
#' hbr_palette_names()
#'
#' @export
hbr_palette_names <- function() names(hbr_palettes)

#' Return `n` colours from a named palette
#'
#' @param name Palette name. See [hbr_palette_names()].
#' @param n Number of colours. If `NULL`, returns the whole palette. If
#'   `n > length(palette)`, the palette is interpolated.
#'
#' @return Character vector of hex colours.
#'
#' @examples
#' hbr_palette("okabe_ito", n = 4)
#' hbr_palette("haplobreed")
#'
#' @export
hbr_palette <- function(name = "haplobreed", n = NULL) {
  if (is.null(name) || !nzchar(name)) name <- "haplobreed"
  if (!name %in% names(hbr_palettes))
    stop("Unknown palette '", name, "'. Available: ",
         paste(hbr_palette_names(), collapse = ", "), call. = FALSE)

  base <- hbr_palettes[[name]]
  if (is.null(n)) return(base)
  n <- as.integer(n)
  if (n < 1L) return(character(0))
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}


# ============================================================================
# 2. Theme registry -- 10 ggplot2 themes
# ============================================================================

#' Named HaploBreedR ggplot2 themes
#'
#' A named list of ten theme-generating functions. Every theme accepts
#' `base_size`, `base_family`, `axis_text_angle`, and `legend_pos`.
#'
#' @format A named list of functions.
#'
#' @noRd
hbr_themes <- list(

  haplobreed = function(base_size = 12, base_family = "sans",
                        axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
        plot.background    = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
        panel.grid.major   = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.3),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.line          = ggplot2::element_line(colour = "#222222", linewidth = 0.5),
        axis.ticks         = ggplot2::element_line(colour = "#222222", linewidth = 0.4),
        axis.text          = ggplot2::element_text(colour = "#222222", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(colour = "#222222", face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30,
                                                   hjust = 0.5, colour = "#222222"),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85,
                                                   hjust = 0.5, colour = "#666666"),
        legend.position    = legend_pos,
        legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.95),
        legend.text        = ggplot2::element_text(size = base_size * 0.85),
        legend.key.size    = ggplot2::unit(0.45, "cm"),
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_rect(fill = "#F5F5F5", colour = NA),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  minimal = function(base_size = 12, base_family = "sans",
                     axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.grid.major   = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.line          = ggplot2::element_blank(),
        axis.ticks         = ggplot2::element_blank(),
        axis.text          = ggplot2::element_text(colour = "#333333", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_blank(),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  classic = function(base_size = 12, base_family = "sans",
                     axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        axis.text          = ggplot2::element_text(colour = "#222222", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_rect(fill = "#F0F0F0", colour = "#333333"),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  bw = function(base_size = 12, base_family = "sans",
                axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.border       = ggplot2::element_rect(colour = "#000000", fill = NA, linewidth = 0.6),
        panel.grid.major   = ggplot2::element_line(colour = "#DDDDDD", linewidth = 0.3),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text          = ggplot2::element_text(colour = "#000000", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#333333"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_rect(fill = "#F0F0F0", colour = "#000000"),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  light = function(base_size = 12, base_family = "sans",
                   axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_light(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
        panel.grid.major   = ggplot2::element_line(colour = "#EAEAEA", linewidth = 0.35),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text          = ggplot2::element_text(colour = "#444444", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_rect(fill = "#F0F0F0", colour = NA),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  dark = function(base_size = 12, base_family = "sans",
                  axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_dark(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#222222", colour = NA),
        plot.background    = ggplot2::element_rect(fill = "#1A1A1A", colour = NA),
        panel.grid.major   = ggplot2::element_line(colour = "#3A3A3A", linewidth = 0.35),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text          = ggplot2::element_text(colour = "#EAEAEA", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(colour = "#EAEAEA", face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30,
                                                   hjust = 0, colour = "#FFFFFF"),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#BBBBBB"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        legend.text        = ggplot2::element_text(colour = "#EAEAEA", size = base_size * 0.85),
        legend.title       = ggplot2::element_text(colour = "#FFFFFF", face = "bold", size = base_size * 0.95),
        strip.background   = ggplot2::element_rect(fill = "#2E2E2E", colour = NA),
        strip.text         = ggplot2::element_text(colour = "#FFFFFF", face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  void = function(base_size = 12, base_family = "sans",
                  axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_void(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0.5),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, hjust = 0.5, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.95),
        legend.text        = ggplot2::element_text(size = base_size * 0.85),
        plot.margin        = ggplot2::margin(4, 6, 4, 6))
  },

  linedraw = function(base_size = 12, base_family = "sans",
                      axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_linedraw(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
        panel.border       = ggplot2::element_blank(),
        panel.grid.major   = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.line          = ggplot2::element_line(colour = "#222222", linewidth = 0.6),
        axis.ticks         = ggplot2::element_line(colour = "#222222", linewidth = 0.4),
        axis.text          = ggplot2::element_text(colour = "#222222", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.30, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        legend.key         = ggplot2::element_blank(),
        strip.background   = ggplot2::element_blank(),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(8, 10, 8, 10))
  },

  economist = function(base_size = 12, base_family = "sans",
                       axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
        plot.background    = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
        panel.grid.major.y = ggplot2::element_line(colour = "#D9D9D9", linewidth = 0.5),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.line          = ggplot2::element_line(colour = "#222222", linewidth = 0.7),
        axis.ticks         = ggplot2::element_blank(),
        axis.text          = ggplot2::element_text(colour = "#222222", size = base_size * 0.9),
        axis.title         = ggplot2::element_text(face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.35, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.85, colour = "#666666"),
        legend.position    = legend_pos,
        legend.background  = ggplot2::element_blank(),
        legend.key         = ggplot2::element_blank(),
        strip.background   = ggplot2::element_blank(),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(10, 12, 10, 12))
  },

  scientific = function(base_size = 12, base_family = "sans",
                        axis_text_angle = 0, legend_pos = "right") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
        plot.background    = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
        panel.grid.major   = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.line          = ggplot2::element_line(colour = "#222222", linewidth = 0.5),
        axis.ticks         = ggplot2::element_line(colour = "#222222", linewidth = 0.4),
        axis.text          = ggplot2::element_text(colour = "#222222", size = base_size * 0.85),
        axis.title         = ggplot2::element_text(colour = "#222222", face = "bold", size = base_size * 1.05),
        axis.text.x        = ggplot2::element_text(angle = axis_text_angle,
                                                   hjust = if (axis_text_angle != 0) 1 else 0.5),
        plot.title         = ggplot2::element_text(face = "bold", size = base_size * 1.20, hjust = 0),
        plot.subtitle      = ggplot2::element_text(size = base_size * 0.80, colour = "#666666"),
        legend.position    = legend_pos,
        legend.title       = ggplot2::element_blank(),
        legend.text        = ggplot2::element_text(size = base_size * 0.85),
        legend.key.size    = ggplot2::unit(0.35, "cm"),
        legend.background  = ggplot2::element_blank(),
        strip.background   = ggplot2::element_blank(),
        strip.text         = ggplot2::element_text(face = "bold", size = base_size * 0.90),
        plot.margin        = ggplot2::margin(6, 8, 6, 8))
  }
)

#' Names of the built-in ggplot2 themes
#'
#' @return A character vector of theme names accepted by [hbr_theme()].
#'
#' @examples
#' hbr_theme_names()
#'
#' @export
hbr_theme_names <- function() names(hbr_themes)

#' Return a named HaploBreedR ggplot2 theme
#'
#' @param name Theme name. See [hbr_theme_names()].
#' @param base_size Base font size in pt.
#' @param base_family Base font family.
#' @param axis_text_angle Rotation for x-axis labels (degrees).
#' @param legend_pos Legend position.
#'
#' @return A `ggplot2` theme object.
#'
#' @examples
#' hbr_theme("minimal")
#' hbr_theme("haplobreed", base_size = 10)
#'
#' @export
hbr_theme <- function(name = "haplobreed",
                      base_size = 12,
                      base_family = "sans",
                      axis_text_angle = 0,
                      legend_pos = "right") {
  if (is.null(name) || !nzchar(name)) name <- "haplobreed"
  if (!name %in% names(hbr_themes))
    stop("Unknown theme '", name, "'. Available: ",
         paste(hbr_theme_names(), collapse = ", "), call. = FALSE)
  hbr_themes[[name]](base_size = base_size,
                     base_family = base_family,
                     axis_text_angle = axis_text_angle,
                     legend_pos = legend_pos)
}

#' Default HaploBreedR theme (backward-compatible alias)
#'
#' @inheritParams hbr_theme
#'
#' @return A `ggplot2` theme object equivalent to
#'   `hbr_theme("haplobreed", ...)`.
#'
#' @noRd
hbr_theme_default <- function(base_size = 12, base_family = "sans",
                              axis_text_angle = 0, legend_pos = "right") {
  hbr_theme("haplobreed",
            base_size = base_size, base_family = base_family,
            axis_text_angle = axis_text_angle, legend_pos = legend_pos)
}


# ============================================================================
# 3. Semantic layer -- palette + theme -> role-based colours
# ============================================================================

#' Test whether a theme is dark
#'
#' @param theme Theme name.
#'
#' @return `TRUE` only for the `"dark"` theme.
#'
#' @noRd
.hbr_theme_is_dark <- function(theme) isTRUE(theme == "dark")

#' Palette anchors for a named palette
#'
#' Returns the primary, secondary, and darkened variants for a named
#' palette: `superior`, `other`, `superior_dark`, `other_dark`, `full`.
#'
#' @param palette Palette name.
#'
#' @return A named list of hex colours.
#'
#' @noRd
.hbr_palette_anchors <- function(palette = "haplobreed") {
  pal <- hbr_palette(palette)
  sup <- pal[1]
  oth <- pal[2]
  darken <- function(hex, f = 0.45) {
    rgb <- grDevices::col2rgb(hex) / 255
    grDevices::rgb(rgb[1] * f, rgb[2] * f, rgb[3] * f)
  }
  list(superior      = sup,
       other         = oth,
       superior_dark = darken(sup, 0.45),
       other_dark    = darken(oth, 0.45),
       full          = pal)
}

#' IUPAC allele cell colours for the haplotype table
#'
#' Returns four hex codes for A / C / G / T, derived from the current
#' palette and theme.
#'
#' @param theme One of [hbr_theme_names()].
#' @param palette One of [hbr_palette_names()].
#'
#' @return Named character vector `c(A=, C=, G=, T=)`.
#'
#' @noRd
hbr_cell_colours <- function(theme = "haplobreed", palette = "haplobreed") {

  if (identical(theme, "bw")) {
    return(c(A = "#F0F0F0", C = "#BFBFBF",
             G = "#7F7F7F", T = "#404040"))
  }

  anchor_h <- c(A = 130, C = 210, G = 35, T = 0)
  anchor_s <- c(A = 0.35, C = 0.40, G = 0.55, T = 0.40)

  pal      <- hbr_palette(palette)
  pal_hsv  <- grDevices::rgb2hsv(grDevices::col2rgb(pal[1]))
  family_h <- pal_hsv["h", 1] * 360

  blend_h <- function(h) {
    d <- family_h - h
    if (d >  180) d <- d - 360
    if (d < -180) d <- d + 360
    w <- 0.30 * exp(-(abs(d) / 120)^2)
    (h + w * d) %% 360
  }

  is_dark <- .hbr_theme_is_dark(theme)
  cfg <- if (is_dark) {
    list(v = 0.58, s_mul = 1.35)
  } else if (identical(theme, "economist")) {
    list(v = 0.90, s_mul = 0.80)
  } else if (identical(theme, "light")) {
    list(v = 0.95, s_mul = 0.85)
  } else if (identical(theme, "classic")) {
    list(v = 0.93, s_mul = 0.95)
  } else {
    list(v = 0.92, s_mul = 1.00)
  }

  vapply(names(anchor_h), function(base) {
    grDevices::hsv(h = blend_h(anchor_h[base]) / 360,
                   s = min(1, anchor_s[base] * cfg$s_mul),
                   v = cfg$v)
  }, character(1), USE.NAMES = TRUE)
}

#' Theme surface colour list
#'
#' Returns a list of surface colours (panel / plot background, text,
#' rules, table colours, grid colours) for the given theme, matching the
#' light or dark visual style.
#'
#' @param theme Theme name.
#'
#' @return A named list of hex colours and a logical `is_dark` flag.
#'
#' @noRd
.hbr_theme_surface <- function(theme = "haplobreed") {
  is_dark <- .hbr_theme_is_dark(theme)
  list(
    is_dark            = is_dark,
    panel_bg           = if (is_dark) "#222222" else "#FFFFFF",
    plot_bg            = if (is_dark) "#1A1A1A" else "#FFFFFF",
    text_primary       = if (is_dark) "#EAEAEA" else "#222222",
    text_secondary     = if (is_dark) "#BBBBBB" else "#555555",
    text_muted         = if (is_dark) "#888888" else "#999999",
    grid_major         = if (is_dark) "#3A3A3A" else "#EEEEEE",
    grid_minor         = if (is_dark) "#2E2E2E" else "#F8F8F8",
    table_bg           = if (is_dark) "#2E2E2E" else "#FAFAFA",
    table_header       = if (is_dark) "#3A3A3A" else "#F0F0F0",
    table_superior_row = if (is_dark) "#4A2E10" else "#FFF3E0",
    rule_strong        = if (is_dark) "#EAEAEA" else "#222222",
    rule_light         = if (is_dark) "#666666" else "#CCCCCC")
}

#' Resolve a semantic role to a hex colour
#'
#' Central lookup used by all plotting helpers.
#'
#' @param role Colour role name.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param n Number of colours for `role = "palette_seq"`.
#'
#' @return Hex colour(s), or `NULL` if the role is unknown.
#'
#' @noRd
hbr_semantic <- function(role,
                         palette = "haplobreed",
                         theme   = "haplobreed",
                         n       = NULL) {
  pal   <- .hbr_palette_anchors(palette)
  srf   <- .hbr_theme_surface(theme)
  cells <- hbr_cell_colours(theme, palette)
  dark  <- srf$is_dark

  switch(role,
         superior          = pal$superior,
         superior_dark     = pal$superior_dark,
         other             = pal$other,
         other_dark        = pal$other_dark,
         accent            = pal$superior,
         palette_seq       = hbr_palette(palette, n = n),

         gwas_lead         = if (dark) "#EF5350" else "#B2182B",
         gwas_window       = if (dark) "#3A2A2A" else "#FCE4D6",
         ld_gradient       = if (dark)
           c("#1A1A1A", "#7A2F00", "#CC5500", "#EF7B30", "#FFC080")
         else
           c("#F2F2F2", "#FEE08B", "#FDAE61", "#F46D43", "#A50026"),
         qq_line           = if (dark) "#64B5F6" else "#2166AC",
         qq_reference      = if (dark) "#777777" else "#666666",
         chrom_alt         = if (dark) c("#3B6EA5", "#A57B2F")
         else c("#2C7FB8", "#F0A830"),
         gene_backbone     = if (dark) "#CCCCCC" else "#3A3A3A",

         text_primary      = srf$text_primary,
         text_secondary    = srf$text_secondary,
         text_muted        = srf$text_muted,
         panel_bg          = srf$panel_bg,
         plot_bg           = srf$plot_bg,
         table_bg          = srf$table_bg,
         table_header      = srf$table_header,
         table_superior    = srf$table_superior_row,
         rule_strong       = srf$rule_strong,
         rule_light        = srf$rule_light,
         grid_major        = srf$grid_major,
         grid_minor        = srf$grid_minor,

         cell_A            = cells[["A"]],
         cell_C            = cells[["C"]],
         cell_G            = cells[["G"]],
         cell_T            = cells[["T"]],
         cell_colours      = cells,

         NULL)
}

#' Consistent p-value formatter for figure annotations
#'
#' Formats a p-value as scientific notation below 1e-3, three
#' significant digits below 0.05, and two digits otherwise.
#'
#' @param p Numeric vector of p-values.
#'
#' @return Character vector, same length as `p`.
#'
#' @noRd
hbr_format_p <- function(p) {
  fmt_one <- function(x) {
    if (!is.finite(x))  return(NA_character_)
    if (x < 1e-3)       return(sprintf("%.1e", x))
    if (x < 0.05)       return(sprintf("%.3f", x))
    sprintf("%.2f", x)
  }
  vapply(p, fmt_one, character(1), USE.NAMES = FALSE)
}

#' Validate palette and theme names
#'
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return `TRUE`, invisibly.
#'
#' @noRd
.hbr_check_palette_theme <- function(palette, theme) {
  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'. Available: ",
         paste(hbr_palette_names(), collapse = ", "), call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'. Available: ",
         paste(hbr_theme_names(), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}


# ============================================================================
# 4. Cached defaults
# ============================================================================

#' Default palette anchors (cached)
#' @noRd
.hbr_pal_default  <- .hbr_palette_anchors("haplobreed")

#' Default theme surface (cached)
#' @noRd
.hbr_srf_default  <- .hbr_theme_surface("haplobreed")

#' Default cell colours (cached)
#' @noRd
.hbr_cell_default <- hbr_cell_colours("haplobreed")

#' Default semantic colour list
#'
#' Named list of role-based colours used across plot functions when no
#' palette or theme argument is supplied.
#'
#' @noRd
hbr_palette_default <- list(
  superior        = .hbr_pal_default$superior,
  superior_light  = .hbr_srf_default$table_superior_row,
  other           = .hbr_pal_default$other,
  other_light     = "#E3F2FD",

  exon            = "#762A83",
  exon_light      = "#F3E5F5",
  utr             = "#999999",
  utr_light       = "#F5F5F5",
  intron          = "#666666",
  promoter        = "#FFFFFF",

  aa_superior     = .hbr_pal_default$superior,
  aa_other        = "#56B4E9",
  aa_synonymous   = "#999999",

  gwas_lead       = "#B2182B",
  gwas_window     = "#F5F5F5",
  gwas_dash       = "#B2182B",

  table_header    = .hbr_srf_default$table_header,
  table_superior  = .hbr_srf_default$table_superior_row,
  table_alternate = .hbr_srf_default$table_bg,
  table_star      = .hbr_pal_default$superior,

  text_primary    = .hbr_srf_default$text_primary,
  text_secondary  = .hbr_srf_default$text_secondary,

  panel_bg        = .hbr_srf_default$panel_bg,
  plot_bg         = .hbr_srf_default$plot_bg,
  grid_major      = .hbr_srf_default$grid_major,
  grid_minor      = .hbr_srf_default$grid_minor,

  chr_odd         = "#2C7FB8",
  chr_even        = "#F0A830",
  chr_odd_light   = "#7FB8D8",
  chr_even_light  = "#F5C670",

  zoom_box        = "#B2182B",
  zoom_shade      = "#FCE4D6",
  zoom_gene       = .hbr_pal_default$superior,
  zoom_snp        = .hbr_pal_default$superior_dark,

  qq_line         = "#2166AC",
  qq_reference    = "#666666",
  qq_ci           = "#E3F2FD",

  cell_A          = .hbr_cell_default[["A"]],
  cell_C          = .hbr_cell_default[["C"]],
  cell_G          = .hbr_cell_default[["G"]],
  cell_T          = .hbr_cell_default[["T"]])
