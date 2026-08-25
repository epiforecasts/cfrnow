## Generate the cfrnow hex logo.
##
## Run from the package root:
##   Rscript data-raw/logo.R
##   inkscape man/figures/logo.svg -o man/figures/logo.png -w 240
##
## The graphic is the line list the model reads: one row per case, running
## from onset to outcome -- died, recovered, or still open. The open rows are
## the whole problem, since a case fatality ratio computed now has to account
## for them.

## --- canvas ----------------------------------------------------------------

width <- 1732L
height <- 2000L
hex <- sprintf(
  "%d,0 %d,500 %d,1500 %d,%d 0,1500 0,500",
  width %/% 2L, width, width, width %/% 2L, height
)

## --- palette ---------------------------------------------------------------

bg <- "#FFFFFF"
bg_foot <- "#FDF2F5"
border <- "#A82A44"
ink <- "#16212A"
died_col <- "#A82A44"
recovered_col <- "#2E7DB8"
open_col <- "#8A98A8"
bar_col <- "#D3DDE7"
rule_col <- "#E2E8EF"

## --- geometry --------------------------------------------------------------

x_start <- 250
x_end <- 1430
ybase <- 1190

flip <- function(y) {
  ybase - y
}

## one row per case: how far it has run, and how it ended
run <- c(0.80, 0.38, 0.62, 0.26, 0.92, 0.50)
fate <- c("died", "recovered", "recovered", "died", "open", "open")
ys <- seq(740, 150, length.out = length(run))

## --- row parts -------------------------------------------------------------

case_bar <- function(x1, x2, y) {
  sprintf(
    paste0(
      '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" ',
      'stroke="%s" stroke-width="34" stroke-linecap="round"/>'
    ),
    x1, flip(y), x2, flip(y), bar_col
  )
}

## died: a filled marker, the only thing that counts in the numerator
died_mark <- function(x, y) {
  sprintf(
    '<circle cx="%.0f" cy="%.0f" r="42" fill="%s"/>',
    x, flip(y), died_col
  )
}

## recovered: an open marker, heavy enough to hold at small sizes
recovered_mark <- function(x, y) {
  sprintf(
    paste0(
      '<circle cx="%.0f" cy="%.0f" r="34" fill="%s" ',
      'stroke="%s" stroke-width="24"/>'
    ),
    x, flip(y), bg, recovered_col
  )
}

## still open: the row simply keeps going
open_mark <- function(x, y) {
  sprintf(
    paste0(
      '<path d="M %.0f %.0f l 86 0 m -34 -34 l 34 34 l -34 34" ',
      'fill="none" stroke="%s" stroke-width="22" ',
      'stroke-linecap="round" stroke-linejoin="round"/>'
    ),
    x, flip(y), open_col
  )
}

rows <- character()
for (i in seq_along(run)) {
  x_stop <- x_start + run[i] * (x_end - x_start)
  mark <- switch(
    fate[i],
    died = died_mark(x_stop, ys[i]),
    recovered = recovered_mark(x_stop, ys[i]),
    open = open_mark(x_stop, ys[i])
  )
  rows <- c(rows, case_bar(x_start, x_stop, ys[i]), mark)
}

## --- assemble --------------------------------------------------------------

svg <- c(
  sprintf(
    paste0(
      '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" ',
      'viewBox="0 0 %d %d">'
    ),
    width, height, width, height
  ),
  "<defs>",
  '<linearGradient id="g" x1="0" y1="0" x2="0" y2="1">',
  sprintf('<stop offset="0" stop-color="%s"/>', bg),
  sprintf('<stop offset="1" stop-color="%s"/>', bg_foot),
  "</linearGradient>",
  sprintf('<clipPath id="clip"><polygon points="%s"/></clipPath>', hex),
  "</defs>",
  sprintf('<polygon points="%s" fill="url(#g)"/>', hex),
  ## clipped so the chart can never bleed past the border
  '<g clip-path="url(#clip)">',
  ## onset, where every row starts
  sprintf(
    paste0(
      '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" ',
      'stroke="%s" stroke-width="12" stroke-linecap="round"/>'
    ),
    x_start - 40, flip(90), x_start - 40, flip(800), rule_col
  ),
  rows,
  "</g>",
  sprintf(
    paste0(
      '<text x="%d" y="1610" text-anchor="middle" ',
      'font-family="Fira Sans, DejaVu Sans, Helvetica, Arial, sans-serif" ',
      'font-size="252" font-weight="600" letter-spacing="2" fill="%s">',
      'cfr<tspan fill="%s">now</tspan></text>'
    ),
    width %/% 2L, ink, died_col
  ),
  sprintf(
    paste0(
      '<polygon points="%s" fill="none" stroke="%s" ',
      'stroke-width="44" stroke-linejoin="round"/>'
    ),
    hex, border
  ),
  "</svg>"
)

writeLines(svg, file.path("man", "figures", "logo.svg"))
