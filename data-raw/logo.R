## Generate the cfrnow hex logo.
##
## Run from the package root:
##   Rscript data-raw/logo.R
##   inkscape man/figures/logo.svg -o man/figures/logo.png -w 240
## Concept: the line list the model reads. One row per case, running from onset
## to outcome -- died, recovered, or still open. The open rows are the whole
## problem: a case fatality ratio computed now has to account for them.

W <- 1732L; H <- 2000L                   # pointy-top hex canvas

## --- palette ---------------------------------------------------------------
bg       <- "#FFFFFF"
bg_foot  <- "#FDF2F5"
border   <- "#A82A44"
ink      <- "#16212A"
died_col <- "#A82A44"
recov_col<- "#2E7DB8"
open_col <- "#8A98A8"
bar_col  <- "#D3DDE7"
rule_col <- "#E2E8EF"

## --- geometry ---------------------------------------------------------------
x0 <- 250; x_end <- 1430
ybase <- 1190
flip <- function(y) ybase - y

## one row per case: how far it has run, and how it ended
len  <- c(0.80, 0.38, 0.62, 0.26, 0.92, 0.50)
fate <- c("died", "recovered", "recovered", "died", "open", "open")
ys   <- seq(740, 150, length.out = length(len))

bar <- function(x1, x2, y) sprintf(
  '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="%s" stroke-width="34" stroke-linecap="round"/>',
  x1, flip(y), x2, flip(y), bar_col)

rows <- character()
for (i in seq_along(len)) {
  xe <- x0 + len[i] * (x_end - x0)
  rows <- c(rows, bar(x0, xe, ys[i]))
  rows <- c(rows, switch(fate[i],
    ## died: a filled marker, the only thing that counts in the numerator
    died = sprintf('<circle cx="%.0f" cy="%.0f" r="42" fill="%s"/>', xe, flip(ys[i]), died_col),
    ## recovered: an open marker, heavy enough to hold at small sizes
    recovered = sprintf('<circle cx="%.0f" cy="%.0f" r="34" fill="%s" stroke="%s" stroke-width="24"/>',
                        xe, flip(ys[i]), bg, recov_col),
    ## still open: the row simply keeps going
    open = sprintf('<path d="M %.0f %.0f l 86 0 m -34 -34 l 34 34 l -34 34" fill="none" stroke="%s" stroke-width="22" stroke-linecap="round" stroke-linejoin="round"/>',
                   xe, flip(ys[i]), open_col)))
}

hex <- sprintf("%d,0 %d,500 %d,1500 %d,%d 0,1500 0,500", W %/% 2L, W, W, W %/% 2L, H)

svg <- c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', W, H, W, H),
  '<defs><linearGradient id="g" x1="0" y1="0" x2="0" y2="1">',
  sprintf('<stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/>', bg, bg_foot),
  '</linearGradient>',
  sprintf('<clipPath id="clip"><polygon points="%s"/></clipPath>', hex),
  '</defs>',
  sprintf('<polygon points="%s" fill="url(#g)"/>', hex),
  '<g clip-path="url(#clip)">',
  ## onset, where every row starts
  sprintf('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="%s" stroke-width="12" stroke-linecap="round"/>',
          x0 - 40, flip(90), x0 - 40, flip(800), rule_col),
  rows,
  sprintf('<text x="%d" y="1610" text-anchor="middle" font-family="Fira Sans, DejaVu Sans, Helvetica, Arial, sans-serif" font-size="252" font-weight="600" letter-spacing="2" fill="%s">cfr<tspan fill="%s">now</tspan></text>',
          W %/% 2L, ink, died_col),
  '</g>',
  sprintf('<polygon points="%s" fill="none" stroke="%s" stroke-width="44" stroke-linejoin="round"/>', hex, border),
  '</svg>')

writeLines(svg, "man/figures/logo.svg")
