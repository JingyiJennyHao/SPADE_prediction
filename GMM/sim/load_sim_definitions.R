# Load definitions only; do not execute example simulations.
lines <- readLines("gmm_sim.Rmd", warn = FALSE)
inside <- FALSE
code <- character()
for (line in lines) {
  if (grepl("^```\\{r", line)) {
    inside <- !grepl("eval\\s*=\\s*FALSE", line)
    next
  }
  if (grepl("^```", line)) { inside <- FALSE; next }
  if (inside) code <- c(code, line)
}
eval(parse(text = code), envir = globalenv())
