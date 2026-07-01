# Debug: rcmdcheck - session 48a66c
pkg_root <- "d:/Sphinx"
log_path <- "d:/Sphinx/debug-48a66c.log"
session_id <- "48a66c"
to_json <- function(x) {
  if (is.logical(x)) return(ifelse(x, "true", "false"))
  if (is.null(x)) return("null")
  if (is.numeric(x) && length(x) == 1) return(as.character(x))
  if (is.character(x) && length(x) == 1) return(sprintf('"%s"', gsub('"', '\\\\"', x)))
  if (is.character(x)) {
    items <- vapply(x, function(v) sprintf('"%s"', gsub('"', '\\\\"', v)), character(1))
    return(sprintf("[%s]", paste(items, collapse = ",")))
  }
  if (is.list(x)) {
    keys <- names(x)
    if (is.null(keys)) keys <- as.character(seq_along(x))
    pairs <- mapply(function(k, v) sprintf('"%s":%s', k, to_json(v)), keys, x, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    return(sprintf("{%s}", paste(pairs, collapse = ",")))
  }
  sprintf('"%s"', as.character(x))
}
write_log <- function(hypothesis_id, message, data = list()) {
  line <- sprintf(
    '{"sessionId":"%s","runId":"cran-prep","hypothesisId":"%s","location":"debug_rcmdcheck.R","message":"%s","data":%s,"timestamp":%s}',
    session_id, hypothesis_id, gsub('"', '\\\\"', message), to_json(data), as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}
chk_path <- pkg_root
built <- tryCatch(
  devtools::build(pkg_root, path = file.path(pkg_root, "build"), quiet = TRUE),
  error = function(e) NULL
)
if (!is.null(built)) chk_path <- built
chk <- rcmdcheck::rcmdcheck(chk_path, args = c("--as-cran", "--no-manual"), quiet = TRUE)
write_log("A", "rcmdcheck summary", list(
  errors = length(chk$errors), warnings = length(chk$warnings), notes = length(chk$notes),
  error_head = if (length(chk$errors)) chk$errors[[1]] else NA,
  warning_head = if (length(chk$warnings)) chk$warnings[[1]] else NA,
  note_head = if (length(chk$notes)) chk$notes[[1]] else NA
))
cat("Errors:", length(chk$errors), "\nWarnings:", length(chk$warnings), "\nNotes:", length(chk$notes), "\n")
for (e in chk$errors) cat("\n--- ERROR ---\n", e, "\n")
for (w in chk$warnings) cat("\n--- WARNING ---\n", w, "\n")
for (n in chk$notes) cat("\n--- NOTE ---\n", n, "\n")
