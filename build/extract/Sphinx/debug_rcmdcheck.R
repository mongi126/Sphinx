# Debug: rcmdcheck note counter - session 48a66c
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

write_log <- function(hypothesis_id, message, data = list(), run_id = "rcmdcheck") {
  line <- sprintf(
    '{"sessionId":"%s","runId":"%s","hypothesisId":"%s","location":"debug_rcmdcheck.R","message":"%s","data":%s,"timestamp":%s}',
    session_id, run_id, hypothesis_id, gsub('"', '\\\\"', message), to_json(data), as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

pkg_root <- normalizePath("d:/Sphinx", winslash = "/")

# #region agent log
write_log("A", "Pre-check config", list(
  vignettes_in_buildignore = any(grepl("^\\^vignettes\\$", readLines(file.path(pkg_root, ".Rbuildignore")))),
  inst_vignettes_ignored = any(grepl("inst/vignettes", readLines(file.path(pkg_root, ".Rbuildignore"))))
))
# #endregion

chk <- rcmdcheck::rcmdcheck(pkg_root, quiet = TRUE)

# #region agent log
write_log("B", "rcmdcheck summary", list(
  errors = length(chk$errors),
  warnings = length(chk$warnings),
  notes = length(chk$notes),
  note_texts = unname(chk$notes)
))
# #endregion

cat("Errors:", length(chk$errors), "\n")
cat("Warnings:", length(chk$warnings), "\n")
cat("Notes:", length(chk$notes), "\n")
if (length(chk$notes) > 0) cat(paste(chk$notes, collapse = "\n\n"), "\n")
