# Debug diagnostic script (NOT part of package) - session 48a66c
log_path <- "d:/Sphinx/debug-48a66c.log"
session_id <- "48a66c"

write_log <- function(hypothesis_id, location, message, data = list()) {
  data_json <- jsonlite::toJSON(data, auto_unbox = TRUE, null = "null")
  line <- sprintf(
    '{"sessionId":"%s","runId":"initial","hypothesisId":"%s","location":"%s","message":"%s","data":%s,"timestamp":%s}',
    session_id,
    hypothesis_id,
    gsub('"', '\\\\"', location),
    gsub('"', '\\\\"', message),
    data_json,
    as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

to_json <- function(x) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
  }
  if (is.logical(x)) return(ifelse(x, "true", "false"))
  if (is.null(x)) return("null")
  if (is.character(x) && length(x) == 1) return(sprintf('"%s"', gsub('"', '\\\\"', x)))
  if (is.character(x)) {
    items <- vapply(x, function(v) sprintf('"%s"', gsub('"', '\\\\"', v)), character(1))
    return(sprintf("[%s]", paste(items, collapse = ",")))
  }
  if (is.list(x)) {
    keys <- names(x)
    if (is.null(keys)) keys <- as.character(seq_along(x))
    pairs <- mapply(function(k, v) {
      sprintf('"%s":%s', k, to_json(v))
    }, keys, x, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    return(sprintf("{%s}", paste(pairs, collapse = ",")))
  }
  sprintf('"%s"', as.character(x))
}

write_log <- function(hypothesis_id, location, message, data = list()) {
  line <- sprintf(
    '{"sessionId":"%s","runId":"initial","hypothesisId":"%s","location":"%s","message":"%s","data":%s,"timestamp":%s}',
    session_id,
    hypothesis_id,
    gsub('"', '\\\\"', location),
    gsub('"', '\\\\"', message),
    to_json(data),
    as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

# #region agent log
write_log("A", "debug_pandoc_check.R:55", "System pandoc path", list(
  sys_which_pandoc = as.character(Sys.which("pandoc")),
  path_has_pandoc = grepl("pandoc", Sys.getenv("PATH"), ignore.case = TRUE)
))
# #endregion

# #region agent log
rmarkdown_ok <- requireNamespace("rmarkdown", quietly = TRUE)
pandoc_info <- if (rmarkdown_ok) {
  tryCatch({
    fp <- rmarkdown::find_pandoc()
    list(
      version = as.character(fp$version),
      dir = as.character(fp$dir),
      name = as.character(fp$name),
      type = as.character(fp$type)
    )
  }, error = function(e) list(error = conditionMessage(e)))
} else {
  list(rmarkdown_missing = TRUE)
}
write_log("B", "debug_pandoc_check.R:68", "R rmarkdown pandoc detection", list(
  rmarkdown_installed = rmarkdown_ok,
  find_pandoc = pandoc_info
))
# #endregion

# #region agent log
write_log("C", "debug_pandoc_check.R:75", "Package vignette config", list(
  vignette_builder = unname(read.dcf("DESCRIPTION")[1, "VignetteBuilder"]),
  inst_vignettes_rmd_count = length(list.files("inst/vignettes", pattern = "\\.Rmd$")),
  inst_vignettes_html_count = length(list.files("inst/vignettes", pattern = "\\.html$")),
  rbuildignore_excludes_vignettes = any(grepl("^\\^vignettes\\$", readLines(".Rbuildignore")))
))
# #endregion

# #region agent log
write_log("D", "debug_pandoc_check.R:84", "R bundled pandoc paths", list(
  r_home = R.home(),
  rmarkdown_installed = rmarkdown_ok,
  bundled_pandoc_available = if (rmarkdown_ok) {
    fp <- rmarkdown::find_pandoc()
    !is.null(fp$dir) && nzchar(fp$dir)
  } else {
    FALSE
  }
))
# #endregion

cat("Diagnostic complete. Log written to:", log_path, "\n")
