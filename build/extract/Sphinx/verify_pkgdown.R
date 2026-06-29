# Post-fix verification - session 48a66c
log_path <- "d:/Sphinx/debug-48a66c.log"
session_id <- "48a66c"

to_json <- function(x) {
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
    pairs <- mapply(function(k, v) sprintf('"%s":%s', k, to_json(v)), keys, x, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    return(sprintf("{%s}", paste(pairs, collapse = ",")))
  }
  sprintf('"%s"', as.character(x))
}

write_log <- function(hypothesis_id, location, message, data = list()) {
  line <- sprintf(
    '{"sessionId":"%s","runId":"post-fix","hypothesisId":"%s","location":"%s","message":"%s","data":%s,"timestamp":%s}',
    session_id, hypothesis_id, gsub('"', '\\\\"', location),
    gsub('"', '\\\\"', message), to_json(data), as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

pkg <- pkgdown::as_pkgdown("d:/Sphinx")
desc <- read.dcf("d:/Sphinx/DESCRIPTION")
desc_urls <- c(desc[1, "URL"], desc[1, "BugReports"])
fp <- if (requireNamespace("rmarkdown", quietly = TRUE)) rmarkdown::find_pandoc() else list(version = "0", dir = NULL)

# #region agent log
write_log("A", "verify_pkgdown.R:38", "URL condition now scalar", list(
  meta_url_length = length(pkg$meta$url),
  meta_url = unname(pkg$meta$url),
  condition_ok = length(!pkg$meta$url %in% desc_urls) == 1
))
# #endregion

# #region agent log
write_log("B", "verify_pkgdown.R:46", "Pandoc availability", list(
  sys_which = as.character(Sys.which("pandoc")),
  find_pandoc_dir = as.character(fp$dir),
  find_pandoc_version = as.character(fp$version)
))
# #endregion

# #region agent log
url_check_err <- tryCatch({
  pkgdown:::check_urls(pkg)
  NULL
}, error = function(e) conditionMessage(e))
write_log("C", "verify_pkgdown.R:57", "check_urls result", list(
  success = is.null(url_check_err),
  error = url_check_err
))
# #endregion

cat("Verification complete.\n")
