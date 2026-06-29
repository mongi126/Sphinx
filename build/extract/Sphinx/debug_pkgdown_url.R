# Debug diagnostic (NOT part of package) - session 48a66c
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

write_log <- function(hypothesis_id, location, message, data = list(), run_id = "pkgdown-url") {
  line <- sprintf(
    '{"sessionId":"%s","runId":"%s","hypothesisId":"%s","location":"%s","message":"%s","data":%s,"timestamp":%s}',
    session_id, run_id, hypothesis_id, gsub('"', '\\\\"', location),
    gsub('"', '\\\\"', message), to_json(data), as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

pkg_root <- "d:/Sphinx"
pkg <- pkgdown::as_pkgdown(pkg_root)
desc <- read.dcf(file.path(pkg_root, "DESCRIPTION"))
desc_urls <- c(desc[1, "URL"], desc[1, "BugReports"])

# #region agent log
write_log("A", "debug_pkgdown_url.R:36", "pkg meta url vector", list(
  meta_url = unname(pkg$meta$url),
  meta_url_length = length(pkg$meta$url)
))
# #endregion

# #region agent log
write_log("B", "debug_pkgdown_url.R:43", "DESCRIPTION urls", list(
  desc_url = unname(desc[1, "URL"]),
  bugreports = unname(desc[1, "BugReports"]),
  desc_urls_length = length(desc_urls)
))
# #endregion

# #region agent log
cond <- tryCatch(
  !pkg$meta$url %in% desc_urls,
  error = function(e) structure(NA, error = conditionMessage(e))
)
write_log("C", "debug_pkgdown_url.R:54", "Failing if-condition evaluation", list(
  condition_length = if (is.logical(cond)) length(cond) else NA,
  condition_values = if (is.logical(cond)) cond else NA,
  error = attr(cond, "error")
))
# #endregion

# #region agent log
write_log("D", "debug_pkgdown_url.R:64", "Post-fix meta url", list(
  meta_url_length = length(pkg$meta$url),
  meta_url = unname(pkg$meta$url)
), run_id = "post-fix")
# #endregion

cat("pkgdown url diagnostic complete.\n")
