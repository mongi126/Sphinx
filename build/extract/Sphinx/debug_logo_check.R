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

write_log <- function(hypothesis_id, location, message, data = list()) {
  line <- sprintf(
    '{"sessionId":"%s","runId":"logo-check","hypothesisId":"%s","location":"%s","message":"%s","data":%s,"timestamp":%s}',
    session_id, hypothesis_id, gsub('"', '\\\\"', location),
    gsub('"', '\\\\"', message), to_json(data), as.integer(Sys.time()) * 1000
  )
  cat(line, "\n", file = log_path, append = TRUE)
}

pkg_root <- "d:/Sphinx"
candidates <- c(
  file.path(pkg_root, "logo.svg"),
  file.path(pkg_root, "man", "figures", "logo.svg"),
  file.path(pkg_root, "logo.png"),
  file.path(pkg_root, "man", "figures", "logo.png")
)

# #region agent log
write_log("A", "debug_logo_check.R:38", "pkgdown find_logo candidate paths", list(
  candidates = candidates,
  exists = file.exists(candidates)
))
# #endregion

# #region agent log
write_log("B", "debug_logo_check.R:45", "User logo location docs/logo.png", list(
  docs_logo_exists = file.exists(file.path(pkg_root, "docs", "logo.png")),
  docs_logo_size = if (file.exists(file.path(pkg_root, "docs", "logo.png"))) {
    file.info(file.path(pkg_root, "docs", "logo.png"))$size
  } else {
    NA
  }
))
# #endregion

# #region agent log
logo_found <- tryCatch(
  pkgdown:::find_logo(pkg_root),
  error = function(e) structure(NA_character_, error = conditionMessage(e))
)
write_log("C", "debug_logo_check.R:59", "pkgdown find_logo result", list(
  found = !is.na(logo_found) && nzchar(logo_found),
  path = if (is.na(logo_found)) NA_character_ else as.character(logo_found),
  error = attr(logo_found, "error")
))
# #endregion

# #region agent log
favicons_err <- tryCatch(
  { pkgdown::build_favicons(pkg_root, overwrite = TRUE); NULL },
  error = function(e) conditionMessage(e)
)
write_log("D", "debug_logo_check.R:70", "build_favicons outcome", list(
  success = is.null(favicons_err),
  error = favicons_err
))
# #endregion

cat("Logo diagnostic complete.\n")
