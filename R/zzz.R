.onAttach <- function(libname, pkgname) {
  version <- as.character(utils::packageVersion("Sphinx"))

  msg <- paste0(
    "+-----------------------------------------------------+\n",
    "| Sphinx v", version, "                                       |\n",
    "| Spatial Proteomics Analysis Toolkit                 |\n",
    "+-----------------------------------------------------+\n\n",
    "Core modules loaded:\n",
    "  - Preprocessing\n",
    "  - Annotation\n",
    "  - Spatial network analysis\n",
    "  - Functional analysis\n",
    "Documentation: help(package = 'Sphinx')\n",
    "Tutorial: vignette('sphinx-intro')\n",
    "Issues: https://github.com/mongi126/Sphinx/issues\n"
  )

  packageStartupMessage(msg)
}

.onUnload <- function(libpath) {
  invisible()
}

utils::globalVariables(c(
  ".", "X", "Y", "Cell_ID", "from", "to", "cluster",
  "Centroid_X", "Centroid_Y", "Protein", "FDR", "Term"
))
