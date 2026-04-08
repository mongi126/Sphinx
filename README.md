# Sphinx: Spatial Proteomics Analysis Toolkit <img src="man/figures/logo.png" width="100" align="right">

[![R-CMD-check](https://github.com/mongi126/Sphinx/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mongi126/Sphinx/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.3-blue)](https://www.r-project.org/)

---

## 📖 Overview

<div align="center">
  <img src="man/figures/workflow.png" width="1000">
</div>

**Sphinx** is a comprehensive R package for spatial proteomics data analysis, providing a one-stop solution from quality control to biological mechanism discovery. The toolkit integrates four core modules:

- **🔬 Data Preprocessing** - Standardized preprocessing supporting CODEX, MIBI, CyCIF and other mainstream platforms
- **🧬 Cell Type Annotation** - Interactive manual annotation based on protein expression gradients
- **🌐 Spatial Neighborhood Modeling** - Integration of Delaunay triangulation, kNN, fixed radius, and dynamic window radius strategies
- **📊 Functional Analysis** - Spatial enrichment analysis linking spatial phenotypes to functional pathways (KEGG, GO, CORUM, Reactome)

---

## 🔧 Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install Sphinx
devtools::install_github("mongi126/Sphinx")
```

### From local source
```r
install.packages("Sphinx_1.0.0.tar.gz", repos = NULL, type = "source")
```

Full documentation and vignettes are available at:
👉 https://mongi126.github.io/Sphinx/
