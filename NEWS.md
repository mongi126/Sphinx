# Sphinx 1.1.1

## Documentation

- Example FOV switched to **reg055_A** from Schurch *et al.* CODEX CRC
  (2020): largest sample (~3.9k cells) with `ClusterName` annotations

# Sphinx 1.1.0

## Visualization

- Soft candy qualitative palette (one hue family per color) with stable
  `assign_celltype_colors()` mapping across plots
- Distance / interaction heatmaps use a soft pink sequential ramp; composition
  heatmaps use candy blue–pink diverging colors
- Larger publication fonts (minimum 8 pt); titles not forced bold
- `visualize_spatial_network()`: opaque points with white halo, optional local
  zoom (`zoom_center` / `zoom_radius` or `xlim`/`ylim`) showing real cell–cell edges
- `visualize_interaction_network()`: legend keeps cell-type dots only
- `visualize_interaction_heatmap()`: smaller cells, numbers off by default
- Volcano plots truncate `-log10(adj.P)` (default 50) with triangle markers

## Documentation

- Tutorials rewritten with annotated public spatial proteomics FOVs
- pkgdown site content refreshed to match the updated API

## Maintenance

- `LazyData` disabled (no packaged datasets)
- Expanded `globalVariables` for data.table NSE helpers
- Package cleanup of build/tooling clutter

# Sphinx 1.0.2

## Changes

- Package maintenance release and rebuild

# Sphinx 1.0.0

## Major changes

- Initial CRAN release
- Four core modules: preprocessing, annotation, spatial network, functional analysis
