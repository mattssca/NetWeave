# NetWeave

> **N**etwork **E**xpression **T**issue **W**eaving and **A**nnotation **V**isualization **E**ngine

A modular R package for building, annotating, and visualising gene interaction networks. Starting from pathway seed genes, NetWeave expands the network using STRING protein-protein interaction data, integrates subtype-specific gene expression, Human Protein Atlas (HPA) protein evidence, and computes network centrality metrics and subtype-specific hub scores.

---

## Installation

```r
# Install Bioconductor dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("biomaRt", "STRINGdb", "clusterProfiler", 
                       "org.Hs.eg.db", "ComplexHeatmap"))

# Install NetWeave from GitHub
devtools::install_github("yourusername/NetWeave")
```

---

## Pipeline Overview

NetWeave follows a linear, modular pipeline. Each function appends new columns to `node_metrics` and passes it forward — making it easy to stop, inspect, or branch at any step.

```
Pathway genes
     │
     ▼
create_and_expand_network()     → edge_data + node_metrics
     │
     ▼
add_subtype_expression()        → mean_expr_<subtype> columns
     │
     ▼
annotate_low_expressed_genes()  → is_low_expression
     │
     ▼
run_anova_subtype_expression()  → anova_p_value, anova_sig, ...
     │
     ▼
rank_subtype_expression()       → rank_<subtype> columns
     │
     ▼
add_hpa_annotations()           → hpa_normal_*, hpa_cancer_* columns
     │
     ▼
derive_hpa_metrics()            → hpa_expression_pattern, hpa_passes_filter, ...
     │
     ▼
filter_network_nodes()          → filtered node_metrics
     │
     ▼
recalculate_network_metrics()   → updated degree, betweenness, community, ...
     │
     ├──▶ run_go_enrichment()              → GO results per community
     ├──▶ calculate_subtype_hub_scores()   → <subtype>_hub_score columns
     │         └──▶ plot_subtype_hub_scores()
     ├──▶ plot_network()                   → ggraph / visNetwork
     └──▶ create_cytoscape_objects()       → TSV + GraphML export
```

---

## Quick Start

```r
library(NetWeave)

# Load bundled example data
data(sjodahl_2017)
data(sjodahl_2017_subtypes)

# 1. Get pathway seed genes
erbb2_genes <- get_pathway_genes(
  gene = "ERBB2",
  reactome_data = Ensembl2Reactome
)

# 2. Build and expand the STRING network
network <- create_and_expand_network(
  expr_data        = sjodahl_2017,
  pathway_genes    = erbb2_genes,
  seed_gene        = "ERBB2",
  max_added_genes  = 20,
  string_data_dir  = "data/string_db"
)

edge_data    <- network$edge_data
node_metrics <- network$node_metrics

# 3. Add subtype expression
node_metrics <- add_subtype_expression(
  node_metrics       = node_metrics,
  expr_data          = sjodahl_2017,
  this_subtype_vector = sjodahl_2017_subtypes
)

# 4. Flag low-expressed genes
node_metrics <- annotate_low_expressed_genes(
  node_metrics       = node_metrics,
  expr_data          = sjodahl_2017,
  low_expr_threshold = 0.25
)

# 5. ANOVA across subtypes
node_metrics <- run_anova_subtype_expression(
  node_metrics        = node_metrics,
  expr_data           = sjodahl_2017,
  this_subtype_vector = sjodahl_2017_subtypes
)

# 6. Rank expression across subtypes
node_metrics <- rank_subtype_expression(node_metrics)

# 7. Add HPA annotations
node_metrics <- add_hpa_annotations(
  node_metrics       = node_metrics,
  normal_hpa         = hpa_normal_data,
  cancer_hpa         = hpa_cancer_data,
  this_tissue_normal = "urinary bladder",
  this_tissue_cancer = "urothelial cancer"
)

node_metrics <- derive_hpa_metrics(node_metrics)

# 8. Filter nodes
filtered <- filter_network_nodes(
  node_metrics         = node_metrics,
  tissue_filter        = "either",
  remove_low_expressed = TRUE,
  filter_anova_sig     = TRUE
)

# 9. Recalculate metrics on the filtered subnetwork
filtered <- recalculate_network_metrics(
  node_metrics = filtered,
  edge_data    = edge_data
)

# 10. GO enrichment per community
go_results <- run_go_enrichment(filtered, ont = "BP")

# 11. Hub scores + heatmap
filtered <- calculate_subtype_hub_scores(filtered)
plot_subtype_hub_scores(filtered, top_n = 20)

# 12. Network plot
plot_network(
  node_metrics = filtered,
  edge_data    = edge_data,
  color_by     = "community",
  size_by      = "degree",
  label_nodes  = "seed"
)

# 13. Export to Cytoscape
create_cytoscape_objects(
  node_metrics = filtered,
  edge_data    = edge_data,
  out_dir      = "output/cytoscape",
  file_prefix  = "erbb2_network"
)
```

---

## Function Reference

| Function | Description |
|---|---|
| `get_pathway_genes()` | Retrieve genes for a pathway or seed gene from Ensembl/Reactome |
| `create_and_expand_network()` | Build STRING PPI network and compute node metrics |
| `add_subtype_expression()` | Add per-subtype mean expression columns |
| `annotate_low_expressed_genes()` | Flag low-expressed genes by percentile or threshold |
| `run_anova_subtype_expression()` | ANOVA + Bonferroni + eta-squared across subtypes |
| `rank_subtype_expression()` | Rank expression across subtypes per gene |
| `add_hpa_annotations()` | Join raw HPA normal tissue and cancer expression fields |
| `derive_hpa_metrics()` | Compute derived HPA classifications and prioritisation flags |
| `get_hpa_tissues()` | List available tissues/cancer types in an HPA data frame |
| `filter_network_nodes()` | Filter by tissue specificity, low expression, ANOVA significance |
| `recalculate_network_metrics()` | Recompute centrality on filtered subnetwork |
| `run_go_enrichment()` | GO over-representation analysis per Louvain community |
| `calculate_subtype_hub_scores()` | Subtype-specific hub scores (centrality × expression specificity) |
| `plot_subtype_hub_scores()` | ComplexHeatmap of subtype hub scores |
| `plot_network()` | Draw network in R (static ggraph or interactive visNetwork) |
| `create_cytoscape_objects()` | Export node/edge TSVs and GraphML for Cytoscape |

---

## Dependencies

**Required:**

| Package | Source | Purpose |
|---|---|---|
| `igraph` | CRAN | Network construction and centrality |
| `dplyr` | CRAN | Data manipulation |
| `STRINGdb` | Bioconductor | STRING PPI data access |
| `biomaRt` | Bioconductor | Gene ID mapping |

**Optional (required for specific functions):**

| Package | Source | Required by |
|---|---|---|
| `clusterProfiler` | Bioconductor | `run_go_enrichment()` |
| `org.Hs.eg.db` | Bioconductor | `run_go_enrichment()` |
| `ComplexHeatmap` | Bioconductor | `plot_subtype_hub_scores()` |
| `circlize` | CRAN | `plot_subtype_hub_scores()` |
| `ggraph` | CRAN | `plot_network()` static mode |
| `ggplot2` | CRAN | `plot_network()` static mode |
| `ggrepel` | CRAN | `plot_network()` label repulsion |
| `visNetwork` | CRAN | `plot_network()` interactive mode |

---

## Bundled Data

- **`sjodahl_2017`** — Microarray expression matrix from the Sjodahl 2017 bladder cancer cohort (genes × samples).
- **`sjodahl_2017_subtypes`** — Named character vector of molecular subtype labels for each sample.
- **`Ensembl2Reactome`** — Human Ensembl gene to Reactome pathway mapping (downloaded from [reactome.org](https://reactome.org)).
- **`hpa_normal_data`** — Human Protein Atlas normal tissue IHC data.
- **`hpa_cancer_data`** — Human Protein Atlas cancer expression data.

---

## Citation

If you use NetWeave in your research, please cite:

> Author (2026). *NetWeave: Protein Interaction Network Analysis with Expression and Annotation Integration*. R package version 0.1.0.

---

## License

MIT © 2026
