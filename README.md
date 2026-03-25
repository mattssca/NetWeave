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
devtools::install_github("https://github.com/mattssca/NetWeave")
```

---

## Pipeline Overview

NetWeave follows a linear, modular pipeline. Each function appends new columns to `node_metrics` and passes it forward — making it easy to stop, inspect, or branch at any step.

```
[run_netweave()  ──────────────────────────────────────────────────────────────]
                                                                               |
get_pathway_genes()             → pathway gene vector                         |
     │                                                                         |
     ▼                                                                         |
create_and_expand_network()     → edge_data + node_metrics                    |
     │                                                                         |
     ▼                                                                         |
add_subtype_expression()        → mean_expr_<subtype> columns                 |
     │                                                                         |
     ▼                                                                         |
annotate_low_expressed_genes()  → is_low_expression                           |
     │                                                                         |
     ▼                                                                         |
run_anova_subtype_expression()  → anova_p_value, anova_sig, ...               |
     │                                                                         |
     ▼                                                                         |
rank_subtype_expression()       → rank_<subtype> columns                      |
     │                                                                         |
     ▼                                                                         |
add_hpa_annotations()           → hpa_normal_*, hpa_cancer_* columns          |
     │                                                                         |
     ▼                                                                         |
derive_hpa_metrics()            → hpa_expression_pattern, hpa_passes_filter, ...|
     │                                                                         |
     ▼                                                                         |
filter_network_nodes()          → filtered node_metrics                       |
     │                                                                         |
     ▼                                                                         |
recalculate_network_metrics()   → updated degree, betweenness, community, ... |
     │                                                                   [─────]
     ├──▶ calculate_subtype_hub_scores()   → <subtype>_hub_score columns
     │         └──▶ plot_subtype_hub_scores()
     ├──▶ plot_network()                   → ggraph / visNetwork
     └──▶ create_cytoscape_objects()       → TSV + GraphML export
```

---

## Quick Start

### Step-by-step

```r
library(NetWeave)

# Load bundled example data
data(sjodahl_2017)          # expression matrix
data(sjodahl_2017_subtypes) # named subtype vector
data(usq_expr)              # independent cohort for low-expression thresholding
data(hpa_normal_data)
data(hpa_cancer_data)

# 1. Get pathway seed genes
erbb2_pathway <- get_pathway_genes(
  this_pathway = "Signaling by ERBB2",
  verbose      = TRUE
)

# 2. Build and expand the STRING network
erbb2_network <- create_and_expand_network(
  expr_data              = sjodahl_2017,
  seed_gene              = "ERBB2",
  pathway_genes          = erbb2_pathway,
  max_added_genes        = 20,
  string_score_threshold = 900,
  string_data_dir        = "data/string_db/"
)

# 3. Add per-subtype mean expression
annotated_node_metrics <- add_subtype_expression(
  this_subtype_vector = sjodahl_2017_subtypes,
  expr_data           = sjodahl_2017,
  node_metrics        = erbb2_network$node_metrics,
  subtype_order       = c("Uro", "GU", "BaSq", "Mes", "ScNE")
)

# 4. Flag low-expressed genes (uses independent cohort)
annotated_node_metrics <- annotate_low_expressed_genes(
  node_metrics       = annotated_node_metrics,
  expr_data          = usq_expr,
  low_expr_threshold = 0.25
)

# 5. ANOVA for subtype-differential expression
annotated_node_metrics <- run_anova_subtype_expression(
  node_metrics        = annotated_node_metrics,
  expr_data           = sjodahl_2017,
  this_subtype_vector = sjodahl_2017_subtypes,
  sig_threshold       = 0.05
)

# 6. Rank expression across subtypes
annotated_node_metrics <- rank_subtype_expression(
  node_metrics = annotated_node_metrics
)

# 7. Add HPA annotations
annotated_node_metrics <- add_hpa_annotations(
  node_metrics       = annotated_node_metrics,
  normal_hpa         = hpa_normal_data,
  cancer_hpa         = hpa_cancer_data,
  this_tissue_normal = "Urinary bladder",
  this_tissue_cancer = "urothelial cancer"
)

# 8. Derive HPA metrics
annotated_node_metrics <- derive_hpa_metrics(
  node_metrics = annotated_node_metrics
)

# 9. Filter nodes
annotated_node_metrics <- filter_network_nodes(
  node_metrics         = annotated_node_metrics,
  tissue_filter        = "either",
  remove_low_expressed = TRUE,
  filter_anova_sig     = TRUE
)

# 10. Recalculate metrics on the filtered subnetwork
annotated_node_metrics <- recalculate_network_metrics(
  node_metrics = annotated_node_metrics,
  edge_data    = erbb2_network$edge_data
)

# 11. Subtype hub scores
annotated_node_metrics <- calculate_subtype_hub_scores(
  node_metrics = annotated_node_metrics
)

# Visualise hub scores
hm_up   <- plot_subtype_hub_scores(annotated_node_metrics, hub_direction = "up")
hm_down <- plot_subtype_hub_scores(annotated_node_metrics, hub_direction = "down")
draw(hm_up)
draw(hm_down)

# Draw network
plot_network(
  node_metrics = annotated_node_metrics,
  edge_data    = erbb2_network$edge_data,
  color_by     = "origin",
  label_nodes  = "all"
)

# Export to Cytoscape
create_cytoscape_objects(
  node_metrics = annotated_node_metrics,
  edge_data    = erbb2_network$edge_data,
  out_dir      = "output/",
  file_prefix  = "erbb2_network"
)
```

### Or use the wrapper

```r
result <- run_netweave(
  seed_gene          = "ERBB2",
  pathway_name       = "Signaling by ERBB2",
  expr_data          = sjodahl_2017,
  subtype_vector     = sjodahl_2017_subtypes,
  low_expr_data      = usq_expr,
  normal_hpa         = hpa_normal_data,
  cancer_hpa         = hpa_cancer_data,
  this_tissue_normal = "Urinary bladder",
  this_tissue_cancer = "urothelial cancer",
  subtype_order      = c("Uro", "GU", "BaSq", "Mes", "ScNE"),
  string_data_dir    = "data/string_db/"
)

# result$node_metrics  — final annotated node table
# result$edge_data     — network edge list
# result$pathway_genes — pathway gene vector from step 1
```

---

## Function Reference

| Function | Description |
|---|---|
| `run_netweave()` | Run the complete pipeline in a single call |
| `get_pathway_genes()` | Retrieve genes for a Reactome pathway or seed gene |
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
- **`usq_expr`** — Independent expression matrix used for low-expression thresholding.
- **`Ensembl2Reactome`** — Human Ensembl gene to Reactome pathway mapping (downloaded from [reactome.org](https://reactome.org)).
- **`hpa_normal_data`** — Human Protein Atlas normal tissue IHC data.
- **`hpa_cancer_data`** — Human Protein Atlas cancer expression data.

---

## Citation

If you use NetWeave in your research, please cite:

> Carl-Adam Mattsson (2026). *NetWeave: Protein Interaction Network Analysis with Expression and Annotation Integration*. R package version 0.1.0.

---

## License

MIT © 2026
