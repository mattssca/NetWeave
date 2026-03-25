#Download reactome and STRING data
#load packages
library(dplyr)
library(biomaRt)

# Create directories if they don't exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/string_db", recursive = TRUE, showWarnings = FALSE)

#========================================
# Download Reactome Data
#========================================
message("Downloading Reactome data...")

download.file(
  url = "https://reactome.org/download/current/Ensembl2Reactome.txt",
  destfile = "data/raw/Ensembl2Reactome.txt"
)

#read the data
Ensembl2Reactome <- read.delim("data/raw/Ensembl2Reactome.txt", 
                               header = FALSE, 
                               stringsAsFactors = FALSE)

colnames(Ensembl2Reactome) <- c("gene_id", "pathway_id", "url", "pathway_name", "evidence", "species")

#filter for human pathways only
Ensembl2Reactome <- Ensembl2Reactome[Ensembl2Reactome$species == "Homo sapiens", ]

#save object
save(Ensembl2Reactome, file = "data/Ensembl2Reactome.Rdata")

message("Reactome data saved!")

#========================================
# Download STRING Database Files
#========================================
string_dir <- "data/string_db"

#download protein aliases
download.file(
  url = "https://stringdb-downloads.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz",
  destfile = file.path(string_dir, "9606.protein.aliases.v11.5.txt.gz"),
  mode = "wb"
)

#download protein info
download.file(
  url = "https://stringdb-downloads.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz",
  destfile = file.path(string_dir, "9606.protein.info.v11.5.txt.gz"),
  mode = "wb"
)

#download protein links
download.file(
  url = "https://stringdb-downloads.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz",
  destfile = file.path(string_dir, "9606.protein.links.v11.5.txt.gz"),
  mode = "wb"
)

#========================================
# Download HPA Data
#========================================
# normal tissue IHC data
download.file(
  url = "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip",
  destfile = "data/raw/normal_hpa_data.tsv.zip",
  mode = "wb"
)

unzip("data/raw/normal_hpa_data.tsv.zip", exdir = "data/raw/")

hpa_normal_data <- read.delim("data/raw/normal_ihc_data.tsv", 
                              header = TRUE, 
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

save(hpa_normal_data, file = "data/hpa_normal_data.Rdata")

# cancer data
download.file(
  url = "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip",
  destfile = "data/raw/cancer_hpa_data.tsv.zip",
  mode = "wb"
)

unzip("data/raw/cancer_hpa_data.tsv.zip", exdir = "data/raw/")

hpa_cancer_data <- read.delim("data/raw/cancer_data.tsv", 
                              header = TRUE, 
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

save(hpa_cancer_data, file = "data/hpa_cancer_data.Rdata")
