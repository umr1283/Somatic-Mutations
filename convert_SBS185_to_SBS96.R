
#!/usr/bin/env Rscript

# Chargement des librairies
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

# ================================
# Fonctions utilitaires
# ================================

# Complément ADN
comp <- c(A="T", T="A", C="G", G="C")

# Conversion SBS185 -> SBS96 (C/T centré)
to_SBS96 <- function(motif) {
  # motif du type X[Y>Z]Z
  m <- str_match(motif, "([ACGT])\\[([ACGT])>([ACGT])\\]([ACGT])")
  if (is.na(m[1,1])) return(NA)
  
  left  <- m[1,2]
  ref   <- m[1,3]
  alt   <- m[1,4]
  right <- m[1,5]
  
  if (ref %in% c("C","T")) {
    return(paste0(left,"[",ref,">",alt,"]",right))
  } else {
    return(paste0(comp[[right]],"[",comp[[ref]],">",comp[[alt]],"]",comp[[left]]))
  }
}

# ================================
# Lecture fichier
# ================================
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_SBS185_to_SBS96.R input.txt output.txt")
}

input_file  <- args[1]
output_file <- args[2]

cat("Lecture du fichier :", input_file, "\n")
df <- read_tsv(input_file, col_types = cols())

# ================================
# Conversion en SBS96
# ================================
cat("Conversion des motifs SBS185 vers SBS96...\n")
df <- df %>%
  mutate(SBS96 = sapply(MutationType, to_SBS96))

# ================================
# Agrégation
# ================================
cat("Agrégation des lignes par motif SBS96...\n")
df_SBS96 <- df %>%
  group_by(SBS96) %>%
  summarise(across(-MutationType, sum, na.rm=TRUE)) %>%
  ungroup()

# Vérif dimensions
cat("Nombre de motifs obtenus :", nrow(df_SBS96), "\n")

# ================================
# Normalisation colonnes (somme = 1)
# ================================
cat("Normalisation des proportions par échantillon...\n")
df_SBS96_norm <- df_SBS96 %>%
  mutate(across(-SBS96, ~ .x / sum(.x)))

# ================================
# Sauvegarde
# ================================
write_tsv(df_SBS96_norm, output_file)
cat("Tableau SBS96 écrit dans :", output_file, "\n")
