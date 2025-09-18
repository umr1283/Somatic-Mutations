#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# ================================
# Canonical SBS96 order
# ================================
bases <- c("A","C","G","T")
contexts <- as.vector(outer(bases, bases, paste0))
subst_order <- c("C>A","C>G","C>T","T>A","T>C","T>G")

canonical_order <- unlist(lapply(subst_order, function(subst) {
  ref <- substr(subst,1,1)
  alt <- substr(subst,3,3)
  sapply(contexts, function(ctx) {
    paste0(substr(ctx,1,1),"[",ref,">",alt,"]",substr(ctx,2,2))
  })
}))

# ================================
# Parsing motif
# ================================
parse_mut <- function(motif) {
  m <- str_match(motif, "([ACGT])\\[([ACGT])>([ACGT])\\]([ACGT])")
  if (is.na(m[1,1])) return(c(NA, NA, NA))
  left  <- m[1,2]; ref <- m[1,3]; alt <- m[1,4]; right <- m[1,5]
  subst <- paste0(ref,">",alt)
  tri   <- paste0(left,ref,right)
  return(c(subst, tri, motif))
}

# ================================
# Format numérique (9 décimales, expo avec E majuscule)
# ================================
fmt_num <- function(x) {
  ifelse(is.na(x), "",
         toupper(formatC(x, format="e", digits=9)))
}

# ================================
# Script principal
# ================================
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert_SBS96.R input.txt output.txt")
}

input_file  <- args[1]
output_file <- args[2]

cat("Lecture du fichier :", input_file, "\n")
df <- read_tsv(input_file, col_types = cols())

# extraction infos
parsed <- t(sapply(df$Type, parse_mut))
df_out <- df %>%
  mutate(`Substitution Type`   = parsed[,1],
         Trinucleotide         = parsed[,2],
         `Somatic Mutation Type` = parsed[,3]) %>%
  select(`Substitution Type`, Trinucleotide, `Somatic Mutation Type`, everything(), -Type)

# complétion des 96 catégories manquantes
df_out <- full_join(
  tibble(`Somatic Mutation Type` = canonical_order),
  df_out,
  by = "Somatic Mutation Type"
) %>%
  mutate(
    `Substitution Type` = ifelse(is.na(`Substitution Type`),
                                 str_extract(`Somatic Mutation Type`, "(C|T)>[ACGT]"),
                                 `Substitution Type`),
    Trinucleotide = ifelse(is.na(Trinucleotide),
                           str_replace_all(`Somatic Mutation Type`, "\\[|\\]|>[ACGT]", ""),
                           Trinucleotide)
  )

# tri selon l’ordre SBS96 canonique
df_out <- df_out %>%
  slice(match(canonical_order, `Somatic Mutation Type`))

# ============
# Filtrage signatures 2→30 + 96
# ============
valid_sigs <- c(paste0("SBS", 1:27), "SBS29", "SBS30", "SBS96")

keep_cols <- c("Substitution Type", "Trinucleotide", "Somatic Mutation Type",
               intersect(valid_sigs, colnames(df_out)))

df_out <- df_out %>% select(all_of(keep_cols))

# mise en forme numérique
df_out <- df_out %>%
  mutate(across(-c(`Substitution Type`, Trinucleotide, `Somatic Mutation Type`),
                ~ fmt_num(.x)))

# écriture
write_tsv(df_out, output_file)
cat("Tableau SBS96 écrit dans :", output_file, "\n")

