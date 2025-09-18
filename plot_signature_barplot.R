library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggtext)

setwd(" ")
# 1️⃣ Lire la matrice SigProfiler
df <- read_tsv(" ")
colnames(df) <- trimws(colnames(df))

# 2️⃣ Renommer les colonnes
colnames(df)[2:3] <- c("SBS12", "SBS15")

# 3️⃣ Format long
df_long <- df %>%
  pivot_longer(
    cols = c(SBS12, SBS15),
    names_to = "Signature",
    values_to = "Value"
  )

# 4️⃣ Extraire substitution et trinucleotide
df_long <- df_long %>%
  mutate(
    SubRaw = str_extract(Type, "\\[.*\\]") %>% str_sub(2, -2),   # ex: "A>T"
    Trinuc  = str_replace_all(Type, "\\[|\\]", "") %>%           # G[A>T]C → GA>T C
              str_replace(">", "") %>%                                   # GATC (mais 4 lettres…)
              str_sub(1,1) %>% paste0(str_sub(Type,3,3), str_sub(Type,5,5))
  )

# Correction : reconstruire correctement le trinuc
df_long <- df_long %>%
  mutate(Trinuc = paste0(substr(Type,1,1),
                         substr(SubRaw,1,1),
                         substr(Type,nchar(Type),nchar(Type))))

# 5️⃣ Re-mapper en 6 classes pyrimidine
df_long <- df_long %>%
  mutate(SubClass = case_when(
    SubRaw %in% c("C>A","G>T") ~ "C>A",
    SubRaw %in% c("C>G","G>C") ~ "C>G",
    SubRaw %in% c("C>T","G>A") ~ "C>T",
    SubRaw %in% c("T>A","A>T") ~ "T>A",
    SubRaw %in% c("T>C","A>G") ~ "T>C",
    SubRaw %in% c("T>G","A>C") ~ "T>G",
    TRUE ~ NA_character_
  ))

# 6️⃣ Couleurs
colors <- c("C>A"="#56B4E9", "C>G"="black", "C>T"="red",
            "T>A"="grey50", "T>C"="#CCFF66", "T>G"="#F5C6C6")

# 7️⃣ Labels avec la base centrale colorée
df_long <- df_long %>%
  mutate(
    TrinucLabel = paste0(
      substr(Trinuc,1,1),
      "<span style='color:", colors[SubClass], "'>", substr(Trinuc,2,2), "</span>",
      substr(Trinuc,3,3)
    )
  )

# 8️⃣ Générer une figure par signature
for(sig in c("SBS12", "SBS15")) {
  p <- ggplot(df_long %>% filter(Signature == sig),
              aes(x = TrinucLabel, y = Value, fill = SubClass)) +
    geom_bar(stat="identity") +
    facet_wrap(~SubClass, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(
      axis.text.x = element_markdown(angle=90, vjust=0.5, hjust=1), # texte enrichi
      legend.position = "none"
    ) +
    labs(x = "Substitution trinucleotide", y = "Mutation proportion",
         title = paste("Signature -", sig))

  ggsave(paste0(sig, "_signature_plot.png"), plot = p, width=12, height=6, dpi=300)
}

