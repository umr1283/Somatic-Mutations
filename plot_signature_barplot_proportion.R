library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# -------------------------------
# 1️⃣ Lecture du fichier
# -------------------------------
df <- read_tsv(" ")

# -------------------------------
# 2️⃣ Calcul de la proportion CH185A pour trier
# -------------------------------
df <- df %>%
  rowwise() %>%
  mutate(PropA = CH185A / sum(c_across(CH185A:CH185B))) %>%
  ungroup()

# -------------------------------
# 3️⃣ Passage en format long pour ggplot
# -------------------------------
df_long <- df %>%
  select(Samples, CH185A, CH185B) %>%
  pivot_longer(cols = CH185A:CH185B,
               names_to = "Signature",
               values_to = "Variants") %>% mutate(Signature = recode(Signature,
                            "CH185A" = "A",
                            "CH185B" = "B"))

# -------------------------------
# 4️⃣ Définir l'ordre des patients par proportion CH185A décroissante
# -------------------------------
patient_order <- df %>%
  arrange(desc(PropA)) %>%
  pull(Samples)

df_long$Samples <- factor(df_long$Samples, levels = patient_order)

# -------------------------------
# 5️⃣ Création du plot empilé
# -------------------------------
ggplot(df_long, aes(x = Samples, y = Variants, fill = Signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("PA"="#56B4E9", "B"="#E69F00")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6),
    legend.position = "top"
  ) +
  labs(
    x = "Patient",
    y = "Nombre de variants",
    fill = "Signature",
    title = "Nombre de variants par signature et par patient"
  )

# -------------------------------
# 6️⃣ Sauvegarde
# -------------------------------
ggsave("variants_per_patient_stacked.png", width = 12, height = 6, dpi = 300)