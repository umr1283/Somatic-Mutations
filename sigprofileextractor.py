import pandas as pd
from SigProfilerExtractor import sigpro as sig

if __name__ == "__main__":
  sig.sigProfilerExtractor(
    input_type="matrix",                     # car tu as déjà généré ta matrice SBS96
    output=" ",                      # dossier parent
    input_data="SBS96_matrix.tsv",           # ton fichier matrice
    reference_genome="GRCh38",               # cohérent avec ton FASTA
    opportunity_genome="GRCh38",             # idem
    context_type="96",                       # standard SBS96
    minimum_signatures=1,                    # borne inférieure testée
    maximum_signatures=5,                   # borne supérieure testée
    nmf_replicates=50,                       # plus tu mets, plus c’est robuste
  )
