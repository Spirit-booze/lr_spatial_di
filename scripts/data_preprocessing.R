library(tidyverse) # для общего парсинга
library(biomaRt) # прицельно для работы с идентификаторами генов

colnames = c('read_id', 'gene_id_ensembl', 'cell_type', 'barcode', 'transcript_id_ensembl')
df_PFC_Hippo = readr::read_tsv('raw_data/P7.PFC+Hipp.tsv.gz', col_names = colnames)
#df_Curio_all = readr::read_tsv('raw_data/Curio.all_samples.allinfo.tsv.gz', col_names = colnames)

# эта регулярка сводит с ума
reads <- df_PFC_Hippo %>%
  mutate(gene_id_clean = sub("\\.\\d+$", "", gene_id_ensembl))

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

unique_genes <- unique(reads$gene_id_clean)

# mgi символ вместо hgnc, потому что МЫШЬ
annot <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique_genes,
  mart       = mart
)

# на всякий случай без дублей
annot = annot %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

# чуть витиевато, чтобы потом толком отлавливать NA
reads = reads %>% 
  left_join(annot, by=c("gene_id_clean" = "ensembl_gene_id")) %>% 
  mutate(
    gene_symbol = ifelse(
      mgi_symbol == "" | is.na(mgi_symbol),
      gene_id_clean,
      mgi_symbol
    )
  )

sum(grepl("^ENSMUSG", reads$gene_symbol))
# ну то есть 3622 случая не получили идентификатора MGI. Бывает.

# теперь преобразуем это в матрицу каунтов
counts_long = reads %>% 
  group_by(barcode, gene_symbol) %>% 
  summarise(count=n(), .groups="drop")

count_wise = counts_long %>% 
  pivot_wider(
    names_from = barcode,
    values_from = count,
    values_fill = 0,
  )

