# CALCOLO FOO, RRA e pulizia dataset
library(tidyverse)

slo <- read_delim("Slovenia wolves FILE DEFINITIVO.csv", delim = ";", 
                   locale = locale(encoding = "UTF-8"))
cro <- read_delim("Croatian wolves FILE DEFINITIVO.csv", delim = ";", 
                   locale = locale(encoding = "UTF-8"))

# Converto da formato largo (specie x campione) a formato lungo
slo_long <- slo %>%
  select(-rank) %>%
  pivot_longer(-scientific_name, names_to = "Sample_ID", values_to = "reads") %>%
  mutate(Popolazione = "Slovenia")

cro_long <- cro %>%
  select(-rank) %>%
  pivot_longer(-scientific_name, names_to = "Sample_ID", values_to = "reads") %>%
  mutate(Popolazione = "Croazia")

grezzi <- bind_rows(slo_long, cro_long)

nrow(grezzi)
n_distinct(grezzi$Sample_ID)
grezzi <- grezzi %>%
  group_by(Sample_ID) %>%
  mutate(totale_grezzo = sum(reads)) %>%
  ungroup() %>%
  mutate(pct_grezza = if_else(totale_grezzo > 0, reads / totale_grezzo * 100, 0))

# Quante osservazioni specie-campione vengono azzerate dalla soglia dell'1%?
# (reads>0 ma sotto la soglia ecologica)
sotto_soglia <- grezzi %>% filter(reads > 0 & pct_grezza < 1)
nrow(sotto_soglia)
print(sotto_soglia)

# Applico la soglia: azzero le read sotto l'1%, poi ricalcolo il totale
# "pulito" e la percentuale finale su quel nuovo totale
grezzi <- grezzi %>%
  mutate(reads_filtrate = if_else(pct_grezza < 1, 0, reads)) %>%
  group_by(Sample_ID) %>%
  mutate(totale_filtrato = sum(reads_filtrate)) %>%
  ungroup() %>%
  mutate(RRA_percentuale = if_else(totale_filtrato > 0, reads_filtrate / totale_filtrato * 100, 0))

# Verifica: per ogni campione con almeno una read valida, l'RRA deve sommare a 100
controllo_somme <- grezzi %>%
  group_by(Sample_ID) %>%
  summarise(somma_RRA = sum(RRA_percentuale))
summary(controllo_somme$somma_RRA)
rra_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(RRA_percentuale = mean(RRA_percentuale), .groups = "drop")

print(rra_pop, n = Inf)
foo_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(FOO_percentuale = mean(RRA_percentuale > 0) * 100, .groups = "drop")

print(foo_pop, n = Inf)
# Il conteggio assoluto delle presenze (mancava, serve per la colonna "Presenza_Assoluta")
presenza_pop <- grezzi %>%
  group_by(Popolazione, scientific_name) %>%
  summarise(Presenza_Assoluta = sum(RRA_percentuale > 0), .groups = "drop")

# Unisco le tre tabelle (presenza + FOO + RRA) in un'unica tabella finale
foo_rra_finale <- presenza_pop %>%
  left_join(foo_pop, by = c("Popolazione", "scientific_name")) %>%
  left_join(rra_pop, by = c("Popolazione", "scientific_name")) %>%
  filter(Presenza_Assoluta > 0) %>%
  rename(Preda = scientific_name) %>%
  mutate(
    FOO_percentuale = round(FOO_percentuale, 2),
    RRA_percentuale = round(RRA_percentuale, 2)
  ) %>%
  arrange(Popolazione, desc(RRA_percentuale))

print(foo_rra_finale, n = Inf)

# Esportazione — sovrascrive il file rotto con quello corretto
write.csv2(foo_rra_finale, "FOO_RRA_Slovenia_and_Croatia.csv", row.names = FALSE)

# calcolo indice di pianka e di levins
rra_matrice <- rra_pop %>%
  pivot_wider(names_from = scientific_name, values_from = RRA_percentuale, values_fill = 0)

levins_finale <- rra_matrice %>%
  rowwise() %>%
  mutate(
    B = 1 / sum(c_across(-Popolazione)^2 / 10000),
    B_A = (B - 1) / (10 - 1)
  ) %>%
  select(Popolazione, B_A)
print(levins_finale)

p_slo <- as.numeric(rra_matrice[rra_matrice$Popolazione=="Slovenia", -1]) / 100
p_cro <- as.numeric(rra_matrice[rra_matrice$Popolazione=="Croazia", -1]) / 100
pianka_finale <- sum(p_slo * p_cro) / sqrt(sum(p_slo^2) * sum(p_cro^2))
print(pianka_finale)




# =========================================================================
# PIPELINE INTEGRATA DIETA LUPO: SLOVENIA VS CROAZIA
# PARTE 1: Grafici di composizione della dieta (FOO% e RRA% - Stile Jura)
# PARTE 2: Creazione Community Matrix (RRA% per singolo campione)
# PARTE 3: Statistica Multivariata (NMDS con vettori envfit, PERMANOVA, SIMPER)
# =========================================================================

# Caricamento delle librerie fondamentali per l'intera pipeline
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(vegan)

# =========================================================================
# PARTE 1: CREAZIONE GRAFICI PUBBLICABILI (FOO% E RRA%)
# =========================================================================

# 1. CARICAMENTO DATI INDICI DIETETICI
df_indici <- read_csv2("FOO&RRA Slovenia and Croatia.csv", show_col_types = FALSE)

# Trasformiamo eventuali valori NA in 0
df_indici[is.na(df_indici)] <- 0

# 2. CATEGORIZZAZIONE E COLORI (TUTTI I TAXA INCLUSI)
# Creiamo una colonna "Categoria" includendo anche le famiglie e sottofamiglie
df_indici <- df_indici %>%
  mutate(Categoria = case_when(
    Preda %in% c("Capreolus capreolus", "Cervus elaphus", "Cervinae", "Sus scrofa", "Rupicapra rupicapra") ~ "wild ungulates",
    Preda %in% c("Ovis aries", "Caprinae", "Ovis", "Capra") ~ "sheep & goats",
    Preda %in% c("Bos") ~ "cattle",
    TRUE ~ "small fauna & other" # Comprende Lepus e altra microfauna
  ))

# Definiamo la palette di colori ESATTA ispirata allo screenshot del paper Jura
colori_jura <- c(
  "wild ungulates" = "#007A4D",       # Verde scuro
  "sheep & goats" = "#E2D19A",        # Beige
  "cattle" = "#B3987E",               # Marroncino chiaro
  "small fauna & other" = "#95C78F"   # Verde chiaro
)

# Fissiamo l'ordine della legenda come nel paper
df_indici$Categoria <- factor(df_indici$Categoria,
                              levels = c("wild ungulates", "sheep & goats", "small fauna & other", "cattle"))

# 3. FUNZIONE GENERATRICE DEL GRAFICO A 2 PANNELLI
# Creiamo una funzione che prende i dati e il nome della nazione e genera il plot
crea_grafico_nazione <- function(dati, nazione) {
  
  # Filtriamo solo i dati della nazione richiesta
  df_nazione <- dati %>% filter(Popolazione == nazione)
  
  # Rimuoviamo le prede che in questa nazione sono a 0 assoluto per non avere barre vuote
  df_nazione <- df_nazione %>% filter(FOO_percentuale > 0 | RRA_percentuale > 0)
  
  # Ordiniamo le prede dalla più consumata alla meno consumata basandoci sulla RRA%
  ordine_prede <- df_nazione %>%
    arrange(RRA_percentuale) %>%
    pull(Preda)
  
  df_nazione$Preda <- factor(df_nazione$Preda, levels = ordine_prede)
  
  # Allunghiamo il dataset per creare i due pannelli (solo FOO% e RRA% come da indicazioni di Marta)
  df_long <- df_nazione %>%
    pivot_longer(
      cols = c(FOO_percentuale, RRA_percentuale),
      names_to = "Metrica",
      values_to = "Valore"
    ) %>%
    mutate(Metrica = recode(Metrica,
                            "FOO_percentuale" = "Frequency of Occurrence",
                            "RRA_percentuale" = "Relative Read Abundance"))
  
  # Fissiamo l'ordine dei pannelli da sinistra a destra
  df_long$Metrica <- factor(df_long$Metrica, levels = c("Frequency of Occurrence", "Relative Read Abundance"))
  
  # Calcoliamo il limite massimo dell'asse X per dare spazio ai numeri alla fine della barra
  max_x <- max(df_long$Valore) * 1.15
  
  # Creazione del Plot
  p <- ggplot(df_long, aes(x = Valore, y = Preda, fill = Categoria)) +
    geom_bar(stat = "identity", width = 0.8) +
    
    # Etichetta numerica alla fine della barra con 1 decimale
    geom_text(aes(label = round(Valore, 1)), hjust = -0.2, size = 3.5, color = "black") +
    
    # Creiamo i due pannelli separati
    facet_wrap(~ Metrica, scales = "free_x") +
    
    # Applichiamo i colori scelti
    scale_fill_manual(values = colori_jura) +
    
    # Espandiamo l'asse X in base al massimo calcolato
    scale_x_continuous(limits = c(0, max_x), expand = c(0, 0)) +
    
    # Titoli ed etichette
    labs(
      title = paste("Diet Composition -", nazione),
      x = "Percentage (%)",
      y = NULL,
      fill = NULL
    ) +
    
    # Stile e formattazione (in linea con lo stile del paper)
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "italic", size = 12, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.5, "cm")
    )
  
  return(p)
}

# 4. GENERAZIONE E VISUALIZZAZIONE GRAFICI PARTE 1
plot_slo <- crea_grafico_nazione(df_indici, "Slovenia")
plot_cro <- crea_grafico_nazione(df_indici, "Croazia")

# Mostra i grafici direttamente su RStudio nel pannello Plots
print(plot_slo)
print(plot_cro)

# Salvataggio opzionale dei file PNG (rimuovere il cancelletto iniziale per salvare)
# ggsave("Grafico_Dieta_Slovenia_Jura.png", plot = plot_slo, width = 10, height = 5, dpi = 300)
# ggsave("Grafico_Dieta_Croazia_Jura.png", plot = plot_cro, width = 10, height = 5, dpi = 300)


# =========================================================================
# PARTE 2: CREAZIONE DELLA COMMUNITY MATRIX PER STATISTICA MULTIVARIATA
# =========================================================================

# 1. CARICAMENTO DEI DATASET PURIFICATI PER CAMPIONE
file_slo <- "Slovenia wolves FILE DEFINITIVO.csv"
file_cro <- "Croatian wolves FILE DEFINITIVO.csv"

df_slo <- read_csv2(file_slo, show_col_types = FALSE)
df_cro <- read_csv2(file_cro, show_col_types = FALSE)

# 2. FUNZIONE CORAZZATA PER TRASPORRE E PREPARARE I DATI
prepara_matrice <- function(df, nome_popolazione) {
  
  # Identifichiamo dinamicamente la colonna delle specie scientifiche
  col_specie <- names(df)[grepl("scientific_name", names(df), ignore.case = TRUE)][1]
  
  # Rinominiamo la colonna specie in "Taxon" per uniformità
  names(df)[names(df) == col_specie] <- "Taxon"
  
  # Rimuoviamo in blocco tutte le possibili colonne di metadati tassonomici
  df <- df %>% select(-any_of(c("rank", "final rank", "species_list")))
  
  # Trasformiamo la tabella da formati orizzontali a verticali per il calcolo
  df_long <- df %>% 
    pivot_longer(cols = -Taxon, names_to = "Sample_ID", values_to = "Reads")
  
  df_wide <- df_long %>% 
    pivot_wider(names_from = Taxon, values_from = Reads, values_fill = 0)
  
  # Aggiungiamo la colonna fondamentale della popolazione di appartenenza
  df_wide$Popolazione <- nome_popolazione
  
  return(df_wide)
}

# 3. APPLICAZIONE DELLA FUNZIONE E UNIONE DELLE MATRICI
mat_slo <- prepara_matrice(df_slo, "Slovenia")
mat_cro <- prepara_matrice(df_cro, "Croazia")

# Uniamo le due nazioni in un unico dataframe
community_matrix <- bind_rows(mat_slo, mat_cro)
community_matrix[is.na(community_matrix)] <- 0

# Identifichiamo le colonne delle prede (escludendo Sample_ID e Popolazione)
colonne_prede <- setdiff(names(community_matrix), c("Sample_ID", "Popolazione"))

# 4. STANDARDIZZAZIONE IN RRA% PER OGNI CAMPIONE
# Trasformiamo i conteggi di reads assolute in percentuali di RRA per riga
community_matrix[colonne_prede] <- t(apply(community_matrix[colonne_prede], 1, function(x) {
  totale <- sum(x)
  if(totale > 0) { 
    return((x / totale) * 100) 
  } else { 
    return(x) 
  }
}))

# Riordiniamo le colonne mettendo ID e Popolazione all'inizio
community_matrix <- community_matrix %>%
  select(Sample_ID, Popolazione, everything())

# Esportiamo il file CSV definitivo della Community Matrix
write.csv2(community_matrix, "Community_Matrix_SloCro.csv", row.names = FALSE)

print("Community Matrix creata con successo!")
print(paste("Totale campioni (righe):", nrow(community_matrix)))
print(paste("Totale specie preda (colonne):", length(colonne_prede)))


# =========================================================================
# =========================================================================
# PARTE 3 (CORRETTA): STATISTICA MULTIVARIATA — PERMANOVA e SIMPER
# Ricalcolati su dati Hellinger-trasformati per essere coerenti con la PCoA
# ufficiale e con il Betadisper già corretto. La versione precedente di
# questo script (con NMDS) usava i dati RRA grezzi — qui invece PERMANOVA
# e SIMPER lavorano sulla stessa identica matrice usata per PCoA/Betadisper.
# =========================================================================
 
library(tidyverse)
library(vegan)
 
# ==============================================================================
# 1. CARICAMENTO DELLA COMMUNITY MATRIX
# ==============================================================================
df <- read_csv2("Community_Matrix_SloCro.csv", show_col_types = FALSE)
 
# Pulizia di sicurezza: rimuoviamo eventuali colonne fantasma o artefatti di join
df <- df %>% select(-any_of("NA"))
 
# NOTA TASSONOMICA (risolta): la distinzione "Cervinae" vs "Cervus elaphus" per
# la Slovenia era dovuta a un rank non aggiornato nei file grezzi. Confermato
# tramite BLASTn che si tratta di Cervus elaphus a livello di specie in
# entrambe le popolazioni — nessuna azione necessaria, dato già coerente.
 
# ==============================================================================
# 2. PREPARAZIONE DEI DATI PER IL PACCHETTO VEGAN
# ==============================================================================
metadati <- df %>% select(Sample_ID, Popolazione)
matrice_prede <- df %>% select(-Sample_ID, -Popolazione)
 
# Forziamo a numerico per evitare conflitti nel calcolo delle distanze
matrice_prede <- as.data.frame(lapply(matrice_prede, as.numeric))
 
# ==============================================================================
# 3. TRASFORMAZIONE DI HELLINGER
# Stessa trasformazione già usata per la PCoA e il Betadisper — necessaria
# per coerenza metodologica tra tutte le analisi multivariate della tesi.
# ==============================================================================
matrice_hellinger <- decostand(matrice_prede, method = "hellinger")
 
# ==============================================================================
# 4. PERMANOVA (su dati Hellinger-trasformati)
# ==============================================================================
set.seed(123)
permanova_risultato <- adonis2(matrice_hellinger ~ Popolazione, data = metadati,
                                 method = "bray", permutations = 9999)
 
print(" ")
print("==================================================")
print("        RISULTATO PERMANOVA (Hellinger)          ")
print("==================================================")
print(permanova_risultato)
# NB: il valore chiave è nella colonna Pr(>F) per il p-value e in R2 per la
# percentuale di varianza spiegata
 
# ==============================================================================
# 5. SIMPER (su dati Hellinger-trasformati)
# ==============================================================================
set.seed(123)
simper_risultato <- simper(matrice_hellinger, group = metadati$Popolazione, permutations = 999)
 
print(" ")
print("==================================================")
print("          RISULTATO SIMPER (Hellinger)            ")
print("==================================================")
print(summary(simper_risultato))
#plotting di SIMPER
library(ggplot2)
library(dplyr)

# Estraggo la tabella completa direttamente dall'oggetto simper già calcolato
# (fonte unica, niente numeri ricopiati a mano)
simper_summary <- summary(simper_risultato)[[1]]
simper_summary$Taxon <- rownames(simper_summary)

simper_summary <- simper_summary %>%
  mutate(
    Taxon = gsub("\\.", " ", Taxon),  # ripristina lo spazio nei nomi (es. "Capreolus capreolus")
    contributo_pct = round(average / sum(average) * 100, 1),
    sig_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    significativo = p < 0.05
  ) %>%
  arrange(desc(contributo_pct))

# Grafico — TUTTI i taxa, non solo i primi 5, colorati per significatività
plot_simper <- ggplot(simper_summary, aes(x = reorder(Taxon, contributo_pct), 
                                            y = contributo_pct, fill = significativo)) +
  geom_col() +
  geom_text(aes(label = paste0(contributo_pct, "% ", sig_label)),
            hjust = -0.1, fontface = "bold", size = 3.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#B22222", "FALSE" = "#003f5c"),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = NULL) +
  expand_limits(y = max(simper_summary$contributo_pct) * 1.15) +
  labs(
    title = "SIMPER Analysis: Species Contributions to Dietary Divergence",
    subtitle = "Slovenia vs Croatia | Bray-Curtis on Hellinger-transformed RRA%\n(* p<0.05, ** p<0.01, *** p<0.001)",
    x = NULL, y = "Contribution to overall Bray-Curtis dissimilarity (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

print(plot_simper)
ggsave("Grafico_SIMPER_corretto.png", plot_simper, width = 10, height = 7, dpi = 300)

# ==============================================================================
# ANALISI DELLA DISPERSIONE MULTIVARIATA (BETADISPER) — versione corretta
# Coerente con la trasformazione di Hellinger già usata nella PCoA/PERMANOVA
# ufficiale della tesi (Sezione 2.5). Autosufficiente: nessuna dipendenza da
# script precedenti o oggetti già in memoria.
# ==============================================================================

library(tidyverse)
library(vegan)

# ==============================================================================
# 1. CARICAMENTO DATI
# ==============================================================================
dataset_lupo <- read_csv2("Community_Matrix_SloCro.csv",
                           locale = locale(decimal_mark = ","),
                           show_col_types = FALSE)

dataset_lupo$Popolazione <- as.factor(dataset_lupo$Popolazione)

taxa_cols <- c("Capreolus capreolus", "Caprinae", "Cervus elaphus", "Ovis aries",
               "Rupicapra rupicapra", "Sus scrofa", "Bos", "Capra", "Lepus", "Ovis")

# Rete di sicurezza per il parsing numerico (ridondante se read_csv2 ha già
# interpretato correttamente il separatore decimale, ma innocua)
for (col in taxa_cols) {
  dataset_lupo[[col]] <- as.numeric(gsub(",", ".", dataset_lupo[[col]]))
}

matrice_prede <- dataset_lupo %>% select(all_of(taxa_cols))

cat("N campioni:", nrow(matrice_prede), "\n")
cat("N taxa:", ncol(matrice_prede), "\n")
cat("Popolazioni:", paste(levels(dataset_lupo$Popolazione), collapse = ", "), "\n")

# ==============================================================================
# 2. TRASFORMAZIONE DI HELLINGER
# Stessa trasformazione della PCoA/PERMANOVA ufficiale: Betadisper deve validare
# le assunzioni dello STESSO test, quindi va calcolato sulla stessa matrice
# trasformata, non sui dati grezzi.
# ==============================================================================
matrice_hellinger <- decostand(matrice_prede, method = "hellinger")

# ==============================================================================
# 3. MATRICE DI DISTANZA E BETADISPER
# ==============================================================================
dist_matrix <- vegdist(matrice_hellinger, method = "bray")

dispersion_mod <- betadisper(dist_matrix, group = dataset_lupo$Popolazione)

# ==============================================================================
# 4. TEST DI PERMUTAZIONE
# ==============================================================================
set.seed(123)
cat("\n--- Risultato Test di Permutazione per Betadisper (Hellinger) ---\n")
test_dispersion <- permutest(dispersion_mod, permutations = 999)
print(test_dispersion)

# ==============================================================================
# 5. DISTANZA MEDIA DAI CENTROIDI (ampiezza di nicchia multivariata)
# ==============================================================================
cat("\n--- Distanza media dai centroidi (più alta = dieta più variabile) ---\n")
print(dispersion_mod$group.distances)

# ==============================================================================
# 6. GRAFICO
# ==============================================================================
png("Betadisper_Plot.png", width = 2000, height = 1600, res = 300)
plot(dispersion_mod, hull = FALSE, ellipse = TRUE,
     main = "Dispersione Multivariata della Dieta (Hellinger + Bray-Curtis)",
     col = c("#009E73", "#CC79A7"),
     lwd = 2, seg.col = "gray80", seg.lwd = 0.5)
legend("topleft", legend = levels(dataset_lupo$Popolazione),
       col = c("#009E73", "#CC79A7"), pch = 16, bty = "n", cex = 1.2)
dev.off()

cat("\nGrafico esportato come 'Betadisper_Plot.png'\n")



# ==============================================================================
# RICOSTRUZIONE DEFINITIVA DEL GRAFICO PCoA
# Bray-Curtis su dati Hellinger-trasformati, stessa metodologia già validata
# per RDA e Betadisper. Stesse 6 categorie di preda aggregate (Cervo, Capriolo,
# Cinghiale, Caprinae, Camoscio, Domestico) per coerenza visiva con la RDA.
# Autosufficiente: nessuna dipendenza da script precedenti.
# ==============================================================================

library(tidyverse)
library(vegan)
library(ggrepel)

# ==============================================================================
# 1. CARICAMENTO DATI E AGGREGAZIONE SPECIE DOMESTICHE
# ==============================================================================
df <- read_csv2("Community_Matrix_SloCro.csv", show_col_types = FALSE)

df <- df %>%
  mutate(Domestico = `Ovis aries` + Bos + Capra + Ovis)

taxa_cols <- c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa",
               "Caprinae", "Rupicapra rupicapra", "Domestico")

for (col in taxa_cols) {
  df[[col]] <- as.numeric(gsub(",", ".", as.character(df[[col]])))
}

metadati <- df %>% select(Sample_ID, Popolazione)
matrice_prede <- df %>% select(all_of(taxa_cols))

cat("N campioni:", nrow(matrice_prede), " N categorie:", ncol(matrice_prede), "\n")

# ==============================================================================
# 2. TRASFORMAZIONE DI HELLINGER E DISTANZA DI BRAY-CURTIS
# ==============================================================================
matrice_hellinger <- decostand(matrice_prede, method = "hellinger")
dist_matrix <- vegdist(matrice_hellinger, method = "bray")

# ==============================================================================
# 3. PCoA (Principal Coordinate Analysis)
# ==============================================================================
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Percentuale di varianza spiegata per i primi due assi (solo autovalori positivi)
eig_positivi <- pcoa_result$eig[pcoa_result$eig > 0]
var_asse1 <- round(100 * pcoa_result$eig[1] / sum(eig_positivi), 2)
var_asse2 <- round(100 * pcoa_result$eig[2] / sum(eig_positivi), 2)
cat("Varianza spiegata: asse 1 =", var_asse1, "% , asse 2 =", var_asse2, "%\n")

site_scores <- as.data.frame(pcoa_result$points)
colnames(site_scores) <- c("PCoA1", "PCoA2")
site_scores$Popolazione <- metadati$Popolazione

# ==============================================================================
# 4. SCORE DELLE SPECIE (medie pesate — wascores, standard per PCoA/NMDS)
# ==============================================================================
species_scores <- as.data.frame(wascores(site_scores[, c("PCoA1", "PCoA2")], matrice_prede))
species_scores$Taxon <- rownames(species_scores)

# ==============================================================================
# 5. GRAFICO
# ==============================================================================
palette_nazioni <- c("Croazia" = "#E69F00", "Slovenia" = "#56B4E9")

plot_pcoa <- ggplot(site_scores, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Popolazione, fill = Popolazione),
             size = 2.5, alpha = 0.6,
             position = position_jitter(width = 0.02, height = 0.02)) +
  stat_ellipse(aes(color = Popolazione, fill = Popolazione),
               geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.8) +
  geom_point(data = species_scores, aes(x = PCoA1, y = PCoA2),
             color = "#B22222", size = 3) +
  geom_text_repel(data = species_scores, aes(x = PCoA1, y = PCoA2, label = Taxon),
                   color = "black", fontface = "bold.italic", size = 4.2,
                   seed = 123, max.overlaps = 20) +
  scale_color_manual(values = palette_nazioni) +
  scale_fill_manual(values = palette_nazioni) +
  labs(
    title = "Dietary Niche Space (Principal Coordinate Analysis - PCoA)",
    subtitle = paste0("Bray-Curtis distance on Hellinger-transformed RRA% | Variance explained: ",
                       var_asse1, "% + ", var_asse2, "%"),
    x = paste0("Coordinate 1 (", var_asse1, "%)"),
    y = paste0("Coordinate 2 (", var_asse2, "%)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

print(plot_pcoa)
ggsave("Grafico_PCoA_Definitivo_Tesi.png", plot_pcoa, width = 10, height = 7.5, dpi = 300)

cat("\nGrafico salvato: Grafico_PCoA_Definitivo_Tesi.png\n")
# Fine del Master Script
# ==============================================================================
# ATTO 1 ESTESO — Selezione della preda principale, con Cinghiale integrato
# Analisi gerarchica a due livelli + grafico finale a due pannelli.
# Autosufficiente: nessuna dipendenza da script precedenti.
# Usa il file Community_Matrix_SloCro.csv aggiornato (colonna "Geographic area").
# ==============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)   # install.packages("patchwork") se non ce l'hai

# ==============================================================================
# 1. CARICAMENTO E PREPARAZIONE DATI
# ==============================================================================
community <- read_csv2("Community_Matrix_SloCro.csv", locale = locale(decimal_mark = ","))
ndvi      <- read_csv2("Wolf_NDVI_Dynamic_Dinaric.csv", locale = locale(decimal_mark = ","))

# Rinomino subito la colonna togliendo lo spazio - elimina alla radice ogni
# problema di backtick/virgolette nei grafici successivi
community <- community %>% rename(Geographic_Area = `Geographic area`)

df <- community %>%
  inner_join(ndvi %>% select(Sample, ELEV_mean, NDVI_mean, NDVI_sd),
             by = c("Sample_ID" = "Sample"))

stopifnot(nrow(df) == 100)

# Standardizzazione predittori (coerente con l'analisi Python già validata)
df$ELEV_mean_z <- as.numeric(scale(df$ELEV_mean))
df$NDVI_mean_z <- as.numeric(scale(df$NDVI_mean))
df$NDVI_sd_z   <- as.numeric(scale(df$NDVI_sd))

# ==============================================================================
# 2. UNIVERSO DI ANALISI — campioni con almeno una tra Cervo/Capriolo/Cinghiale
#    come specie dominante (max %RRA tra le tre)
# ==============================================================================
df_sub <- df %>%
  filter(`Cervus elaphus` > 0 | `Capreolus capreolus` > 0 | `Sus scrofa` > 0) %>%
  rowwise() %>%
  mutate(
    dominante = c("Cervus elaphus", "Capreolus capreolus", "Sus scrofa")[
      which.max(c_across(c(`Cervus elaphus`, `Capreolus capreolus`, `Sus scrofa`)))
    ]
  ) %>%
  ungroup()

cat("N campioni nell'universo a 3 prede:", nrow(df_sub), "\n")
cat("Distribuzione preda dominante:\n")
print(table(df_sub$dominante))

# ==============================================================================
# 3. LIVELLO 1 — Cinghiale-dominante vs Cervidi-dominante (Cervo+Capriolo insieme)
# ==============================================================================
df_sub$is_boar_dominant <- as.integer(df_sub$dominante == "Sus scrofa")

modello_livello1 <- glm(
  is_boar_dominant ~ ELEV_mean_z + NDVI_mean_z,
  data = df_sub,
  family = binomial(link = "logit")
)

cat("\n", strrep("=", 70), "\n")
cat("LIVELLO 1 — Cinghiale-dominante vs Cervidi-dominante\n")
cat(strrep("=", 70), "\n")
print(summary(modello_livello1))
cat("\nOdds Ratio e IC 95%:\n")
print(exp(cbind(OR = coef(modello_livello1), confint(modello_livello1))))

# ==============================================================================
# 4. LIVELLO 2 — dentro i Cervidi: Cervo vs Capriolo
# ==============================================================================
df_cervidi <- df_sub %>% filter(dominante != "Sus scrofa")
df_cervidi$is_cervo_dominant <- as.integer(df_cervidi$dominante == "Cervus elaphus")

modello_livello2 <- glm(
  is_cervo_dominant ~ ELEV_mean_z + NDVI_mean_z,
  data = df_cervidi,
  family = binomial(link = "logit")
)

cat("\n", strrep("=", 70), "\n")
cat("LIVELLO 2 — dentro i Cervidi: Cervo vs Capriolo\n")
cat(strrep("=", 70), "\n")
print(summary(modello_livello2))
cat("\nOdds Ratio e IC 95%:\n")
print(exp(cbind(OR = coef(modello_livello2), confint(modello_livello2))))

# ==============================================================================
# 5. CONTROLLO SUPPLEMENTARE — NDVI_sd come predittore del Cinghiale
#    (documentato, non nel modello finale: non significativo)
# ==============================================================================
modello_check_ndvi_sd <- glm(
  as.integer(`Sus scrofa` > 0) ~ NDVI_sd_z,
  data = df,
  family = binomial(link = "logit")
)
cat("\n", strrep("=", 70), "\n")
cat("CONTROLLO — presenza Cinghiale ~ NDVI_sd (per documentazione, non nel modello finale)\n")
cat(strrep("=", 70), "\n")
print(summary(modello_check_ndvi_sd))

# ==============================================================================
# 6. GRAFICO — due pannelli affiancati
# ==============================================================================
palette_area <- c("Croazia" = "#C73E1D", "Slovenia" = "#2E86AB")
# Stessa palette del grafico HHH. Se la RDA usa ancora "#E69F00"/"#56B4E9",
# allineala anche lei a questi due colori per coerenza tra le tre figure.

elev_range <- range(df_sub$ELEV_mean)
grid_elev <- seq(elev_range[1], elev_range[2], length.out = 200)
grid_elev_z <- (grid_elev - mean(df_sub$ELEV_mean)) / sd(df_sub$ELEV_mean)

predici_curva <- function(modello, grid_z, ndvi_z_fisso = 0) {
  nd <- data.frame(ELEV_mean_z = grid_z, NDVI_mean_z = ndvi_z_fisso)
  pred <- predict(modello, newdata = nd, type = "link", se.fit = TRUE)
  p <- plogis(pred$fit)
  lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  data.frame(ELEV_mean = grid_elev, p = p, lower = lower, upper = upper)
}

curva1 <- predici_curva(modello_livello1, grid_elev_z)
curva2 <- predici_curva(modello_livello2, grid_elev_z)

# --- Pannello A: Livello 1 ---
panelA <- ggplot() +
  geom_ribbon(data = curva1, aes(x = ELEV_mean, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.4) +
  geom_line(data = curva1, aes(x = ELEV_mean, y = p), color = "black", linewidth = 1) +
  geom_jitter(data = df_sub,
              aes(x = ELEV_mean, y = is_boar_dominant, color = Geographic_Area),
              width = 0, height = 0.03, size = 2.2, alpha = 0.7) +
  scale_color_manual(values = palette_area, name = "Geographic Area") +
  scale_y_continuous(breaks = c(0, 1),
                      labels = c("Cervid-dominated", "Boar-dominated")) +
  labs(title = "A) Wild boar vs. cervids",
       subtitle = "p = 0.011 (elevation)",
       x = "Elevation (m)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# --- Pannello B: Livello 2 ---
panelB <- ggplot() +
  geom_ribbon(data = curva2, aes(x = ELEV_mean, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.4) +
  geom_line(data = curva2, aes(x = ELEV_mean, y = p), color = "black", linewidth = 1) +
  geom_jitter(data = df_cervidi,
              aes(x = ELEV_mean, y = is_cervo_dominant, color = Geographic_Area),
              width = 0, height = 0.03, size = 2.2, alpha = 0.7) +
  scale_color_manual(values = palette_area, name = "Geographic Area") +
  scale_y_continuous(breaks = c(0, 1),
                      labels = c("Roe deer-dominated", "Red deer-dominated")) +
  labs(title = "B) Red deer vs. roe deer (within cervids)",
       subtitle = "p = 0.581 (elevation, not significant)",
       x = "Elevation (m)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# --- Figura finale ---
figura_finale <- panelA + panelB +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(figura_finale)
ggsave("Grafico_Gerarchico_Cinghiale_Cervidi.png", figura_finale,
       width = 12, height = 6, dpi = 300)

cat("\nGrafico salvato: Grafico_Gerarchico_Cinghiale_Cervidi.png\n")

# ==============================================================================
# HABITAT HETEROGENEITY HYPOTHESIS (HHH) — versione finale
# Inglese, colonna Geographic_Area (senza spazi), palette coerente con i
# grafici RDA e Cinghiale/Cervidi (Croazia = #C73E1D, Slovenia = #2E86AB).
# Autosufficiente: nessuna dipendenza da script precedenti.
# ==============================================================================

library(readr)
library(dplyr)
library(ggplot2)

# ==============================================================================
# 1. CARICAMENTO DATI
# ==============================================================================
community <- read_csv2("Community_Matrix_SloCro.csv", locale = locale(decimal_mark = ","))
ndvi      <- read_csv2("Wolf_NDVI_Dynamic_Dinaric.csv", locale = locale(decimal_mark = ","))

# Rinomino subito la colonna, stessa logica usata per gli altri script
community <- community %>% rename(Geographic_Area = `Geographic area`)

species_cols <- c("Capreolus capreolus", "Caprinae", "Cervus elaphus", "Ovis aries",
                   "Rupicapra rupicapra", "Sus scrofa", "Bos", "Capra", "Lepus", "Ovis")

community$Richness_per_sample <- rowSums(community[species_cols] > 0)

df <- community %>%
  inner_join(ndvi %>% select(Sample, ELEV_mean, NDVI_mean, NDVI_sd),
             by = c("Sample_ID" = "Sample"))

stopifnot(nrow(df) == 100)

# ==============================================================================
# 2. VARIABILE RISPOSTA BINARIA E STANDARDIZZAZIONE
# ==============================================================================
df$Dieta_mista <- as.integer(df$Richness_per_sample > 1)

df$NDVI_sd_z   <- as.numeric(scale(df$NDVI_sd))
df$ELEV_mean_z <- as.numeric(scale(df$ELEV_mean))
df$Geographic_Area <- relevel(factor(df$Geographic_Area), ref = "Slovenia")

# ==============================================================================
# 3. MODELLO LOGISTICO
# ==============================================================================
modello_hhh <- glm(Dieta_mista ~ NDVI_sd_z + ELEV_mean_z + Geographic_Area,
                    data = df, family = binomial(link = "logit"))

print(summary(modello_hhh))
cat("\nOdds Ratio e IC 95%:\n")
print(exp(cbind(OR = coef(modello_hhh), confint(modello_hhh))))

# ==============================================================================
# 4. GRAFICO
# ==============================================================================
palette_area <- c("Croazia" = "#C73E1D", "Slovenia" = "#2E86AB")
# Stessa palette usata nei grafici RDA e Cinghiale/Cervidi - coerenza tra
# tutte e tre le figure ambientali della tesi.

p <- ggplot(df, aes(x = NDVI_sd, y = Dieta_mista)) +
  geom_jitter(aes(color = Geographic_Area), width = 0, height = 0.03, size = 2.5, alpha = 0.7) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, color = "black", fill = "grey70") +
  scale_color_manual(values = palette_area, name = "Geographic Area") +
  scale_y_continuous(breaks = c(0, 1), labels = c("Pure\n(1 taxon)", "Mixed\n(>1 taxon)")) +
  labs(title = "Habitat Heterogeneity Hypothesis (HHH)",
       subtitle = "Probability of mixed diet as a function of environmental heterogeneity",
       x = "Habitat Heterogeneity (NDVI Standard Deviation)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p)
ggsave("HHH_Logistic_Plot_Tesi.png", plot = p, width = 10, height = 6.5, dpi = 200)

cat("\nGrafico salvato: HHH_Logistic_Plot_Tesi.png\n")

