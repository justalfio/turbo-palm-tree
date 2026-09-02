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
# PARTE 3: STATISTICA MULTIVARIATA (NMDS, VETTORI ENVFIT, PERMANOVA, SIMPER)
# =========================================================================

# 1. CARICAMENTO DELLA COMMUNITY MATRIX APENA CREATA
df <- read_csv2("Community_Matrix_SloCro.csv", show_col_types = FALSE)

# Pulizia di sicurezza: rimuoviamo eventuali colonne fantasma o artefatti di join
df <- df %>% select(-any_of("NA"))

# NOTA ECOLOGICA SULLA TASSONOMIA:
# In questa fase manteniamo "Cervinae" (per la Slovenia) distinto da "Cervus elaphus" (per la Croazia).
# Essendo "Cervinae" una sottofamiglia (che potenzialmente maschera reads irrisolte di capriolo),
# adottiamo un approccio conservativo e manteniamo il dato puro del database senza forzare unioni.

# 2. PREPARAZIONE DEI DATI PER IL PACCHETTO VEGAN
# Il pacchetto vegan richiede una matrice di soli numeri senza metadati di colonna
metadati <- df %>% select(Sample_ID, Popolazione)
matrice_prede <- df %>% select(-Sample_ID, -Popolazione)

# Forziamo a numerico per evitare conflitti nel calcolo delle distanze
matrice_prede <- as.data.frame(lapply(matrice_prede, as.numeric))



# 6. PERMANOVA (Test Statistico Definitivo sulla dissimilarità trofica)
# Verifico se la diversità spaziale osservata nell'NMDS è statisticamente significativa lungo il gradiente
set.seed(123)
permanova_risultato <- adonis2(matrice_prede ~ Popolazione, data = metadati, method = "bray", permutations = 9999)

print(" ")
print("==================================================")
print("              RISULTATO PERMANOVA                 ")
print("==================================================")
print(permanova_risultato)
# NB: Il valore chiave si trova nella colonna Pr(>F) per il p-value e in R2 per la percentuale di varianza spiegata

# 7. SIMPER (Similarity Percentages)
# Ottengo la classifica quantitativa esatta delle prede responsabili della diversità tra Slovenia e Croazia
simper_risultato <- simper(matrice_prede, group = metadati$Popolazione)

print(" ")
print("==================================================")
print("                RISULTATO SIMPER                  ")
print("==================================================")
print(summary(simper_risultato))
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

# Fine del Master Script
# ==============================================================================
# ATTO 1 ESTESO — Selezione della preda principale, con Cinghiale integrato
# Analisi gerarchica a due livelli, come discusso con Marta.
# Autosufficiente: nessuna dipendenza da script precedenti.
# ==============================================================================

library(readr)
library(dplyr)

# ==============================================================================
# 1. CARICAMENTO E PREPARAZIONE DATI
# ==============================================================================
community <- read_csv2("Community_Matrix_SloCro.csv", locale = locale(decimal_mark = ","))
ndvi      <- read_csv2("Wolf_NDVI_Dynamic_Dinaric.csv", locale = locale(decimal_mark = ","))

df <- community %>%
  inner_join(ndvi %>% select(Sample, ELEV_mean, NDVI_mean, NDVI_sd),
             by = c("Sample_ID" = "Sample"))

stopifnot(nrow(df) == 100)

# Standardizzazione predittori (coerente con l'analisi Python gia' validata)
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
# 4. LIVELLO 2 — dentro i Cervidi: Cervo vs Capriolo (gia' noto, ricalcolato
#    qui sul sotto-campione corretto per coerenza con la struttura gerarchica)
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
# 5. CONTROLLO SUPPLEMENTARE (documentato, non nel modello principale) —
#    NDVI_sd come predittore del Cinghiale: testato su richiesta, risultato
#    non significativo (p~0.93-0.96), NDVI_sd medio identico con/senza
#    cinghiale in dieta. Non incluso nel modello finale: l'eterogeneita'
#    ambientale non discrimina il cinghiale, a differenza dell'elevazione.
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
