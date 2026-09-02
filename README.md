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

# 3. CALCOLO NMDS (Non-metric Multidimensional Scaling)
# Calcolo la distanza di Bray-Curtis (standard per la dieta) e proietto in 2D (k=2)
# Impostiamo autotransform = FALSE per non alterare le nostre RRA% già standardizzate
set.seed(123) 
print("Calcolo dell'NMDS e dei vettori envfit in corso...")
nmds_risultato <- metaMDS(matrice_prede, distance = "bray", k = 2, trymax = 100, trace = FALSE, autotransform = FALSE)

# Estraggo le coordinate spaziali dei campioni nello spazio di ordinamento
coordinate_nmds <- as.data.frame(scores(nmds_risultato, display = "sites"))
coordinate_nmds$Popolazione <- metadati$Popolazione

# 4. CALCOLO DEI VETTORI DELLE PREDE (ENVFIT - STILE KUNZ ET AL. / JURA)
# Calcoliamo in quale direzione specifica le prede tirano le diete dei lupi
set.seed(123)
fit_prede <- envfit(nmds_risultato, matrice_prede, permutations = 999)

# Estraggo le coordinate delle frecce tenendo solo le prede statisticamente significative (p < 0.05)
frecce_prede <- as.data.frame(scores(fit_prede, display = "vectors"))
frecce_prede$Preda <- rownames(frecce_prede)
frecce_prede$p_value <- fit_prede$vectors$pvals
frecce_prede <- frecce_prede %>% filter(p_value < 0.05) 

# Moltiplichiamo le coordinate per allungare visivamente le frecce nel grafico per maggiore leggibilità
moltiplicatore <- 1.5
frecce_prede$NMDS1 <- frecce_prede$NMDS1 * moltiplicatore
frecce_prede$NMDS2 <- frecce_prede$NMDS2 * moltiplicatore

# 5. VISUALIZZAZIONE NMDS CON VETTORI TRAMITE GGPLOT2
palette_nazioni <- c("Slovenia" = "#56B4E9", "Croazia" = "#E69F00")

plot_nmds <- ggplot(coordinate_nmds, aes(x = NMDS1, y = NMDS2)) +
  # Uso geom_jitter per sparpagliare leggermente i campioni con dieta monospecifica (es. 100% capriolo) 
  # evitando che si sovrappongano perfettamente in un solo punto invisibile
  geom_jitter(aes(color = Popolazione, fill = Popolazione), size = 2.5, alpha = 0.5, width = 0.05, height = 0.05) +
  
  # Ellissi di confidenza al 95% per visualizzare l'area di nicchia delle due popolazioni
  stat_ellipse(aes(color = Popolazione, fill = Popolazione), geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.8) + 
  
  # Frecce direzionali delle prede (Vettori significativi)
  geom_segment(data = frecce_prede, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "black", linewidth = 0.8) +
  # Etichette delle frecce tassonomiche
  geom_text(data = frecce_prede, aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = Preda), 
            color = "black", fontface = "italic", size = 4.5) +

  scale_color_manual(values = palette_nazioni) +
  scale_fill_manual(values = palette_nazioni) +
  labs(
    title = "NMDS - Diet Composition Driven by Key Prey",
    subtitle = paste("2D Stress:", round(nmds_risultato$stress, 3), "| Arrows show significant prey vectors (p < 0.05)"), 
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 11, hjust = 0.5, margin = margin(b=10)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Mostro il grafico NMDS a schermo e lo esporto ad alta risoluzione
print(plot_nmds)
ggsave("Grafico_NMDS_Frecce_Paper.png", plot = plot_nmds, width = 9, height = 7, dpi = 300)

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
