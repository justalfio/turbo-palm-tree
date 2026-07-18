 =========================================================================
# CREAZIONE GRAFICI PUBBLICABILI# Ripristinato dataset completo + Visualizzazione in RStudio
# =========================================================================
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
CARICAMENTO DATI
df_indici <- read_csv2("FOO&RRA Slovenia and Croatia.csv", show_col_types = FALSE)

# Trasformiamo i NA in 0
df_indici[is.na(df_indici)] <- 0

# STREAMING_CHUNK: CATEGORIZZAZIONE E COLORI (TUTTI I TAXA INCLUSI)
# Creiamo una colonna "Categoria" includendo anche le famiglie e sottofamiglie
df_indici <- df_indici %>%
  mutate(Categoria = case_when(
    Preda %in% c("Capreolus capreolus", "Cervus elaphus", "Cervinae", "Sus scrofa", "Rupicapra rupicapra") ~ "wild ungulates",
    Preda %in% c("Ovis aries", "Caprinae", "Ovis", "Capra") ~ "sheep & goats",
    Preda %in% c("Bos") ~ "cattle",
    TRUE ~ "small fauna & other" # Comprende Lepus
  ))

# Definiamo la palette di colori ESATTA ispirata allo screenshot del paper
colori_jura <- c(
  "wild ungulates" = "#007A4D",       # Verde scuro
  "sheep & goats" = "#E2D19A",        # Beige
  "cattle" = "#B3987E",               # Marroncino chiaro
  "small fauna & other" = "#95C78F"   # Verde chiaro
)

# Fissiamo l'ordine della legenda come nel paper
df_indici$Categoria <- factor(df_indici$Categoria,
                              levels = c("wild ungulates", "sheep & goats", "small fauna & other", "cattle"))

# STREAMING_CHUNK: FUNZIONE GENERATRICE DEL GRAFICO
# Creiamo una funzione che prende i dati e il nome della nazione e "sputa" il grafico perfetto
crea_grafico_nazione <- function(dati, nazione) {
  
  # Filtriamo solo i dati della nazione richiesta
  df_nazione <- dati %>% filter(Popolazione == nazione)
  
  # Rimuoviamo le prede che in questa nazione sono a 0 assoluto (per non avere barre vuote)
  df_nazione <- df_nazione %>% filter(FOO_percentuale > 0 | RRA_percentuale > 0)
  
  # Ordiniamo le prede dalla più mangiata alla meno mangiata basandoci sulla RRA
  ordine_prede <- df_nazione %>%
    arrange(RRA_percentuale) %>%
    pull(Preda)
  
  df_nazione$Preda <- factor(df_nazione$Preda, levels = ordine_prede)
  
  # "Allunghiamo" il dataset per creare i due pannelli (Facet)
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
  
  # Calcoliamo il limite massimo dell'asse X per dare spazio al numerino finale
  max_x <- max(df_long$Valore) * 1.15
  
  # Creazione del Plot
  p <- ggplot(df_long, aes(x = Valore, y = Preda, fill = Categoria)) +
    # Le barre
    geom_bar(stat = "identity", width = 0.8) +
    
    # I NUMERINI alla fine della barra (con 1 decimale, come nel paper)
    geom_text(aes(label = round(Valore, 1)), hjust = -0.2, size = 3.5, color = "black") +
    
    # Creiamo i due pannelli separati
    facet_wrap(~ Metrica, scales = "free_x") +
    
    # Applichiamo i colori scelti
    scale_fill_manual(values = colori_jura) +
    
    # Espandiamo l'asse X in base al massimo
    scale_x_continuous(limits = c(0, max_x), expand = c(0, 0)) +
    
    # Titoli
    labs(
      title = paste("Diet Composition -", nazione),
      x = "Percentage (%)",
      y = NULL,
      fill = NULL
    ) +
    
    # Stile e formattazione (esattamente come lo screenshot)
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "italic", size = 12, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      
      # Titolo del grafico (Slovenia o Croazia)
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
      
      # Box dei titoli dei pannelli
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      
      # Griglia (solo verticale, come nel paper)
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Legenda in basso
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.5, "cm")
    )
  
  return(p)
}

# STREAMING_CHUNK: GENERAZIONE E VISUALIZZAZIONE
plot_slo <- crea_grafico_nazione(df_indici, "Slovenia")
plot_cro <- crea_grafico_nazione(df_indici, "Croazia")

# MOSTRA I GRAFICI DIRETTAMENTE IN RSTUDIO (Nel pannello "Plots" in basso a destra)
# Quando esegui queste due righe, i grafici appariranno uno dopo l'altro nel visualizzatore.
# Puoi usare le freccette "<- ->" nel pannello Plots per scorrere tra Slovenia e Croazia.
print(plot_slo)
print(plot_cro)

# SALVATAGGIO OPZIONALE (Togli il '#' se vuoi salvare i file PNG)
# ggsave("Grafico_Dieta_Slovenia_Jura.png", plot = plot_slo, width = 10, height = 5, dpi = 300)
# ggsave("Grafico_Dieta_Croazia_Jura.png", plot = plot_cro, width = 10, height = 5, dpi = 300)

print("Codice eseguito: controlla il pannello 'Plots' in basso a destra su RStudio!")


#COMMUNITY MATRIX per statistica multivariata
library(dplyr)
> library(tidyr)
> library(readr)
> 
> # 1. CARICAMENTO DEI DATASET PURIFICATI
> # Usiamo i file esportati alla fine della pipeline bioinformatica
> file_slo <- "Slovenia wolves FILE DEFINITIVO.csv"
> file_cro <- "Croatian wolves FILE DEFINITIVO.csv"
> 
> # Se hai salvato il file sloveno con un nome diverso (es. "Slovenia_Dataset_Purificato.csv"),
> # correggi il nome qui sopra!
> df_slo <- read_csv2(file_slo, show_col_types = FALSE)
ℹ Using "','" as decimal and "'.'" as grouping mark. Use `read_delim()` for more control.

[1mindexing[0m [34mSlovenia wolves FILE DEFINITIVO.csv[0m [=============] [32m53.91MB/s[0m, eta: [36m 0s[0m
                                                                                                                   
> df_cro <- read_csv2(file_cro, show_col_types = FALSE)
ℹ Using "','" as decimal and "'.'" as grouping mark. Use `read_delim()` for more control.

[1mindexing[0m [34mCroatian wolves FILE DEFINITIVO.csv[0m [============] [32m111.68MB/s[0m, eta: [36m 0s[0m
                                                                                                                   
> 
> # 2. FUNZIONE CORAZZATA PER TRASPORRE E PREPARARE I DATI
> prepara_matrice <- function(df, nome_popolazione) {
+   
+   # Identifichiamo dinamicamente la colonna delle specie (risolve il problema tra scientific_name e final scientific_name)
+   col_specie <- names(df)[grepl("scientific_name", names(df), ignore.case = TRUE)][1]
+   
+   # Rinominiamo la colonna specie in "Taxon" per uniformità
+   names(df)[names(df) == col_specie] <- "Taxon"
+   
+   # Rimuoviamo in blocco tutte le possibili colonne di metadati inutili (che si chiamino rank o final rank)
+   df <- df %>% select(-any_of(c("rank", "final rank", "species_list")))
+   
+   # Trasformiamo la tabella: da "larga" a "lunga" e poi di nuovo "larga" al contrario
+   df_long <- df %>% 
+     pivot_longer(cols = -Taxon, names_to = "Sample_ID", values_to = "Reads")
+   
+   df_wide <- df_long %>% 
+     pivot_wider(names_from = Taxon, values_from = Reads, values_fill = 0)
+   
+   # Aggiungiamo la colonna fondamentale per il test statistico!
+   df_wide$Popolazione <- nome_popolazione
+   
+   return(df_wide)
+ }
> 
> # Applichiamo la funzione ai due dataset
> mat_slo <- prepara_matrice(df_slo, "Slovenia")
> mat_cro <- prepara_matrice(df_cro, "Croazia")
> mat_slo
# A tibble: 49 × 9
   Sample_ID `Capreolus capreolus` Caprinae Cervinae `Ovis aries`
   <chr>                     <dbl>    <dbl>    <dbl>        <dbl>
 1 M1XTJ                      8430        0     3955            0
 2 M2CAA                     14333        0        0            0
 3 M2CCK                      5705        0        0            0
 4 M2E0E                      3438        0        0            0
 5 M2E30                     15513        0        0            0
 6 MSV013                        0        0    17968            0
 7 MSV045                     2415        0        0            0
 8 MSV04C                     3975        0        0            0
 9 MSV04E                     9215        0        0            0
10 MSV066                        0    31976        0            0
# ℹ 39 more rows
# ℹ 4 more variables: `Rupicapra rupicapra` <dbl>, `Sus scrofa` <dbl>,
#   `NA` <dbl>, Popolazione <chr>
# ℹ Use `print(n = ...)` to see more rows
> community_matrix <- bind_rows(mat_slo, mat_cro)
> community_matrix[is.na(community_matrix)] <- 0
> colonne_prede <- setdiff(names(community_matrix), c("Sample_ID", "Popolazione"))
> community_matrix[colonne_prede] <- t(apply(community_matrix[colonne_prede], 1, function(x) {
+   totale <- sum(x)
+   if(totale > 0) { 
+     return((x / totale) * 100) 
+   } else { 
+     return(x) 
+   }
+ }))
> community_matrix <- community_matrix %>%
+   select(Sample_ID, Popolazione, everything())
> write.csv2(community_matrix, "Community_Matrix_SloCro.csv", row.names = FALSE)
> 
> print("Community Matrix creata con successo!")
[1] "Community Matrix creata con successo!"
> print(paste("Totale campioni (righe):", nrow(community_matrix)))
[1] "Totale campioni (righe): 100"
> print(paste("Totale specie preda (colonne):", length(colonne_prede)))
[1] "Totale specie preda (colonne): 12"


# =========================================================================
# FASE 3: STATISTICA MULTIVARIATA (NMDS, PERMANOVA, SIMPER)
# Obiettivo: Testare la significatività statistica della differenza di dieta 
# lungo il gradiente ambientale Slovenia-Croazia.
# =========================================================================

# 0. Setup iniziale dei pacchetti (decommentare se serve installarli)
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("readr")

library(vegan)
library(ggplot2)
library(dplyr)
library(readr)

# 1. CARICAMENTO DELLA COMMUNITY MATRIX
# Importo la matrice unita (campioni x prede) che ho generato nello script precedente
df <- read_csv2("Community_Matrix_SloCro.csv", show_col_types = FALSE)

# Piccolo fix di sicurezza: levo un'eventuale colonna fantasma chiamata "NA" 
# che a volte si crea durante l'unione dei dataset se ci sono mismatch nei nomi
df <- df %>% select(-any_of("NA"))

# 2. PREPARAZIONE DEI DATI PER VEGAN
# Il pacchetto vegan è rigido: la matrice deve contenere SOLO dati numerici (le abbondanze). 
# Quindi stacco la colonna con la nazione (metadati) per usarla poi come fattore di raggruppamento.
metadati <- df %>% select(Sample_ID, Popolazione)
matrice_prede <- df %>% select(-Sample_ID, -Popolazione)

# Forzo la matrice a numerico per evitare errori fastidiosi durante il calcolo delle distanze
matrice_prede <- as.data.frame(lapply(matrice_prede, as.numeric))

# 3. NMDS (Non-metric Multidimensional Scaling)
# Calcolo la distanza di Bray-Curtis e proietto i miei campioni in uno spazio 2D (k=2)
# Imposto un seed fisso (123) così l'ordinamento spaziale mi viene sempre uguale ogni volta che runno lo script
set.seed(123) 
print("Calcolo dell'NMDS in corso... (potrebbe richiedere qualche secondo)")
nmds_risultato <- metaMDS(matrice_prede, distance = "bray", k = 2, trymax = 100, trace = FALSE)

# Estraggo le coordinate X e Y dei singoli campioni per poter fare il plot personalizzato con ggplot2
coordinate_nmds <- as.data.frame(scores(nmds_risultato, display = "sites"))
coordinate_nmds$Popolazione <- metadati$Popolazione

# 4. VISUALIZZAZIONE NMDS CON GGPLOT2 (Stile Paper)
# Palette personalizzata per mantenere la coerenza con i colori dello studio
palette_nazioni <- c("Slovenia" = "#56B4E9", "Croazia" = "#E69F00")

plot_nmds <- ggplot(coordinate_nmds, aes(x = NMDS1, y = NMDS2, color = Popolazione, fill = Popolazione)) +
  geom_point(size = 2.5, alpha = 0.8) +
  # Aggiungo le ellissi di confidenza al 95% per far vedere bene i due cluster
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.8) + 
  scale_color_manual(values = palette_nazioni) +
  scale_fill_manual(values = palette_nazioni) +
  labs(
    title = "NMDS - Diet Composition (Slovenia vs Croatia)",
    # Stampo il valore di stress nel sottotitolo per confermare l'affidabilità del 2D (se < 0.20 è ottimo)
    subtitle = paste("2D Stress value:", round(nmds_risultato$stress, 3)), 
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

# Visualizzo il grafico e lo esporto in alta risoluzione per la tesi
print(plot_nmds)
ggsave("Grafico_NMDS_Dieta.png", plot = plot_nmds, width = 8, height = 6, dpi = 300)

# 5. PERMANOVA (Il test statistico definitivo)
# Questo è il cuore dell'analisi: mi serve il p-value per confermare alla commissione 
# che la spaccatura ecologica tra le due diete è statisticamente significativa.
set.seed(123)
permanova_risultato <- adonis2(matrice_prede ~ Popolazione, data = metadati, method = "bray", permutations = 9999)

print(" ")
print("==================================================")
print("             RISULTATO PERMANOVA                  ")
print("==================================================")
print(permanova_risultato)
# NB: Guardo la colonna "Pr(>F)" per il p-value e "R2" per la varianza spiegata dal gradiente.

# 6. SIMPER (Similarity Percentages)
# Ok, le diete sono diverse. Ma quali prede creano la spaccatura? 
# Uso SIMPER per ottenere la classifica delle specie chiave (spoiler: saranno Cinghiale e Capriolo/Cervo).
simper_risultato <- simper(matrice_prede, group = metadati$Popolazione)

print(" ")
print("==================================================")
print("               RISULTATO SIMPER                   ")
print("==================================================")
print(summary(simper_risultato))

# =========================================================================
# FASE 3: STATISTICA MULTIVARIATA (NMDS, PERMANOVA, SIMPER)
# Obiettivo: Testare la significatività statistica della differenza di dieta 
# lungo il gradiente ambientale (Slovenia vs Croazia) e visualizzare i vettori
# delle prede chiave in stile paper (es. Kunz et al. / Jura Mountains).
# =========================================================================

# 0. SETUP INIZIALE E LIBRERIE
# Verifico se i pacchetti sono installati, altrimenti li installo e li carico
pacchetti_necessari <- c("vegan", "ggplot2", "dplyr", "readr")
for (pkg in pacchetti_necessari) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 1. CARICAMENTO DELLA COMMUNITY MATRIX
# Importo la matrice (campioni sulle righe, prede sulle colonne)
df <- read_csv2("Community_Matrix_SloCro.csv", show_col_types = FALSE)

# Pulizia di sicurezza: rimuovo eventuali colonne fantasma o artefatti di join
df <- df %>% select(-any_of("NA"))

# NOTA ECOLOGICA SULLA TASSONOMIA:
# In questa fase manteniamo "Cervinae" (per la Slovenia) distinto da "Cervus elaphus" (per la Croazia).
# Anche se ecologicamente simili, uniamo i taxa solo se la bioinformatica ci dà certezza di specie.
# Essendo "Cervinae" una sottofamiglia (che potenzialmente maschera reads irrisolte di capriolo),
# adottiamo un approccio conservativo e manteniamo il dato puro del database.

# 2. PREPARAZIONE DEI DATI PER VEGAN
# vegan richiede una matrice di soli numeri. Stacco i metadati (Nazione).
metadati <- df %>% select(Sample_ID, Popolazione)
matrice_prede <- df %>% select(-Sample_ID, -Popolazione)

# Forzo a numerico per evitare conflitti nel calcolo delle distanze
matrice_prede <- as.data.frame(lapply(matrice_prede, as.numeric))

# 3. NMDS (Non-metric Multidimensional Scaling)
# Calcolo la distanza di Bray-Curtis (standard per la dieta) e proietto in 2D (k=2).
# autotransform = FALSE per non alterare le nostre RRA% già calcolate.
set.seed(123) 
print("Calcolo dell'NMDS e dei vettori in corso...")
nmds_risultato <- metaMDS(matrice_prede, distance = "bray", k = 2, trymax = 100, trace = FALSE, autotransform = FALSE)

# Estraggo le coordinate spaziali dei campioni
coordinate_nmds <- as.data.frame(scores(nmds_risultato, display = "sites"))
coordinate_nmds$Popolazione <- metadati$Popolazione

# 4. CALCOLO DEI VETTORI DELLE PREDE (ENVFIT - STILE PAPER)
# Calcolo in quale direzione specifica le prede "tirano" le diete dei lupi
set.seed(123)
fit_prede <- envfit(nmds_risultato, matrice_prede, permutations = 999)

# Estraggo le coordinate delle frecce (tengo solo le prede significative p < 0.05)
frecce_prede <- as.data.frame(scores(fit_prede, display = "vectors"))
frecce_prede$Preda <- rownames(frecce_prede)
frecce_prede$p_value <- fit_prede$vectors$pvals
frecce_prede <- frecce_prede %>% filter(p_value < 0.05) 

# Moltiplico le coordinate per allungare visivamente le frecce nel plot
moltiplicatore <- 1.5
frecce_prede$NMDS1 <- frecce_prede$NMDS1 * moltiplicatore
frecce_prede$NMDS2 <- frecce_prede$NMDS2 * moltiplicatore

# 5. VISUALIZZAZIONE NMDS CON GGPLOT2
palette_nazioni <- c("Slovenia" = "#56B4E9", "Croazia" = "#E69F00")

plot_nmds <- ggplot(coordinate_nmds, aes(x = NMDS1, y = NMDS2)) +
  # Uso geom_jitter per sparpagliare i campioni con dieta monospecifica (es. 100% capriolo) 
  # evitando che si sovrappongano in un solo punto invisibile.
  geom_jitter(aes(color = Popolazione, fill = Popolazione), size = 2.5, alpha = 0.5, width = 0.05, height = 0.05) +
  
  # Ellissi di confidenza (95%)
  stat_ellipse(aes(color = Popolazione, fill = Popolazione), geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.8) + 
  
  # Frecce direzionali delle prede (Vettori)
  geom_segment(data = frecce_prede, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "black", linewidth = 0.8) +
  # Etichette dei vettori
  geom_text(data = frecce_prede, aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = Preda), 
            color = "black", fontface = "italic", size = 4.5) +

  scale_color_manual(values = palette_nazioni) +
  scale_fill_manual(values = palette_nazioni) +
  labs(
    title = "NMDS - Diet Composition Driven by Key Prey",
    subtitle = paste("2D Stress:", round(nmds_risultato$stress, 3), "| Arrows show significant prey vectors"), 
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

# Mostro il grafico e lo esporto ad alta risoluzione
print(plot_nmds)
ggsave("Grafico_NMDS_Frecce_Paper.png", plot = plot_nmds, width = 9, height = 7, dpi = 300)

# 6. PERMANOVA (Test Statistico Definitivo)
# Verifico se la spaccatura spaziale osservata nell'NMDS è statisticamente significativa
set.seed(123)
permanova_risultato <- adonis2(matrice_prede ~ Popolazione, data = metadati, method = "bray", permutations = 9999)

print(" ")
print("==================================================")
print("             RISULTATO PERMANOVA                  ")
print("==================================================")
print(permanova_risultato)

# 7. SIMPER (Similarity Percentages)
# Ottengo la classifica esatta delle prede che contribuiscono maggiormente 
# alla diversità tra le due popolazioni
simper_risultato <- simper(matrice_prede, group = metadati$Popolazione)

print(" ")
print("==================================================")
print("               RISULTATO SIMPER                   ")
print("==================================================")
print(summary(simper_risultato))
