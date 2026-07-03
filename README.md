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
