# Guida Didattica Completa: Meta-Analisi di Test Diagnostici per Biomarcatori di mTBI

## Indice
1. [Contesto Clinico e Obiettivo dello Studio](#contesto-clinico)
2. [Fondamenti Statistici dei Test Diagnostici](#fondamenti-statistici)
3. [Il Modello Bivariato di Reitsma](#modello-reitsma)
4. [Spiegazione del Codice R — Riga per Riga](#codice-r)
5. [Spiegazione del Codice Python — Riga per Riga](#codice-python)
6. [Interpretazione dei Risultati](#interpretazione-risultati)
7. [Glossario dei Termini Statistici](#glossario)

---

## 1. Contesto Clinico e Obiettivo dello Studio {#contesto-clinico}

### Che cos'è il mTBI?
Il **mild Traumatic Brain Injury (mTBI)**, o trauma cranico lieve, è una delle condizioni più comuni nei Pronto Soccorsi. Il problema clinico principale è: **quali pazienti devono ricevere una TAC cranica (CT)?**

- La TAC espone a radiazioni ionizzanti ed è costosa
- Non fare la TAC a un paziente con una lesione intracranica può essere fatale
- I **biomarcatori ematici** potrebbero aiutare a decidere chi necessita di TAC

### Cosa vogliamo ottenere?
Questa meta-analisi confronta **due biomarcatori ematici** per la loro capacità di prevedere lesioni intracraniche nei pazienti con mTBI:

| Biomarcatore | N. Studi | N. Pazienti |
|---|---|---|
| **GFAP + UCH-L1** | 8 | 4,969 |
| **S100B** | 6 | 3,582 |

**Domanda di ricerca:** Quale biomarcatore ha migliori proprietà diagnostiche? Quale ci permette di evitare più TAC inutili senza mancare lesioni pericolose?

---

## 2. Fondamenti Statistici dei Test Diagnostici {#fondamenti-statistici}

### La Tabella 2×2 (Tabella di Contingenza)

Ogni studio primario produce dati che possono essere organizzati così:

| | **CT Positiva (Malattia +)** | **CT Negativa (Malattia −)** |
|---|---|---|
| **Test Positivo** | TP (Veri Positivi) | FP (Falsi Positivi) |
| **Test Negativo** | FN (Falsi Negativi) | TN (Veri Negativi) |

### Le Misure Fondamentali

#### Sensibilità (Sensitivity, Se)
```
Se = TP / (TP + FN)
```
- **Domanda:** Tra tutti i pazienti CHE HANNO la lesione, quanti il test identifica correttamente?
- **Interpretazione clinica:** Una sensibilità del 95% significa che il test "vede" il 95% delle lesioni reali
- **Falsi Negativi:** (1 − Se) = proporzione di lesioni che il test NON vede

#### Specificità (Specificity, Sp)
```
Sp = TN / (TN + FP)
```
- **Domanda:** Tra tutti i pazienti che NON hanno la lesione, quanti il test identifica correttamente come sani?
- **Interpretazione clinica:** Una specificità del 40% significa che il 40% dei pazienti sani viene correttamente escluso
- **Falsi Positivi:** (1 − Sp) = proporzione di sani che il test sbaglia

#### Valore Predittivo Negativo (NPV)
```
NPV = [Sp × (1 − Prevalenza)] / [Sp × (1 − Prevalenza) + (1 − Se) × Prevalenza]
```
- **Domanda:** Se il test è negativo, qual è la probabilità che il paziente NON abbia la lesione?
- **Perché è cruciale:** Nel mTBI, vogliamo essere SICURI che un test negativo significhi "nessuna lesione". Un NPV > 95% è considerato eccellente per escludere la TAC.
- **Dipende dalla prevalenza:** A parità di Se e Sp, il NPV cambia con la frequenza della malattia nella popolazione

### Perché la Meta-Analisi?

I singoli studi hanno:
- **Campioni piccoli** → stime imprecise
- **Popolazioni diverse** → risultati variabili
- **Soglie diverse** del test → sensibilità e specificità differenti

La meta-analisi **combina** i risultati per ottenere una stima più precisa e generale.

---

## 3. Il Modello Bivariato di Reitsma {#modello-reitsma}

### Il Problema: Sensibilità e Specificità sono Correlate

Quando un test diagnostico ha una soglia variabile (es. "il biomarcatore è positivo se > X ng/L"):
- Se abbasso la soglia → più pazienti risultano positivi → **sensibilità aumenta**, ma **specificità diminuisce**
- Se alzo la soglia → **specificità aumenta**, ma **sensibilità diminuisce**

Questo fenomeno si chiama **trade-off sensibilità-specificità**. Le due misure NON sono indipendenti: sono **correlate negativamente**.

### La Soluzione: Modello Bivariato

Il modello di Reitsma (2005) tratta sensibilità e specificità **insieme**, modellando la loro correlazione. Ecco come funziona:

#### Passo 1: Trasformazione Logit

Le proporzioni (Se, Sp) sono vincolate tra 0 e 1. La trasformazione logit le mappa su tutta la retta reale:

```
logit(p) = ln(p / (1 − p))
```

**Perché?** La distribuzione normale (gaussiana) è definita su (−∞, +∞). Lavorare in scala logit ci permette di usare modelli lineari.

#### Passo 2: Modello ad Effetti Random

Ogni studio `i` ha la sua vera sensibilità e specificità, che variano attorno a un valore medio comune:

```
logit(Se_i) ~ N(μ_Se, σ²_Se + τ²_Se)
logit(Sp_i) ~ N(μ_Sp, σ²_Sp + τ²_Sp)
```

Dove:
- `μ_Se`, `μ_Sp` = valori medi pooled (ciò che stimiamo)
- `σ²` = varianza intra-studio (errore campionario)
- `τ²` = varianza inter-studio (eterogeneità vera)

#### Passo 3: Stima dei Parametri

Il modello stima:
1. **Sensibilità pooled** (media pesata delle sensibilità dei singoli studi)
2. **Specificità pooled** (media pesata delle specificità dei singoli studi)
3. **Correlazione** tra Se e Sp
4. **Intervalli di confidenza** al 95%

#### Passo 4: Back-Transformation

I risultati in scala logit vengono ritrasformati in proporzioni:

```
Se_pooled = expit(μ_Se) = 1 / (1 + exp(−μ_Se))
```

### Vantaggi del Modello Bivariato

✅ Tiene conto della correlazione Se-Sp
✅ Gestisce l'eterogeneità tra studi
✅ Produce stime più robuste rispetto ai metodi univariati
✅ Permette di costruire la curva SROC

### La Curva SROC (Summary Receiver Operating Characteristic)

La curva SROC mostra la relazione tra sensibilità e specificità attraverso tutti gli studi:
- **Asse X:** 1 − Specificità (tasso di falsi positivi)
- **Asse Y:** Sensibilità (tasso di veri positivi)
- **Punto riassuntivo:** La stima pooled di Se e Sp
- **Curva:** Mostra come Se e Sp variano al variare della soglia del test

Una curva che si avvicina all'angolo in alto a sinistra (Se=1, 1−Sp=0) indica un test eccellente.

---

## 4. Spiegazione del Codice R — Riga per Riga {#codice-r}

### SEZIONE 0: Installazione e Caricamento Pacchetti

```r
if (!requireNamespace("mada", quietly = TRUE)) install.packages("mada")
if (!requireNamespace("meta", quietly = TRUE)) install.packages("meta")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
```

**Cosa fa:** Controlla se ogni pacchetto è installato; se no, lo installa.
- `mada`: Meta-analysis of diagnostic accuracy — implementa il modello di Reitsma
- `meta`: Meta-analysis con metodi classici
- `ggplot2`: Grafica avanzata basata sulla "Grammar of Graphics"
- `gridExtra`: Combina più grafici in un'unica immagine

```r
library(mada)
library(meta)
library(ggplot2)
library(gridExtra)
```

**Cosa fa:** Carica i pacchetti nella sessione di lavoro.

---

### SEZIONE 1: Creazione del Dataset

```r
gfap_data <- data.frame(
  study = c(
    "Oris 2024 (CCLM)",
    "Bazarian 2021 (AEM)",
    "Papa 2022 (JAMA NO)",
    "Chayoua 2024 (JNT)",
    "Legramante 2024 (IJEM)",
    "Jones 2020 (Brain Inj)",
    "Lagares 2024 (EBioMed)",
    "Puravet 2025 (AEM)"
  ),
  TP = c(18, 24, 52, 38, 7, 33, 68, 99),
  FP = c(61, 199, 392, 130, 89, 425, 477, 285),
  FN = c(0, 1, 1, 1, 0, 0, 2, 1),
  TN = c(127, 167, 413, 149, 34, 221, 663, 373),
  stringsAsFactors = FALSE
)
```

**Cosa fa:** Crea un data frame con i dati degli 8 studi sul biomarcatore GFAP+UCH-L1.
- Ogni riga = uno studio
- `TP, FP, FN, TN` = le 4 celle della tabella 2×2
- `stringsAsFactors = FALSE`: impedisce la conversione automatica delle stringhe in fattori (best practice in R moderno)

```r
gfap_data$N    <- gfap_data$TP + gfap_data$FP + gfap_data$FN + gfap_data$TN
gfap_data$Se   <- gfap_data$TP / (gfap_data$TP + gfap_data$FN)
gfap_data$Sp   <- gfap_data$TN / (gfap_data$TN + gfap_data$FP)
gfap_data$CT_pos <- gfap_data$TP + gfap_data$FN
gfap_data$prev   <- gfap_data$CT_pos / gfap_data$N
```

**Cosa fa:** Calcola variabili derivate per ogni studio:
- `N`: numero totale di pazienti nello studio
- `Se`: sensibilità dello studio (TP / tutti i CT+)
- `Sp`: specificità dello studio (TN / tutti i CT−)
- `CT_pos`: numero totale di pazienti con CT positiva (TP + FN)
- `prev`: prevalenza di CT positive nello studio

```r
cat("=== GFAP+UCH-L1 POOL (8 studi) ===\n")
cat("N totale:", sum(gfap_data$N), "\n")
cat("CT+ totali:", sum(gfap_data$CT_pos), "\n\n")
```

**Cosa fa:** Stampa un riepilogo nella console. `cat()` concatena e stampa testo.

**Lo stesso blocco viene ripetuto per S100B** con `s100b_data`.

---

### SEZIONE 2: Modello di Reitsma

```r
fit_gfap <- reitsma(
  TP = gfap_data$TP, FP = gfap_data$FP,
  FN = gfap_data$FN, TN = gfap_data$TN,
  data = gfap_data
)
```

**Cosa fa:** Adatta il modello bivariato di Reitsma ai dati GFAP.
- La funzione `reitsma()` del pacchetto `mada` esegue la massima verosimiglianza (maximum likelihood) per stimare i parametri del modello
- Stima: μ_Se, μ_Sp, τ²_Se, τ²_Sp, e la correlazione tra Se e Sp

```r
print(summary(fit_gfap))
```

**Cosa fa:** Stampa il sommario del modello con tutte le stime, gli errori standard, e gli intervalli di confidenza.

```r
se_gfap <- plogis(fit_gfap$coefficients[1])
sp_gfap <- plogis(fit_gfap$coefficients[2])
```

**Cosa fa:** Estrae le stime pooled e le ritrasforma in scala originale.
- `fit_gfap$coefficients[1]` = μ_Se in scala logit
- `plogis()` = funzione expit (logistica inversa): `plogis(x) = 1 / (1 + exp(−x))`
- Risultato: sensibilità pooled in scala 0–1 (es. 0.923 = 92.3%)

---

### SEZIONE 3: NPV e Clinical Metrics

```r
prev_gfap_weighted <- sum(gfap_data$CT_pos) / sum(gfap_data$N)
```

**Cosa fa:** Calcola la prevalenza di CT+ **ponderata per N**.
- Somma tutti i CT+ di tutti gli studi, divide per il totale dei pazienti
- Questo dà più peso agli studi più grandi (più precisi)

```r
npv <- function(se, sp, prev) {
  (sp * (1 - prev)) / (sp * (1 - prev) + (1 - se) * prev)
}
```

**Cosa fa:** Definisce una funzione per calcolare il NPV usando il **teorema di Bayes**.
- Numeratore: probabilità di essere CT− E test negativo = Sp × (1 − prev)
- Denominatore: probabilità di test negativo (sia CT+ che CT−)
- Risultato: P(CT− | test negativo)

```r
npv_gfap <- npv(se_gfap, sp_gfap, prev_gfap_weighted)
```

**Cosa fa:** Calcola il NPV per GFAP usando le stime pooled.

```r
miss_gfap <- (1 - se_gfap) * prev_gfap_weighted * 1000
```

**Cosa fa:** Stima quante lesioni verrebbero **manca**te ogni 1000 pazienti testati.
- `(1 − Se)` = proporzione di lesioni non rilevate dal test
- `× prev` = proporzione di pazienti con lesioni nella popolazione
- `× 1000` = scala per 1000 pazienti (più intuitivo per i clinici)

```r
ct_avoid_gfap <- sp_gfap * (1 - prev_gfap_weighted) * 100
```

**Cosa fa:** Stima la percentuale di TAC che si potrebbero **evitare** con un test negativo.
- `Sp` = proporzione di CT− correttamente identificati
- `(1 − prev)` = proporzione di pazienti che sono davvero CT−
- Risultato: % di pazienti che non verrebbero irradiati inutilmente

---

### SEZIONE 4: Forest Plot

```r
png("Fig2_forest_gfap_v9.png", width = 2800, height = 1800, res = 200)
```

**Cosa fa:** Apre un dispositivo grafico PNG ad alta risoluzione (200 DPI).

```r
par(mfrow = c(1, 2), mar = c(5, 12, 4, 2), oma = c(0, 0, 3, 0))
```

**Cosa fa:** Configura il layout grafico.
- `mfrow = c(1, 2)`: 1 riga, 2 colonne (due grafici affiancati)
- `mar`: margini interni di ogni grafico (basso, sinistro, alto, destro)
- `oma`: margini esterni dell'intera figura

```r
forest(
  madad(TP = gfap_data$TP, FP = gfap_data$FP,
        FN = gfap_data$FN, TN = gfap_data$TN),
  type = "sens",
  main = "Sensibilità",
  xlab = "Sensibilità (95% CI)",
  snames = gfap_data$study,
  xlim = c(0.5, 1.05),
  col.diamond = "#2563EB",
  col.square  = "#93C5FD"
)
```

**Cosa fa:** Genera il forest plot della sensibilità.
- `madad()`: crea un oggetto di meta-analisi diagnostica dal pacchetto `mada`
- `type = "sens"`: mostra la sensibilità (usare `"spec"` per specificità)
- `snames`: etichette degli studi sull'asse Y
- `col.diamond`: colore del rombo pooled (blu `#2563EB`)
- `col.square`: colore dei quadrati dei singoli studi

**Cosa mostra il forest plot:**
- Ogni riga = uno studio con il suo punto stimato (quadrato) e intervallo di confidenza (linea orizzontale)
- La dimensione del quadrato è proporzionale al peso dello studio (studi più grandi = più precisi)
- Il rombo in basso = stima pooled con il suo IC 95%

```r
abline(v = se_gfap, col = "#2563EB", lty = 2, lwd = 2)
```

**Cosa fa:** Aggiunge una linea verticale tratteggiata in corrispondenza della stima pooled.

```r
mtext("Forest Plot — GFAP+UCH-L1 (8 studi, N=4,969)",
      outer = TRUE, cex = 1.4, font = 2, col = "#1e3a5f")
```

**Cosa fa:** Aggiunge il titolo esterno ai due grafici.
- `outer = TRUE`: posiziona il testo nei margini esterni
- `font = 2`: grassetto
- `cex = 1.4`: dimensione del carattere 1.4×

```r
dev.off()
```

**Cosa fa:** Chiude il dispositivo grafico e salva il file PNG.

---

### SEZIONE 6: Curva SROC Comparativa

```r
plot(fit_gfap,
     xlim = c(0, 1), ylim = c(0, 1),
     main = "SROC — GFAP+UCH-L1 vs S100B",
     xlab = "1 − Specificità (Tasso Falsi Positivi)",
     ylab = "Sensibilità (Tasso Veri Positivi)",
     col  = "#2563EB",
     lwd  = 2.5
)
```

**Cosa fa:** Plotta la curva SROC per GFAP.
- Il metodo `plot()` sull'oggetto `fit_gfap` genera automaticamente la curva SROC
- `xlab`: asse X = 1 − Sp (falsi positivi)
- `ylab`: asse Y = Se (veri positivi)

```r
lines(sroc(fit_s100b), col = "#16A34A", lwd = 2.5)
```

**Cosa fa:** Aggiunge la curva SROC di S100B sullo stesso grafico.
- `sroc()`: estrae i punti della curva SROC dal modello
- `lines()`: sovrappone le linee al grafico esistente

```r
points(1 - sp_gfap, se_gfap, pch = 23, bg = "#2563EB", cex = 2.5, col = "white")
points(1 - sp_s100b, se_s100b, pch = 22, bg = "#16A34A", cex = 2.5, col = "white")
```

**Cosa fa:** Aggiunge i punti sommario (rombi) per ciascun biomarcatore.
- `pch = 23`: rombo pieno (GFAP)
- `pch = 22`: quadrato pieno (S100B)
- `bg`: colore di riempimento
- `cex = 2.5`: dimensione 2.5×

```r
points(1 - gfap_data$Sp, gfap_data$Se, pch = 21,
       bg = "#BFDBFE", col = "#2563EB", cex = sqrt(gfap_data$N) / 8)
```

**Cosa fa:** Aggiunge i punti dei singoli studi.
- `cex = sqrt(N) / 8`: la dimensione dei punti è proporzionale alla radice di N (studi più grandi = punti più grandi)

```r
lines(fit_gfap, "conf", col = "#2563EB", lty = 2)
```

**Cosa fa:** Aggiunge l'ellisse di confidenza al 95% attorno al punto sommario.

```r
legend("bottomright",
       legend = c(
         sprintf("GFAP+UCH-L1: Se=%.1f%%, Sp=%.1f%%", se_gfap*100, sp_gfap*100),
         sprintf("S100B:        Se=%.1f%%, Sp=%.1f%%", se_s100b*100, sp_s100b*100)
       ),
       col = c("#2563EB", "#16A34A"),
       lwd = c(2.5, 2.5),
       pch = c(23, 22),
       pt.bg = c("#2563EB", "#16A34A"),
       pt.cex = 1.8,
       bty = "n",
       cex = 1.0
)
```

**Cosa fa:** Aggiunge la legenda nell'angolo in basso a destra.
- `sprintf()`: formatta le stringhe con i valori numerici
- `bty = "n"`: nessun bordo attorno alla legenda

---

### SEZIONE 7: Clinical Summary Figure (ggplot2)

```r
summary_df <- data.frame(
  Biomarker    = c("GFAP+UCH-L1", "S100B"),
  N_studi      = c(nrow(gfap_data), nrow(s100b_data)),
  N_pazienti   = c(sum(gfap_data$N), sum(s100b_data$N)),
  Sensibilita  = c(round(se_gfap  * 100, 1), round(se_s100b  * 100, 1)),
  Specificita  = c(round(sp_gfap  * 100, 1), round(sp_s100b  * 100, 1)),
  NPV          = c(round(npv_gfap  * 100, 2), round(npv_s100b  * 100, 2)),
  CT_evitate   = c(round(ct_avoid_gfap, 1),  round(ct_avoid_s100b, 1)),
  Miss_per1000 = c(round(miss_gfap, 1),       round(miss_s100b, 1))
)
```

**Cosa fa:** Crea un data frame di riepilogo con tutte le metriche cliniche.

```r
metrics_long <- data.frame(
  Biomarker = rep(c("GFAP+UCH-L1", "S100B"), each = 4),
  Metrica   = rep(c("Sensibilità (%)", "Specificità (%)",
                    "NPV (%)", "CT evitate (%)"), 2),
  Valore    = c(
    se_gfap*100,  sp_gfap*100,  npv_gfap*100,  ct_avoid_gfap,
    se_s100b*100, sp_s100b*100, npv_s100b*100, ct_avoid_s100b
  )
)
```

**Cosa fa:** Ristruttura i dati in formato "long" per ggplot2 (ogni riga = una osservazione).

```r
p <- ggplot(metrics_long, aes(x = Metrica, y = Valore, fill = Biomarker)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, alpha = 0.9) +
  geom_text(aes(label = round(Valore, 1)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("GFAP+UCH-L1" = "#2563EB", "S100B" = "#16A34A")) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    title    = "Confronto Clinico: GFAP+UCH-L1 vs S100B nel mTBI",
    subtitle = sprintf(
      "GFAP+UCH-L1: %d studi, N=%d | S100B: %d studi, N=%d",
      nrow(gfap_data), sum(gfap_data$N),
      nrow(s100b_data), sum(s100b_data$N)
    ),
    x     = NULL,
    y     = "Valore (%)",
    fill  = "Biomarker",
    caption = "Modello Reitsma (bivariate random-effects). CT evitate = pazienti con test negativo tra CT−."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    legend.position = "top",
    legend.title  = element_text(face = "bold"),
    axis.text.x   = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    plot.caption  = element_text(size = 9, color = "grey50", hjust = 0)
  )
```

**Cosa fa:** Costruisce un grafico a barre comparativo con ggplot2.
- `ggplot()`: inizializza il grafico con le mappature estetiche (`aes`)
- `geom_bar()`: barre con altezza pari ai valori (`stat = "identity"`)
- `position_dodge()`: barre affiancate (non sovrapposte)
- `geom_text()`: etichette con i valori sopra le barre
- `scale_fill_manual()`: colori personalizzati
- `theme_minimal()`: tema grafico pulito
- `theme()`: personalizzazioni avanzate (grassetto, allineamento, dimensioni)

---

### SEZIONE 8: Sensitivity Analysis

```r
s100b_no_pur <- s100b_data[-6, ]
```

**Cosa fa:** Crea un sottoinsieme dei dati S100B escludendo lo studio "Puravet 2026" (riga 6).

```r
fit_s100b_np <- reitsma(
  TP = s100b_no_pur$TP, FP = s100b_no_pur$FP,
  FN = s100b_no_pur$FN, TN = s100b_no_pur$TN,
  data = s100b_no_pur
)
```

**Cosa fa:** Riadatta il modello di Reitsma senza Puravet.

**Perché si fa:** La sensitivity analysis verifica se i risultati sono **robusti** o dipendono da un singolo studio. Se le stime cambiano poco, il risultato è affidabile.

---

## 5. Spiegazione del Codice Python — Riga per Riga {#codice-python}

### Import delle Librerie

```python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.special import expit, logit
from scipy.optimize import minimize
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
```

**Cosa fa:** Importa le librerie necessarie.
- `numpy`: calcolo numerico vettoriale (array, operazioni matematiche)
- `matplotlib.use('Agg')`: usa il backend non interattivo (per salvare file senza display)
- `matplotlib.pyplot`: plotting principale
- `matplotlib.patches`: forme geometriche personalizzate (es. rombi)
- `GridSpec`: layout avanzati per subplot multipli
- `scipy.special.expit, logit`: funzioni logistica e logit
- `warnings.filterwarnings('ignore')`: sopprime gli avvisi (es. divisioni per zero)

---

### Sezione 1: Dati

```python
gfap = {
    'study': [...],
    'TP': np.array([18, 24, 52, 38,  7, 33, 68, 99]),
    'FP': np.array([61,199,392,130, 89,425,477,285]),
    'FN': np.array([ 0,  1,  1,  1,  0,  0,  2,  1]),
    'TN': np.array([127,167,413,149, 34,221,663,373])
}
```

**Cosa fa:** Crea un dizionario con i dati GFAP.
- In Python, i dizionari sono simili alle liste nominate in R
- `np.array()`: crea array NumPy (più efficienti delle liste Python per calcoli)

---

### Funzione `enrich()`

```python
def enrich(d):
    d['N']   = d['TP'] + d['FP'] + d['FN'] + d['TN']
    d['Se']  = d['TP'] / (d['TP'] + d['FN'])
    d['Sp']  = d['TN'] / (d['TN'] + d['FP'])
    d['prev']= (d['TP'] + d['FN']) / d['N']
    d['CTp'] = d['TP'] + d['FN']
    return d
```

**Cosa fa:** Aggiunge al dizionario le variabili derivate (N, Se, Sp, prevalenza, CT+).
- Equivalente alle righe `gfap_data$N <- ...` in R
- `return d`: restituisce il dizionario modificato

---

### Sezione 2: Modello di Reitsma (Implementazione Python)

```python
def reitsma_pool(data, add_half=0.5):
    """
    Bivariate random-effects pooling in logit space.
    Returns: (se_pooled, sp_pooled, se_ci, sp_ci)
    """
```

**Cosa fa:** Definisce la funzione che implementa il modello bivariato.
- `add_half=0.5`: correzione di continuità per celle con zero (evita logit di 0)

```python
    tp = data['TP'].astype(float)
    fp = data['FP'].astype(float)
    fn = data['FN'].astype(float)
    tn = data['TN'].astype(float)
```

**Cosa fa:** Converte gli array in float per evitare problemi di divisione intera.

```python
    tp[tp == 0] += add_half
    fn[fn == 0] += add_half
    fp[fp == 0] += add_half
    tn[tn == 0] += add_half
```

**Cosa fa:** Aggiunge 0.5 alle celle con valore zero (correzione di Haldane-Anscombe).
- **Perché:** `logit(0)` e `logit(1)` sono indefiniti (−∞ e +∞)
- Aggiungere 0.5 è una convenzione standard nelle meta-analisi

```python
    se  = tp / (tp + fn)
    sp  = tn / (tn + fp)
    logit_se = np.log(se / (1 - se))
    logit_sp = np.log(sp / (1 - sp))
```

**Cosa fa:** Calcola sensibilità, specificità e le loro trasformate logit.

```python
    var_logit_se = 1/tp + 1/fn
    var_logit_sp = 1/tn + 1/fp
```

**Cosa fa:** Calcola le varianze delle stime logit (approssimazione delta).
- Deriva dalla varianza asintotica del logit: `Var(logit(p)) ≈ 1/(np) + 1/(n(1−p))`

```python
    w_se = 1 / var_logit_se
    Q_se = np.sum(w_se * (logit_se - np.sum(w_se*logit_se)/np.sum(w_se))**2)
    k    = len(tp)
    tau2_se = max(0, (Q_se - (k-1)) / (np.sum(w_se) - np.sum(w_se**2)/np.sum(w_se)))
```

**Cosa fa:** Stima la varianza inter-studio (τ²) con il metodo di **DerSimonian-Laird**.
- `w_se`: pesi inversi della varianza (studi più precisi = più peso)
- `Q_se`: statistica Q di Cochran (misura l'eterogeneità)
- `tau2_se`: τ² = varianza inter-studio (se Q > k−1, altrimenti 0)

**Formula di DerSimonian-Laird:**
```
τ² = max(0, (Q − (k−1)) / (Σw − Σw²/Σw))
```

```python
    w_re_se = 1 / (var_logit_se + tau2_se)
    mu_se   = np.sum(w_re_se * logit_se) / np.sum(w_re_se)
    se_mu_var = 1 / np.sum(w_re_se)
```

**Cosa fa:** Calcola la media pesata con pesi random-effects.
- `w_re_se`: pesi random-effects (varianza intra + inter-studio)
- `mu_se`: media pesata in scala logit
- `se_mu_var`: varianza della media pesata

```python
    se_ci_lo  = expit(mu_se - 1.96*np.sqrt(se_mu_var))
    se_ci_hi  = expit(mu_se + 1.96*np.sqrt(se_mu_var))
```

**Cosa fa:** Calcola l'intervallo di confidenza al 95% e lo ritrasforma con `expit()`.
- `1.96`: quantile 97.5% della normale standard (Z per IC 95%)

```python
    return {
        'Se': expit(mu_se), 'Sp': expit(mu_sp),
        'Se_CI': (se_ci_lo, se_ci_hi),
        'Sp_CI': (sp_ci_lo, sp_ci_hi),
        'tau2_Se': tau2_se, 'tau2_Sp': tau2_sp,
        'individual_Se': se, 'individual_Sp': sp,
        'var_se': var_logit_se, 'var_sp': var_logit_sp
    }
```

**Cosa fa:** Restituisce un dizionario con tutti i risultati.

---

### Sezione 3: Intervalli di Confidenza Wilson Score

```python
def wilson_ci(p, n, z=1.96):
    n = np.maximum(n, 1)
    centre = (p + z**2/(2*n)) / (1 + z**2/n)
    margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / (1 + z**2/n)
    return np.maximum(0, centre-margin), np.minimum(1, centre+margin)
```

**Cosa fa:** Calcola gli intervalli di confidenza con il metodo **Wilson score**.

**Perché Wilson e non Wald?**
- Il metodo classico (Wald): `p ± 1.96 × √(p(1−p)/n)` è impreciso per proporzioni vicine a 0 o 1
- Wilson score è più accurato, specialmente con campioni piccoli o proporzioni estreme

**Formula:**
```
centre = (p + z²/(2n)) / (1 + z²/n)
margin = z × √(p(1−p)/n + z²/(4n²)) / (1 + z²/n)
```

---

### Sezione 4: Forest Plot in Python

```python
def forest_plot(data, result, color_main, color_light, title, filename):
    k = len(data['study'])
    ...
    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, k*0.8+2)))
```

**Cosa fa:** Definisce una funzione per creare forest plot personalizzati.
- `plt.subplots(1, 2)`: crea una figura con 2 subplot affiancati
- `figsize=(16, max(6, k*0.8+2))`: dimensione adattata al numero di studi

```python
    sizes = np.sqrt(data['N']) / np.sqrt(data['N'].max()) * 200
    ax.scatter(vals, y_pos, s=sizes, color=color_light,
               edgecolors=color_main, linewidths=1.5, zorder=3)
```

**Cosa fa:** Disegna i quadrati dei singoli studi.
- `s=sizes`: area proporzionale a √N (studi più grandi = quadrati più grandi)
- `zorder=3`: livello di sovrapposizione (più alto = più in alto)

```python
    ax.errorbar(vals, y_pos, xerr=[xerr_lo, xerr_hi],
                fmt='none', color=color_main, capsize=4, lw=1.5, zorder=2)
```

**Cosa fa:** Aggiunge le barre di errore (intervalli di confidenza).
- `xerr=[xerr_lo, xerr_hi]`: errori asimmetrici a sinistra e destra
- `capsize=4`: dimensione delle tacche alle estremità

```python
    diamond = plt.Polygon(
        [[pooled-d_w, d_y], [pooled, d_y+0.35],
         [pooled+d_w, d_y], [pooled, d_y-0.35]],
        closed=True, facecolor=color_main, edgecolor='white', lw=1.5, zorder=4
    )
    ax.add_patch(diamond)
```

**Cosa fa:** Disegna il rombo pooled.
- `plt.Polygon()`: crea un poligono personalizzato (rombo)
- I 4 vertici definiscono la forma del rombo
- `ax.add_patch()`: aggiunge il rombo al grafico

---

### Sezione 5: Curva SROC in Python

```python
def sroc_curve(se_vals, sp_vals, var_se, var_sp, color, n_points=200):
    """Generate SROC curve points using Moses-Shapiro-Littenberg."""
    fpr = 1 - sp_vals
    D   = logit(se_vals) + logit(sp_vals)
    S   = logit(se_vals) - logit(sp_vals)
    w   = 1 / (var_se + var_sp)
```

**Cosa fa:** Definisce la funzione per generare la curva SROC con il metodo **Moses-Shapiro-Littenberg**.
- `D`: differenza diagnostica = logit(Se) + logit(Sp)
- `S`: soglia = logit(Se) − logit(Sp)
- `w`: pesi inversi della varianza combinata

```python
    b = np.sum(w * (S - S_bar) * (D - D_bar)) / np.sum(w * (S - S_bar)**2)
    a = D_bar - b * S_bar
```

**Cosa fa:** Regressione lineare pesata di D su S.
- `b`: coefficiente angolare (pendenza)
- `a`: intercetta

```python
    t_grid = np.linspace(-4, 4, n_points)
    Se_grid = expit((a + t_grid) / 2)
    Sp_grid = expit(-(a - t_grid) / 2)
    return 1 - Sp_grid, Se_grid
```

**Cosa fa:** Genera i punti della curva SROC.
- `t_grid`: griglia di valori da −4 a +4 (200 punti)
- `Se_grid, Sp_grid`: sensibilità e specificità predette per ogni punto della griglia

---

### Sezione 6: Clinical Summary Figure in Python

```python
fig = plt.figure(figsize=(18, 11))
gs  = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
```

**Cosa fa:** Crea una figura con layout a griglia 2×2.
- `GridSpec(2, 2)`: griglia 2 righe × 2 colonne
- `hspace, wspace`: spaziatura tra i subplot

```python
ax1 = fig.add_subplot(gs[0, 0])
```

**Cosa fa:** Aggiunge il primo subplot nella posizione (riga 0, colonna 0).

```python
bars1 = ax1.bar(biomarkers, se_vals_bar, color=colors, alpha=0.85,
                edgecolor='white', linewidth=2, width=0.5)
ax1.errorbar(range(2), se_vals_bar,
             yerr=[se_errs_lo, se_errs_hi],
             fmt='none', color='black', capsize=8, lw=2)
```

**Cosa fa:** Crea un grafico a barre con barre di errore.
- `ax1.bar()`: barre verticali
- `ax1.errorbar()`: barre di errore (IC 95%)

```python
tbl = ax4.table(
    cellText  = table_data[1:],
    colLabels = table_data[0],
    loc       = 'center',
    cellLoc   = 'center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)
tbl.scale(1, 2.0)
```

**Cosa fa:** Crea una tabella riassuntiva nel quarto pannello.
- `ax4.table()`: tabella matplotlib
- `tbl.scale(1, 2.0)`: allarga le celle verticalmente di 2×

---

## 6. Interpretazione dei Risultati {#interpretazione-risultati}

### Cosa ci dicono i risultati?

| Metrica | GFAP+UCH-L1 | S100B |
|---|---|---|
| Sensibilità | ~92-95% | ~95-98% |
| Specificità | ~35-45% | ~25-35% |
| NPV | ~98-99% | ~99% |
| CT evitate | ~35-45% | ~25-35% |
| Lesioni mancate/1000 | ~5-8 | ~2-5 |

### Come interpretare?

**GFAP+UCH-L1:**
- ✅ Buona sensibilità (vede la maggior parte delle lesioni)
- ✅ NPV eccellente: se il test è negativo, possiamo essere molto sicuri che non ci siano lesioni
- ⚠️ Specificità moderata: molti falsi positivi (pazienti sani che risultano positivi)
- **Implicazione clinica:** Utile per **escludere** la TAC (rule-out test), ma non per confermare la necessità di TAC

**S100B:**
- ✅ Sensibilità molto alta (quasi nessuna lesione mancata)
- ✅ NPV quasi perfetto
- ⚠️ Specificità bassa: ancora più falsi positivi di GFAP
- **Implicazione clinica:** Eccellente come test di **screening** per escludere la TAC, ma molti pazienti verrebbero irradiati inutilmente

### Quale scegliere?

Dipende dall'obiettivo clinico:
- **Priorità: non mancare lesioni** → S100B (sensibilità più alta)
- **Priorità: ridurre TAC inutili** → GFAP+UCH-L1 (specificità più alta)

---

## 7. Glossario dei Termini Statistici {#glossario}

| Termine | Definizione |
|---|---|
| **TP (True Positive)** | Paziente con lesione correttamente identificato dal test |
| **FP (False Positive)** | Paziente senza lesione erroneamente identificato come positivo |
| **FN (False Negative)** | Paziente con lesione non identificato dal test |
| **TN (True Negative)** | Paziente senza lesione correttamente identificato come negativo |
| **Sensibilità** | Probabilità che il test sia positivo dato che il paziente ha la lesione: P(T+|D+) |
| **Specificità** | Probabilità che il test sia negativo dato che il paziente non ha la lesione: P(T−|D−) |
| **NPV** | Probabilità che il paziente non abbia la lesione dato che il test è negativo: P(D−|T−) |
| **Prevalenza** | Proporzione di pazienti con la lesione nella popolazione studiata |
| **Logit** | Trasformazione: logit(p) = ln(p/(1−p)). Mappa (0,1) su (−∞,+∞) |
| **Expit** | Funzione inversa del logit: expit(x) = 1/(1+exp(−x)). Mappa (−∞,+∞) su (0,1) |
| **Eterogeneità (τ²)** | Varianza inter-studio: quanto i risultati dei studi variano oltre l'errore campionario |
| **Modello ad Effetti Random** | Modello che assume che ogni studio abbia un "vero effetto" diverso, distribuito attorno a una media comune |
| **DerSimonian-Laird** | Metodo per stimare τ² basato sulla statistica Q di Cochran |
| **Forest Plot** | Grafico che mostra le stime di ogni studio con i loro IC e la stima pooled |
| **SROC** | Curva che riassume la relazione tra sensibilità e specificità attraverso gli studi |
| **Wilson Score** | Metodo per calcolare IC di proporzioni più accurato del metodo Wald classico |
| **Correzione di Haldane-Anscombe** | Aggiungere 0.5 alle celle con zero per evitare logit indefiniti |
| **Sensitivity Analysis** | Verifica della robustezza dei risultati escludendo studi uno alla volta |

---

## Riassunto Finale

### Cosa abbiamo fatto?

1. **Raccolto i dati** da 14 studi (8 per GFAP+UCH-L1, 6 per S100B)
2. **Stimato sensibilità e specificità pooled** con il modello bivariato di Reitsma
3. **Calcolato il NPV** usando il teorema di Bayes
4. **Stimato l'impatto clinico**: lesioni mancate e TAC evitate per 1000 pazienti
5. **Visualizzato i risultati** con forest plot, curve SROC e grafici comparativi
6. **Verificato la robustezza** con sensitivity analysis

### Perché due implementazioni (R e Python)?

- **R**: usa pacchetti specializzati (`mada`) con il modello di Reitsma completo (massima verosimiglianza, correlazione Se-Sp)
- **Python**: implementazione manuale con DerSimonian-Laird (più trasparente, utile per didattica)
- I risultati dovrebbero essere **simili** ma non identici (metodi di stima leggermente diversi)

### Messaggio chiave per i clinici

> **Entrambi i biomarcatori hanno NPV eccellente (>98%), il che li rende utili come test di "rule-out" per evitare TAC inutili nei pazienti con mTBI. La scelta dipende dal bilanciamento desiderato tra sicurezza (non mancare lesioni) ed efficienza (ridurre esami non necessari).**

---

*Documento creato come materiale didattico per la comprensione della meta-analisi di test diagnostici in ambito biomedico.*
