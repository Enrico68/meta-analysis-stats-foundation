# Guida alla Meta-Analisi Diagnostica per Biomarcatori mTBI
## Spiegazione semplice del codice e dei concetti statistici

---

## Indice
1. Introduzione: cosa e una meta-analisi
2. I dati: la tabella di contingenza 2x2
3. Il codice R spiegato riga per riga
4. Il codice Python spiegato
5. Spiegazione dei grafici
6. Concetti statistici fondamentali

---

## 1. Introduzione: cosa e una meta-analisi

### L idea di base

Immagina di voler sapere se un nuovo test medico funziona bene. Un singolo studio potrebbe avere pochi pazienti o risultati sfortunati. La **meta-analisi** e come fare una media intelligente di molti studi diversi per ottenere una risposta piu affidabile.

**Metafora**: Se chiedi a 3 amici se un ristorante e buono, potresti avere opinioni diverse. Se leggi 100 recensioni su TripAdvisor, hai un idea molto piu precisa!

### Meta-analisi diagnostica

In questo progetto analizziamo **biomarcatori** (sostanze nel sangue) per diagnosticare il **trauma cranico lieve (mTBI)**. I due biomarcatori sono:
- **GFAP + UCH-L1**: Proteine cerebrali rilasciate dopo un trauma
- **S100B**: Un altra proteina di origine cerebrale

---

## 2. I dati: la tabella di contingenza 2x2

### La tabella della verita

Ogni studio confronta il test con la verita (es. TC o risonanza magnetica). I risultati si organizzano in 4 caselle:

```
                         RISULTATO REALE
                    (secondo la TC/Risonanza)
                    +------------------------+
                    |  Positivo |  Negativo  |
     +--------------+-----------+------------+
     |   Positivo   |    TP     |     FP     |
TEST |   (test +)   |  (veri +) | (falsi +)  |
     +--------------+-----------+------------+
     |   Negativo   |    FN     |     TN     |
     |   (test -)   |  (falsi -)|  (veri -)  |
     +--------------+-----------+------------+
```

### Significato delle 4 caselle

| Simbolo | Nome | Significato | Esempio |
|---------|------|-------------|---------|
| **TP** | True Positives | Test dice SI e il paziente ha davvero il problema | Il test rileva il trauma e la TC conferma |
| **TN** | True Negatives | Test dice NO e il paziente e davvero sano | Il test e negativo e la TC e negativa |
| **FP** | False Positives | Test dice SI ma il paziente e sano (falso allarme) | Il test e positivo ma la TC e negativa |
| **FN** | False Negatives | Test dice NO ma il paziente ha il problema (pericoloso!) | Il test e negativo ma la TC mostra un trauma |

### Esempio pratico dai dati

```r
# Studio Bazarian 2021 per GFAP+UCH-L1
TP = 115   # 115 pazienti: test +, TC +
FP = 1061  # 1061 pazienti: test +, TC - (falsi allarmi)
FN = 5     # 5 pazienti: test -, TC + (perditi!)
TN = 720   # 720 pazienti: test -, TC -
```

---

## 3. Il codice R spiegato riga per riga

### 3.1 Caricamento delle librerie

```r
if (!require(mada)) install.packages(mada, repos=https://cran.r-project.org)
if (!require(metafor)) install.packages(metafor, repos=https://cran.r-project.org)
library(mada); library(metafor)
```

**Spiegazione**: 
- `mada` = Meta-Analysis of Diagnostic Accuracy (meta-analisi per test diagnostici)
- `metafor` = Meta-Analysis FOR (libreria generale per meta-analisi)
- Il codice controlla se sono installate, e se no le installa automaticamente

### 3.2 Inserimento dei dati

```r
gfap <- data.frame(
  study=c(Bazarian 2021,Lagares 2024,Faisal 2023,Legramante 2024,
          Ladang 2025,Milevoj K. 2025,Curran 2025,Campagna 2025),
  TP=c(115,176,12,7,112,107,12,8), 
  FP=c(1061,946,75,89,148,529,56,59),
  FN=c(5,3,0,0,1,5,0,0), 
  TN=c(720,313,43,34,101,181,21,30))
gfap$N <- gfap$TP+gfap$FP+gfap$FN+gfap$TN
```

**Spiegazione**:
- Si crea una tabella con i nomi degli 8 studi su GFAP+UCH-L1
- Ogni riga e uno studio con i suoi 4 numeri (TP, FP, FN, TN)
- `gfap$N` calcola il numero totale di pazienti per studio

### 3.3 Il modello di Reitsma (il cuore dell analisi)

```r
fit_gfap <- reitsma(gfap, correction.control=all, correction=0.5)
```

**Cosa fa**: Il **modello bivariato di Reitsma** e il metodo standard d oro per meta-analisi diagnostiche. 

**Perche bivariato?** Perche analizza contemporaneamente:
1. **Sensibilita** (quanto e bravo a trovare i malati)
2. **Specificita** (quanto e bravo a escludere i sani)

**La correzione di 0.5**: Quando uno studio ha FN=0 o TP=0, alcune formule si rompono (divisione per zero). Si aggiunge 0.5 a tutte le caselle per evitare problemi matematici.

### 3.4 Metriche cliniche calcolate

```r
npv <- sum(gfap$TN)/(sum(gfap$TN)+sum(gfap$FN))
ct_red <- (sum(gfap$FN)+sum(gfap$TN))/sum(gfap$N)*100
```

| Metrica | Formula | Significato clinico |
|---------|---------|---------------------|
| **NPV** | TN / (TN + FN) | Se il test e negativo, qual e la probabilita che il paziente sia davvero sano? |
| **CT-red** | (FN + TN) / N x 100 | Percentuale di TC che si possono evitare usando il test |

**Esempio**: NPV = 99% significa che se il test e negativo, c e il 99% di probabilita che il paziente NON abbia un trauma (quindi la TC non serve).

### 3.5 Il Forest Plot (grafico a foresta)

La funzione `make_forest` crea il grafico piu importante. Per ogni studio mostra:
- Un **quadrato**: la stima puntuale (sensibilita o specificita)
- Una **linea orizzontale**: l intervallo di confidenza (margine di errore)
- Un **rombo rosso**: la stima pooled (combinata di tutti gli studi)

### 3.6 Il test di Deeks (funnel plot asimmetria)

```r
gfap$DOR <- ((gfap$TP+0.5)*(gfap$TN+0.5))/((gfap$FP+0.5)*(gfap$FN+0.5))
gfap$lnDOR <- log(gfap$DOR)
deeks <- lm(lnDOR~inv_sqrt_ESS, data=gfap, weights=1/var_lnDOR)
```

**Cosa e il DOR?** Diagnostic Odds Ratio = (TPxTN)/(FPxFN). Misura quanto e bravo il test a discriminare malati da sani.

**Il test di Deeks**: Controlla se ci sono studi "mancanti" (publication bias). Se p > 0.10, non c e asimmetria.

---

## 4. Il codice Python spiegato

### 4.1 Struttura generale

Il file `mTBI_DTA_Plots_v2.py` crea grafici professionali partendo dai risultati della meta-analisi fatta in R.

### 4.2 Intervallo di confidenza di Wilson

```python
def wilson_ci(x, n, confidence=0.95):
    z = stats.norm.ppf(1 - (1 - confidence) / 2)  # 1.96 per 95%
    p = x / n  # proporzione osservata
    # Formula complessa che da un intervallo piu accurato
```

**Cosa e?** L intervallo di confidenza dice: "con il 95% di probabilita, il valore vero sta tra X e Y".

**Esempio**: Se uno studio ha sensibilita 95% con IC 90-98%, significa che il valore vero e probabilmente tra 90% e 98%.

### 4.3 Sensibilita e Specificita

```python
def calculate_sens_spec(data, correction=0.5):
    sens = tp / (tp + fn)  # Quanti malati ho trovato?
    spec = tn / (tn + fp)  # Quanti sani ho escluso?
```

| Misura | Formula | Significato |
|--------|---------|-------------|
| **Sensibilita** | TP / (TP + FN) | Percentuale di malati correttamente identificati |
| **Specificita** | TN / (TN + FP) | Percentuale di sani correttamente esclusi |

**Esempio**: Sensibilita 96% = su 100 malati, il test ne trova 96 (4 li perde!).

---

## 5. Spiegazione dei grafici

### 5.1 Forest Plot (forest_GFAP_UCH_L1_python.png)

Il Forest Plot e il grafico PIU IMPORTANTE della meta-analisi.

#### Come leggerlo:

```
SENSITIVITA (sinistra)              SPECIFICITA (destra)
  0.70    0.85    1.00               0.0     0.3     0.5
    |       |       |                  |       |       |
Campagna 2025 -----[==]---- 94.4%      -----[==]---- 33.9%
Curran 2025    --[====]--- 96.2%       ---[==]----- 27.6%
   ...
Pooled           [<>]      98.1%         [<>]      31.9%
                   ^                         ^
              (rombo rosso)             (rombo rosso)
```

#### Elementi del grafico:

| Elemento | Forma | Colore | Significato |
|----------|-------|--------|-------------|
| Stima singolo studio | Quadrato | Blu | Sensibilita/Specificita di QUELLO studio |
| Intervallo confidenza | Linea orizzontale | Blu | Margine di errore dello studio |
| Stima pooled (combinata) | Rombo | Rosso | La media di TUTTI gli studi |
| Linea verticale | Tratteggiata | Rossa | Il valore pooled di riferimento |

#### Cosa ci dice questo grafico:

**Sensibilita (pooled: 98.1%)**:
- Il test e OTTIMO a trovare i pazienti con trauma
- Tutti gli studi sono concordi (linee si sovrappongono)
- Quasi nessun falso negativo (molto sicuro!)

**Specificita (pooled: 31.9%)**:
- Il test e SCARSO a escludere i sani
- Molti falsi positivi (falsi allarmi)
- Questo e voluto: meglio un falso allarme che un paziente perduto!

---

### 5.2 Deeks Funnel Plot (deeks_funnel_GFAP_python.png)

Il Funnel Plot serve per controllare se ci sono studi "mancanti" (publication bias).

#### Come leggerlo:

```
ln(DOR) - Odds Ratio Diagnostico
  4.0 |                    Ladang 2025
      |                         o
  3.0 |
      |  o Lagares 2024
  2.5 | o Bazarian 2021
      |      o Milevoj K. 2025
      |                 |              2.0 |             \  <-- linea di regressione
      |                    |               o Faisal 2023
      |        o Legramante 2024
      |                   o Curran 2025
      |                    o Campagna 2025
      +------------------------------------
      0.01   0.02   0.03   0.04   0.05
           1/sqrt(Effective Sample Size)
```

#### Elementi del grafico:

| Elemento | Forma | Significato |
|----------|-------|-------------|
| Punti | Cerchi blu | Ogni studio, posizionato in base alla precisione |
| Linea | Tratteggiata rossa | Regressione lineare (trend) |
| Asse X | 1/sqrt(ESS) | Precisione dello studio (piu piccolo = piu preciso) |
| Asse Y | ln(DOR) | Accuratezza diagnostica dello studio |

#### Interpretazione:

- **p = 0.431** (dal titolo): Molto maggiore di 0.10
- **Conclusione**: NON c e asimmetria significativa
- **Significato**: Non sembra che studi "scomodi" siano stati nascosti

**Metafora**: E come controllare se in una classe i voti bassi sono stati nascosti. Se i voti sono distribuiti simmetricamente, probabilmente tutti sono stati mostrati!

---

### 5.3 SROC Plot (SROC_S100B_python.png)

SROC = Summary Receiver Operating Characteristic

#### Come leggerlo:

```
Sensitivity (True Positive Rate)
  1.0 |                    o Oris (2024)
      |                   D Summary (Se=0.98, Sp=0.19)
  0.9 |              o Oris (2021)
      |                o Hopman (2023)
  0.8 |
      |
  ... |
      |
  0.0 +------------------------------------
      0.0    0.2    0.4    0.6    0.8    1.0
           False Positive Rate (1 - Specificity)
```

#### Elementi del grafico:

| Elemento | Forma | Colore | Significato |
|----------|-------|--------|-------------|
| Studi individuali | Cerchi | Verde | Posizione di ogni studio nello spazio Se-Sp |
| Stima pooled | Rombo | Verde scuro | La media combinata di tutti gli studi |
| Angolo in alto a sinistra | - | - | Zona IDEALE (alta Se, bassa FPR) |

#### Interpretazione per S100B:

- **Sensibilita = 0.98 (98%)**: Ottima! Trova quasi tutti i traumi
- **Specificita = 0.19 (19%)**: Molto bassa! Molti falsi positivi
- **Posizione**: In alto a destra = buona sensibilita, ma scarsa specificita

#### Confronto ideale:

```
        Specificita alta
              ^
              |
   Specifica  |    IDEALE
   ma poco    |    (alto sinistra)
   sensibile  |
   (basso     |
   sinistra)  |
   <----------+----------> Specificita bassa
              |           ma molto
              |           sensibile
              |           (basso destra)
              v
        Sensibilita bassa
```

---

## 6. Concetti statistici fondamentali

### 6.1 Sensibilita vs Specificita: il compromesso

```
              | Test + | Test -
--------------+--------+--------
Ha il trauma  |   TP   |   FN   <- Sensibilita = TP/(TP+FN)
              |        |
Non ha trauma |   FP   |   TN   <- Specificita = TN/(TN+FP)
              |        |
```

**Compromesso**: Se alzo il "livello di allarme" del test:
- Aumenta la specificita (meno falsi positivi)
- Ma diminuisce la sensibilita (piu falsi negativi - PERICOLOSO!)

Per il trauma cranico, preferiamo **alta sensibilita** anche a costo di falsi allarmi.

### 6.2 Intervallo di confidenza (IC)

**Definizione semplice**: "Siamo sicuri al 95% che il valore vero sia tra X e Y"

**Esempio**: 
- Stima: Sensibilita = 96%
- IC 95%: 94% - 98%
- Significato: Il valore vero e probabilmente tra 94% e 98%

**Regola pratica**: Se l'IC e stretto, lo studio e preciso. Se e largo, c e incertezza.

### 6.3 Modello bivariato di Reitsma

**Perche bivariato?** Perche sensibilita e specificita sono CORRELATE:
- Studi con pazienti piu gravi tendono ad avere ALTA sensibilita e BASSA specificita
- Studi con pazienti meno gravi tendono ad avere BASSA sensibilita e ALTA specificita

Il modello di Reitsma tiene conto di questa correlazione.

### 6.4 P-value e significativita

| Valore p | Interpretazione |
|----------|-----------------|
| p < 0.05 | Significativo (qualcosa sta succedendo) |
| p > 0.10 | Non significativo (probabilmente casuale) |
| 0.05-0.10 | Zona grigia |

**Nel test di Deeks**: p = 0.431 > 0.10 = Nessun bias di pubblicazione evidente

### 6.5 Diagnostic Odds Ratio (DOR)

```
        TP x TN      (veri positivi x veri negativi)
DOR = ----------- = ----------------------------------
        FP x FN      (falsi positivi x falsi negativi)
```

**Interpretazione**:
- DOR = 1: Il test e inutile (come lanciare una moneta)
- DOR = 10: Buon test
- DOR = 100: Ottimo test

**Esempio**: DOR = 50 significa che e 50 volte piu probabile avere un risultato corretto che sbagliato.

---

## Riepilogo finale

### Cosa abbiamo imparato:

1. **I dati**: Ogni studio fornisce 4 numeri (TP, FP, FN, TN)

2. **Il modello**: Reitsma combina tutti gli studi tenendo conto della correlazione Se-Sp

3. **Le metriche**:
   - Sensibilita: quanti malati trovo (vogliamo ALTA)
   - Specificita: quanti sani escludo (accettiamo BASSA per sicurezza)
   - NPV: probabilita di essere sani se test negativo
   - CT-red: quante TC evitiamo

4. **I grafici**:
   - Forest Plot: mostra risultati di ogni studio + stima combinata
   - Funnel Plot: controlla se ci sono studi mancanti
   - SROC: mostra il compromesso Se-Sp nello spazio bidimensionale

5. **Il risultato pratico**: 
   - GFAP+UCH-L1: Se=98.1%, Sp=31.9% 
   - Ottimo per ESCLUDERE traumi (test negativo = niente TC)
   - Molti falsi positivi, ma questo e accettabile in emergenza

---

*Questa guida e pensata per studenti e ricercatori che si avvicinano per la prima volta alla meta-analisi diagnostica. Per approfondimenti, consultare la letteratura specialistica sui modelli bivariati di Reitsma.*
