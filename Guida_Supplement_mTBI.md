# 📗 Guida al Supplemento del Manoscritto mTBI
## "Diagnostic Accuracy of Blood Biomarkers GFAP, UCH-L1, and S100B"

---

## 📋 **PAGINA 1: eMethods - Detailed Statistical Analysis**

### Cos'è una sezione eMethods?
È la versione "estesa" dei metodi. Nel manoscritto principale i metodi sono sintetici, qui troviamo **tutti i dettagli tecnici** per chi vuole replicare l'analisi.

---

### 🔬 Meta-analytic Model (Modello Meta-analitico)

**Cosa dice:**
> *"Pooled sensitivity and specificity were estimated separately using univariate random-effects models with restricted maximum likelihood (REML) estimation"*

**Traduzione semplice:**
- I ricercatori hanno calcolato una **media ponderata** di sensibilità e specificità dai diversi studi
- "Random-effects" significa che i risultati dei singoli studi possono variare (non sono identici)
- "REML" è un metodo statistico sofisticato per stimare questa variazione

**Perché REML e non semplice media?**
- I singoli studi hanno dimensioni diverse (da 130 a 1.901 pazienti)
- REML dà più peso agli studi più grandi (più precisi)
- È il metodo standard per meta-analisi diagnostiche

---

### 🔄 Logit Transformation (Trasformazione Logit)

**Il problema:** Le proporzioni (tipo 96,7%) non sono distribuite normalmente (una campana di Gauss).
- Non possono essere < 0% o > 100%
- Sono "schiacciate" ai due estremi

**La soluzione - Formula del Logit:**
```
logit(p) = ln(p / (1-p))

Esempio:
Se p = 0,967 (96,7%)
logit(0,967) = ln(0,967 / 0,033) = ln(29,3) ≈ 3,38
```

**Cosa fa R:**
1. Trasforma tutte le sensibilità in logit
2. Fa la media pesata su scala logit
3. Riconverte indietro con **inverse logit** (plogis in R)

**Perché è importante?**
Senza questa trasformazione, la media di 90% e 99% darebbe un risultato sbagliato (94,5% invece del corretto ~96%).

---

### 🛠️ Continuity Correction (Correzione di Continuità)

**Il problema:**
- Se in uno studio TP = 0 (nessun vero positivo), la sensibilità è 0/qualcosa = 0
- logit(0) = -∞ (meno infinito!) → impossibile da calcolare

**La soluzione:**
> *"A continuity correction of 0.5 was added to all four cells of that study before logit transformation"*

**Esempio pratico:**
```
Studio Papa 2022: TP=23, FN=0, TN=81, FP=245
Prima della correzione: Se = 23/(23+0) = 100% → logit = +∞
Dopo la correzione:   TP=23,5, FN=0,5 → Se = 23,5/24 = 97,9% → logit calcolabile!
```

**Nota:** Questa correzione introduce un **piccolo bias**, ma è necessaria. Per questo Papa 2022 risulta Se = 100,0% nel manoscritto ma ha valori leggermente diversi nei calcoli interni.

---

### 📊 Heterogeneity Assessment (Valutazione Eterogeneità)

**Tre misure usate:**

| Misura | Cosa misura | Interpretazione |
|--------|-------------|-----------------|
| **Q (Cochran)** | Se la variabilità è casuale | p < 0,10 = c'è eterogeneità reale |
| **I²** | % di variabilità "vera" vs casuale | 0-25% bassa, 25-50% moderata, 50-75% sostanziale, >75% considerevole |
| **τ² (tau-squared)** | Varianza tra studi su scala logit | Più alto = più diversi tra loro |

**Nel nostro studio:**
```
Sensibilità GFAP:  I² = 0%, τ² = 0,000, Q p = 0,77
→ Tutti gli studi concordano! Nessuna eterogeneità.

Specificità GFAP:  I² = 93%, τ² = 0,093, Q p < 0,001
→ Gli studi discordano molto! L'età spiega questa variabilità.
```

---

### 📈 Summary ROC Curves (Curve ROC Riassuntive)

**Cosa sono:**
Curve che mostrano il trade-off tra sensibilità e specificità.

**Modello Bivariato vs Univariato:**
| Modello | Descrizione | Uso |
|---------|-------------|-----|
| **Univariato** (principale) | Analizza Se e Sp separatamente | Stime pooled di Se e Sp |
| **Bivariato** (SROC) | Analizza Se e Sp insieme, considerando la loro correlazione | Curve ROC e ellissi di confidenza |

**Perché la correlazione è importante?**
- Negli studi diagnostici, Se e Sp sono spesso correlati
- Se uno studio ha metodologia più rigorosa, tende ad avere sia Se più alta che Sp più alta
- Il modello bivariato tiene conto di questo

**Funzione R:** `mada::reitsma()` implementa il modello di Reitsma (modello bivariato standard per meta-analisi diagnostiche).

---

### 🏥 Clinical Utility Metrics (Metriche di Utilità Clinica)

**Formule utilizzate:**

**1. NPV (Negative Predictive Value):**
```
NPV = [Sp × (1 - prev)] / [Sp × (1 - prev) + (1 - Se) × prev]

Dove:
- Sp = specificità pooled
- Se = sensibilità pooled  
- prev = prevalenza CT+ (9,3% per GFAP, 7,8% per S100B)

Esempio GFAP:
NPV = [0,271 × (1 - 0,093)] / [0,271 × 0,907 + (1 - 0,967) × 0,093]
    = 0,246 / (0,246 + 0,003)
    = 0,246 / 0,249
    = 0,9875 = 98,75%
```

**2. CT Scans Avoided (TC Evitate):**
```
% CT avoided = Sp × (1 - prev) × 100

Esempio GFAP:
0,271 × 0,907 × 100 = 24,6%
```
Questa formula dice: quanti dei pazienti sani (1-prev) vengono correttamente identificati come tali (Sp)?

**3. Missed Lesions (Lesioni Perse):**
```
Missed per 1.000 = (1 - Se) × prev × 1.000

Esempio GFAP:
(1 - 0,967) × 0,093 × 1.000 = 0,033 × 0,093 × 1.000 = 3,1
```
Questo è il numero di lesioni che "passano" attraverso il test perché falsi negativi.

---

### 🔍 Publication Bias (Bias di Pubblicazione)

**Test di Egger:**
```
Regressione: Standardized Effect vs Precision

Standardized Effect = (logit Se dello studio - logit Se pooled) / SE dello studio
Precision = 1 / SE dello studio
```

**Interpretazione:**
- Se l'intercetta è significativamente diversa da 0 (p < 0,05) → c'è asimmetria (possible bias)
- Se p > 0,05 → nessuna evidenza di bias

**Nel nostro studio:**
```
GFAP+UCH-L1: p = 0,430 → ✅ No bias
S100B:        p = 0,224 → ✅ No bias
```

**⚠️ Limitazione importante:**
> *"Given the small number of studies (k=7 and k=4), the power of this test is limited"*

Con pochi studi, il test ha bassa potenza statistica (alto rischio di falsi negativi). Per questo i funnel plot sono importanti (valutazione visiva).

---

### ✅ Cross-validation (Convalida Incrociata)

**Cosa hanno fatto:**
- Replicato TUTTE le analisi in Python indipendentemente da R
- Metodo diverso: DerSimonian-Laird (DL) invece di REML

**Risultato:**
> *"Point estimates were concordant within 0.1 percentage points"*

**Perché è importante:**
- Conferma che i risultati non dipendono dal software usato
- DerSimonian-Laird è più semplice ma meno efficiente di REML
- Se i due metodi danno lo stesso risultato, possiamo fidarci

---

## 📄 **PAGINA 2: eTable 1 - Complete Search Strategies**

### Perché mostrare le query di ricerca?
Perché la **replicabilità** è fondamentale in scienza. Chiunque deve poter rifare la stessa ricerca e trovare gli stessi articoli.

### Database utilizzati:

| Database | Query A (GFAP) | Query B (S100B) | Query C (Combined) | Totale |
|----------|---------------|-----------------|-------------------|--------|
| **PubMed** | 109 | 122 | - | 231 |
| **Scopus** | 149 | 169 | 177 | 495 |
| **Cochrane** | - | - | 8 | 8 |
| **Semantic Scholar** | 39 | 16 | - | 55 |
| **Reference lists** | - | - | 2 | 2 |

**Totale:** 1.131 record
**Dopo rimozione duplicati:** 420 record unici

### Strategie di ricerca booleane:

**Query A (GFAP/UCH-L1):**
```
(GFAP OR "glial fibrillary acidic protein" OR UCH-L1 OR 
"ubiquitin C-terminal hydrolase") 
AND 
("mild traumatic brain injury" OR mTBI OR "minor head injury" OR 
concussion OR "minor head trauma") 
AND 
("computed tomography" OR "CT scan" OR intracranial OR "rule out" 
OR "diagnostic accuracy")
```

**Spiegazione:**
- **OR** tra sinonimi (trova almeno uno)
- **AND** tra concetti diversi (deve trovarli tutti)
- Parole chiave: biomarcatore + condizione clinica + outcome

**Query B (S100B):**
Simile ma con termini specifici per S100B.

**Query C (Scopus - Combined):**
Query più ampia per catturare studi che menzionano entrambi i biomarcatori.

---

## 📄 **PAGINA 3-4: eTable 2 - Extracted 2×2 Contingency Tables**

### Cosa sono le tabelle 2×2?
Sono la "matrice di confusione" di ogni studio:

```
                    TC+ (Lesione presente)    TC- (Lesione assente)
                    ┌─────────────────────┬─────────────────────┐
Test POSITIVO       │       TP            │        FP           │
                    │  (Vero Positivo)    │   (Falso Positivo)  │
                    ├─────────────────────┼─────────────────────┤
Test NEGATIVO       │       FN            │        TN           │
                    │  (Falso Negativo)   │   (Vero Negativo)   │
                    └─────────────────────┴─────────────────────┘
```

### Analisi degli studi:

#### **GFAP+UCH-L1 (7 studi - Analisi Primaria):**

| Studio | N | TP | FP | FN | TN | Se% | Sp% | Note |
|--------|---|----|----|----|----|-----|-----|------|
| **Chayoua 2024** | 253 | 57 | 157 | 2 | 37 | 96,6% | 19,1% | Specificità più bassa |
| **Bazarian 2021** | 1.901 | 115 | 1.061 | 5 | 720 | 95,8% | 40,4% | Studio più grande, Sp più alta |
| **Papa 2022** | 349 | 23 | 245 | 0 | 81 | 100,0% | 24,8% | FN=0 → Se=100% |
| **Lagares 2024** | 1.438 | 176 | 945 | 3 | 314 | 98,3% | 24,9% | Studio europeo |
| **Puravet 2026** | 1.010 | 50 | 700 | 0 | 260 | 100,0% | 27,1% | FN=0 → Se=100% |
| **Milevoj 2025** | 822 | 107 | 529 | 5 | 181 | 95,5% | 25,5% | Popolazione anziana |
| **Legramante 2024** | 130 | 7 | 89 | 0 | 34 | 100,0% | 27,6% | Studio più piccolo |

**Totale GFAP:** N=5.903, CT+=535 (9,3%)

#### **S100B (4 studi - Analisi Primaria):**

| Studio | N | TP | FP | FN | TN | Se% | Sp% | Note |
|--------|---|----|----|----|----|-----|-----|------|
| **Puravet 2026** | 1.010 | 48 | 720 | 2 | 240 | 96,0% | 25,0% | Stessa coorte di GFAP |
| **Seidenfaden 2021** | 566 | 32 | 366 | 0 | 168 | 100,0% | 31,5% | Campione in ospedale |
| **Rogan 2023** | 133 | 15 | 81 | 1 | 36 | 93,8% | 30,8% | 1 falso negativo |
| **Hopman 2023** | 495 | 69 | 355 | 5 | 66 | 93,2% | 15,7% | Dati da abstract |

**Totale S100B:** N=2.204, CT+=164 (7,8%)

#### **Studi Analisi di Sensitività (con asterisco \*):**

**GFAP (3 studi):**
- **Ladang 2025**: N=362, Se=99,1%, Sp=40,6% → prevalenza 31%! (pazienti anticoagulati)
- **Spaziani 2026**: N=441, Se=96,0%, Sp=24,3%
- **Osmić-Husni 2026**: N=102, Se=100%, Sp=32,8%

**S100B (3 studi):**
- **Jones 2020**: N=679, Se=84,6% → sensibilità più bassa, popolazione diversa
- **Oris 2021**: N=1.172, Se=100%, Sp=19,7% → coorte sovrapposta (esclusa)
- **Asadollahi 2015**: N=158, Se=94,9%, Sp=35,4%

**⚠️ Perché questi studi sono nell'analisi di sensitività?**
1. **Prevalenza molto diversa** (17-50% vs 8-9%)
2. **Popolazioni selezionate** (solo anticoagulati, solo neurochirurgia)
3. **Risk of bias elevato** (vedi eTable 3)

---

## 📄 **PAGINA 5: eTable 3 - QUADAS-2 Risk of Bias Assessment**

### Cos'è QUADAS-2?
È lo strumento standard per valutare la qualità degli studi diagnostici. Ha 4 domini:

| Dominio | Cosa valuta | Domanda chiave |
|---------|-------------|----------------|
| **D1: Patient Selection** | Selezione pazienti | I pazienti sono rappresentativi? |
| **D2: Index Test** | Test in esame | Il test è stato eseguito correttamente? |
| **D3: Reference Standard** | Gold standard | La TC è stata fatta/interpretata correttamente? |
| **D4: Flow & Timing** | Flusso e tempi | Tutti hanno ricevuto entrambi i test? |

### Legenda colori:
- 🟢 **Low (Basso)** - Nessun problema
- 🟡 **Unclear (Non chiaro)** - Informazioni insufficienti
- 🔴 **High (Alto)** - Problema serio

### Analisi risultati:

#### **Studi con Risk of Basso (tutti i domini Low):**
- Chayoua 2024
- Bazarian 2021
- Lagares 2024
- Puravet 2026
- Milevoj 2025
- Seidenfaden 2021

#### **Studi con Problemi:**

| Studio | Problema | Perché |
|--------|----------|--------|
| **Papa 2022** | D4: Unclear | Non tutti i pazienti negativi hanno fatto TC |
| **Legramante 2024** | D1: High | Popolazione ristretta (solo pazienti anziani con deficit cognitivi) |
| **Rogan 2023** | D4: Unclear | Tempistiche non chiare |
| **Hopman 2023** | D1, D4: Unclear | Dati solo da abstract, popolazione non chiara |

#### **Studi Analisi di Sensitività (tutti con problemi):**
- **Ladang, Spaziani**: D1 High (popolazioni selezionate), Applicability High
- **Osmić-Husni**: OK per RoB ma prevalenza alta
- **Jones 2020**: D2 High (test non standard), Applicability High
- **Oris 2021**: D4 High (follow-up non appropriato)
- **Asadollahi 2015**: D1, D2 High (retrospettivo, test non standard)

### Colonna "Applicability Concerns":
Valuta se lo studio è applicabile alla pratica clinica reale:
- **Low**: Sì, rappresentativo
- **High**: No, popolazione troppo selezionata (es: solo anticoagulati, solo neurochirurgia)

---

## 📄 **PAGINA 6: eTable 4 - Sensitivity Analysis Results**

### Cosa sono le analisi di sensitività?
Sono analisi "Cosa succederebbe se..." per testare la robustezza dei risultati.

### Tipi di analisi di sensitività presentate:

#### **1. Espansione del pool GFAP (+3 studi):**
```
Primary pool:     7 studi, N=5.903, Se=96,7%, Sp=27,1%
Expanded pool:   10 studi, N=6,808, Se=?, Sp=?
```
**Domanda:** I risultati cambiano se includiamo anche i 3 studi "borderline"?

#### **2. Espansione del pool S100B (+3 studi):**
```
Primary pool:     4 studi, N=2.204, Se=94,6%, Sp=25,0%
Expanded pool:    7 studi, N=4.213, Se=?, Sp=?
```

#### **3. Analisi sottogruppo Abbott (già presentata nel manoscritto):**
```
Abbott only:      5 studi, N=4.116, Se=96,0%, Sp=27,8%
```
**Domanda:** I risultati dipendono dalla piattaforma di misurazione?

### Perché [to compute] (da calcolare)?
Il supplemento mostra la struttura ma i calcoli finali non sono stati ancora eseguiti (o non sono riportati qui). Questo è comune nelle bozze di manoscritti.

---

## 📄 **PAGINA 7: eTable 5 - Studies Excluded at Full-Text Review**

### Categorie di esclusione:

#### **1. Meta-analisi/Revisioni (5 studi):**
- Puravet 2025 (Ann Emerg Med)
- Karamian 2025 (Eur J Trauma Emerg Surg)
- Amoo 2022 (Neurosurg Rev)
- Mondello 2021 (J Neurotrauma)
- Unden 2010 (J Head Trauma Rehabil)

**Perché esclusi:** Non sono studi originali, non hanno dati primari.

#### **2. Cohort Sovrapposte (4 studi):**
- **Oris 2024 (IJMS)** e **Oris 2024 (CCLM)**: Sottoinsieme di Puravet 2026 (stessi pazienti di Clermont-Ferrand)
- **Lapić 2024**: N=62, sovrapposizione con Milevoj 2025 (coorte Zagreb)
- **Welch 2025**: Stessi campioni ALERT-TBI già usati in Bazarian 2021

**Perché esclusi:** Evitare doppio conteggio degli stessi pazienti.

#### **3. N < 100 (3 studi):**
- Jalali 2025: N=123 (solo 7 CT+), Se=83,3% inaffidabile
- Curran 2025: N=89
- Biberthaler 2004: N=75

**Perché esclusi:** Troppo pochi pazienti, stime imprecise.

#### **4. Protocolli (1 studio):**
- Richard 2021: Solo protocollo BRAINI, risultati pubblicati poi in Lagares 2024

#### **5. Outcome Diverso (3 studi):**
- Blais Lécuyer 2021: Outcome = ICH clinicamente significativo (non lesione alla TC)
- Faisal 2023: Outcome = adesione a linee guida
- Boucher 2023: Outcome = PPCS (sintomi post-concussionali) a 90 giorni

#### **6. No 2×2 Data (2 studi):**
- Papa 2023: Solo AUC (area sotto la curva), senza dati 2×2 al cutoff predefinito
- Calcagnile 2012: Studio interventistico, dati non estraibili

#### **7. Altro (1 studio):**
- de Kruijk 2001: Nessun gold standard TC, solo confronto con controlli sani

---

## 📄 **PAGINA 8: PRISMA-DTA Checklist**

### Cos'è PRISMA-DTA?
**PR**eferred **R**eporting **I**tems for **S**ystematic reviews and **M**eta-**A**nalyses of **D**iagnostic **T**est **A**ccuracy

È una checklist di 27 item che devono essere riportati in ogni revisione sistematica di test diagnostici.

### Esempi di item:

| Item | Contenuto | Dove trovarlo |
|------|-----------|---------------|
| 1 | Titolo identifica come revisione DTA | Titolo del manoscritto |
| 5 | Registrazione protocollo | Methods - Protocol |
| 8 | Strategia di ricerca completa | eTable 1 (questo supplemento!) |
| 12 | Valutazione rischio di bias | Methods + eTable 3 |
| 15 | Valutazione rischio di bias across studies | Results - Publication Bias |

**Nota:** "To be completed" significa che la checklist deve essere compilata con i numeri di pagina prima della sottomissione finale.

---

## 🧮 **Riepilogo Statistico del Supplemento**

### Metodi chiave utilizzati:
1. **REML** su scala logit per pooling
2. **I², τ², Q** per eterogeneità
3. **Modello bivariato Reitsma** per SROC
4. **Egger test** per publication bias
5. **Cross-validation** R vs Python

### Dati grezzi disponibili:
- 16 studi con tabelle 2×2 complete
- 10 nell'analisi primaria, 6 in analisi di sensitività

### Qualità degli studi:
- 6 studi con Risk of Basso su tutti i domini
- Alcuni studi con limitazioni (popolazioni selezionate, dati da abstract)

---

## 💡 **Perché il Supplemento è Importante?**

### Per i ricercatori:
- Permette di **replicare** l'analisi
- Mostra tutti i dettagli metodologici
- Fornisce i dati grezzi per meta-analisi future

### Per i revisori:
- Verifica che i metodi siano appropriati
- Controlla la qualità degli studi inclusi
- Valuta la trasparenza della ricerca

### Per i lettori interessati:
- Capisce perché certi studi sono stati esclusi
- Vede la qualità individuale di ogni studio
- Comprende le limitazioni dell'analisi

---

*Fine della guida al supplemento*
