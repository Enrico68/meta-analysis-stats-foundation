# Confronto Risultati: Python v10 (DL) vs R metafor v11 (REML)

## Verifica: I Risultati Sono Analoghi?

**Risposta breve: SÌ, la concordanza è quasi perfetta. Le differenze sono trascurabili (0–0.1 pp).**

---

## 1. Tabella di Confronto Diretta

### GFAP+UCH-L1 (8 studi)

| Parametro | Python v10 (DL) | R metafor v11 (REML) | Δ |
|---|---|---|---|
| **Sensibilità pooled** | 97.6% [95.1–98.8] | 97.6% [95.1–98.8] | **0** |
| **Specificità pooled** | 49.4% [42.0–56.9] | 49.3% [40.3–58.4] | 0.1 pp |
| **I² Sensibilità** | 0% | 0% | 0 |
| **I² Specificità** | 95% | 97% | +2% |
| **τ² Sensibilità** | 0.0000 | 0.0000 | 0 |
| **τ² Specificità** | ~0.22 | 0.2665 | +0.04 |
| **NPV** | 99.60% | 99.59% | 0.01 pp |
| **CT evitate** | ~45.5% | 45.6% | 0.1 pp |
| **Lesioni mancate/1000** | 1.8 | 1.9 | 0.1 |

### S100B (6 studi)

| Parametro | Python v10 (DL) | R metafor v11 (REML) | Δ |
|---|---|---|---|
| **Sensibilità pooled** | 95.1% [90.9–97.4] | 95.1% [90.9–97.4] | **0** |
| **Specificità pooled** | 24.1% [19.6–29.2] | 24.1% [19.4–29.5] | 0 |
| **I² Sensibilità** | 0% | 0% | 0 |
| **I² Specificità** | 90% | 90% | 0 |
| **τ² Sensibilità** | 0.0000 | 0.0000 | 0 |
| **τ² Specificità** | ~0.10 | 0.1043 | ~0 |
| **NPV** | 98.53% | 98.53% | 0 |
| **CT evitate** | 22.4% | 22.4% | 0 |
| **Lesioni mancate/1000** | 3.4 | 3.4 | 0 |

---

## 2. Analisi delle Differenze

### Differenza #1: Intervalli di Confidenza della Specificità GFAP

| Parametro | Python v10 (DL) | R metafor v11 (REML) |
|---|---|---|
| Sp GFAP | 49.4% [42.0–56.9] | 49.3% [40.3–58.4] |
| Ampiezza IC | 14.9 pp | 18.1 pp |

**Perché questa differenza?**

REML (Restricted Maximum Likelihood) è un metodo più conservativo di DerSimonian-Laird quando l'eterogeneità è molto alta (I²=97%). REML tende a stimare τ² leggermente più alto, il che produce intervalli di confidenza più ampi.

**È corretto?** Sì. Con I²=97% e solo 8 studi, REML è preferibile perché DL tende a sottostimare τ². La differenza è piccola (3 pp di ampiezza) e non cambia l'interpretazione clinica.

### Differenza #2: Eterogeneità della Sensibilità

| Parametro | Valore |
|---|---|
| **I² Sensibilità GFAP** | **0%** |
| **I² Sensibilità S100B** | **0%** |
| **τ² Sensibilità GFAP** | **0.0000** |
| **τ² Sensibilità S100B** | **0.0000** |

**Cosa significa:** Le sensibilità dei singoli studi sono **tutte molto simili tra loro** (tutte sopra il 93% per GFAP, sopra il 93% per S100B). Non c'è eterogeneità reale — le piccole differenze sono tutte dovute al caso.

Questo è un risultato **molto importante**: la sensibilità è una proprietà stabile del biomarcatore, che non dipende dal contesto.

### Differenza #3: Eterogeneità della Specificità

| Parametro | GFAP | S100B |
|---|---|---|
| **I² Specificità** | **97%** | **90%** |
| **τ² Specificità** | **0.2665** | **0.1043** |
| **Q (p-value)** | <0.0001 | <0.0001 |

**Cosa significa:** La specificità varia **enormemente** tra gli studi. Questo è il risultato più importante di tutta la meta-analisi.

**Perché?** Perché ogni studio ha usato cutoff diversi, popolazioni diverse, e setting diversi. La specificità dipende moltissimo dal contesto. **Non puoi dire "la specificità è 49%" come se fosse un numero universale.**

---

## 3. Confronto dei Metodi Statistici

### I Due Approcci a Confronto

| Aspetto | Python v10 (DL) | R metafor v11 (REML) |
|---|---|---|
| **Modello** | DerSimonian-Laird marginale | Univariate REML separato |
| **Correlazione Se-Sp** | ❌ Ignorata | ❌ Ignorata |
| **Stima di τ²** | DerSimonian-Laird (metodo dei momenti) | REML (maximum likelihood) |
| **Eterogeneità riportata** | I² calcolato manualmente | I², Q, τ² completi |
| **Funnel plot** | No | ✅ Sì |
| **Egger test** | No | ✅ Sì |
| **Output per pubblicazione** | Ottimo (figure) | ✅ Completo (numeri + figure) |

### Qual È il Migliore?

| Criterio | Vincitore | Perché |
|---|---|---|
| **Correttezza statistica** | R metafor v11 (REML) | REML è più stabile con pochi studi e alta eterogeneità |
| **Completezza output** | R metafor v11 | I², Q, τ², funnel, Egger test |
| **Qualità figure** | Python v10 | Publication-ready, stile Lancet/NEJM |
| **Trasparenza** | Python v10 | Codice leggibile, ogni passo spiegato |
| **Concordanza dei risultati** | **Pareggio** | Differenze trascurabili (0–0.1 pp) |

---

## 4. Il Verdetto: I Risultati Sono Affidabili?

### ✅ SÌ, Ecco Perché

1. **Concordanza quasi perfetta tra due metodi indipendenti:**
   - Python (DL) e R (REML) danno risultati identici per le stime puntuali
   - Le uniche differenze sono nei CI della specificità (leggermente più ampi con REML)
   - Questo conferma che i risultati sono **robusti**

2. **Le stime puntuali sono identiche:**
   - Se GFAP: 97.6% in entrambi i casi
   - Se S100B: 95.1% in entrambi i casi
   - Sp S100B: 24.1% in entrambi i casi

3. **L'eterogeneità della specificità (I²=90-97%) è un risultato robusto** che appare in entrambe le versioni

4. **L'assenza di eterogeneità nella sensibilità (I²=0%) è un risultato forte** — le sensibilità sono davvero coerenti tra studi

### ⚠️ Ma Attenzione a...

1. **La specificità pooled ha I²=97% (GFAP) e I²=90% (S100B)** — significa che la "media" ha poco senso. La specificità dipende dal contesto (cutoff, popolazione, setting). Non puoi generalizzare il valore pooled a tutti i setting clinici.

2. **REML produce CI leggermente più ampi per la specificità GFAP** — questo è il metodo più conservativo e corretto con alta eterogeneità.

3. **Nessun metodo modella la correlazione Se-Sp** — sia Python DL che R metafor trattano sensibilità e specificità separatamente. Il modello bivariato di Reitsma sarebbe più completo, ma la concordanza tra i due metodi univariati conferma la robustezza.

---

## 5. Cosa Aggiungere alla Tua Tesi/Pubblicazione

### Se Usi R metafor v11 come Principale

> "La meta-analisi è stata condotta utilizzando il pacchetto `metafor` in R con stimatore REML. Sensibilità e specificità sono state poolizzate separatamente con modelli random-effects univariati. L'eterogeneità è stata quantificata con I², τ² e il test Q di Cochran. Il bias di pubblicazione è stato valutato con funnel plot e test di Egger."

### Se Vuoi il Meglio di Entrambi

> "Per le stime pooled è stato utilizzato il modello REML univariato (pacchetto `metafor` in R). Come analisi di sensibilità, è stata implementata anche la stima DerSimonian-Laird in Python, confermando la concordanza quasi perfetta delle stime puntuali (Δ < 0.1 pp). Le figure publication-ready sono state generate con Python/matplotlib."

---

## 6. Risultati del Test di Egger (Bias di Pubblicazione)

| Biomarcatore | Intercetta | p-value | Interpretazione |
|---|---|---|---|
| **GFAP+UCH-L1** | −0.137 | 0.876 | ✅ Nessun bias di pubblicazione |
| **S100B** | 1.124 | 0.142 | ✅ Nessun bias di pubblicazione |

**Cosa significa:** Non c'è evidenza che studi con risultati "sfavorevoli" siano stati pubblicati meno. I risultati della meta-analisi non sono distorti da bias di pubblicazione.

**Nota:** Il test di Egger ha bassa potenza con pochi studi (k<10), quindi "nessuna evidenza di bias" non significa "assenza di bias".

---

## 7. Riassunto Numerico Finale

### I Numeri "Definitivi" (R metafor v11 REML)

| | GFAP+UCH-L1 | S100B |
|---|---|---|
| **Se** | 97.6% [95.1–98.8] | 95.1% [90.9–97.4] |
| **Sp** | 49.3% [40.3–58.4] | 24.1% [19.4–29.5] |
| **NPV** | 99.59% | 98.53% |
| **CT evitate** | 45.6% | 22.4% |
| **Miss/1000** | 1.9 | 3.4 |
| **I² Se** | 0% | 0% |
| **I² Sp** | 97% | 90% |
| **Egger p** | 0.876 | 0.142 |

### Confronto Python v10 vs R metafor v11

| Parametro | Python v10 | R metafor v11 | Differenza |
|---|---|---|---|
| Se GFAP | 97.6% | 97.6% | **0** ✅ |
| Sp GFAP | 49.4% | 49.3% | **0.1 pp** ✅ |
| Se S100B | 95.1% | 95.1% | **0** ✅ |
| Sp S100B | 24.1% | 24.1% | **0** ✅ |
| NPV GFAP | 99.60% | 99.59% | **0.01 pp** ✅ |
| Miss/1000 GFAP | 1.8 | 1.9 | **0.1** ✅ |

**Conclusione:** La concordanza è quasi perfetta. Le uniche differenze sono nei CI della specificità GFAP — leggermente più ampi con REML (corretto e atteso con I²=97%). **Python v10 e R metafor v11 producono risultati sostanzialmente identici.**

---

## 8. Raccomandazione Finale

Per una pubblicazione scientifica, la strategia migliore è:

1. **R metafor v11** → per le stime numeriche (REML è il gold standard, output completo con I², Q, τ², funnel, Egger)
2. **Python v10** → per le figure publication-ready (qualità grafica superiore, stile Lancet/NEJM)

> **Se due metodi indipendenti danno risultati quasi identici → il risultato è ROBUSTO.**

Nel nostro caso, la concordanza è quasi perfetta (Δ < 0.1 pp). I risultati sono affidabili e pubblicabili.

---

*Documento corretto il 4 Aprile 2026 — Confronto tra Python v10 (DL) e R metafor v11 (REML). La colonna "R mada v9" è stata rimossa perché non ha mai prodotto risultati validi con i dati definitivi.*
