# `getP()` in `shipVertex.py` — Gleichungen mit Begründung

Funktion: `shipVertex.py`, verschachtelt in `TwoTrackVertex()`.
Aufruf: `P, covP = getP(values, emat, m1, m2)`

Zweck: aus dem **9‑Parameter‑Vertexfit** (TMinuit) den **4‑Impuls des
Mutterteilchens** (HNL‑Kandidat) und dessen **4×4‑Kovarianzmatrix** berechnen.

---

## 0. Eingaben

`values` = die 9 Fit‑Parameter aus `gMinuit.GetParameter(...)`:

| Index | Name (Minuit) | Bedeutung | Startwert im Code |
|------|----------------|-----------|-------------------|
| 0,1,2 | `X/Y/Z pos` | Vertexposition (in `getP` **nicht** benutzt) | `HNLPos` |
| 3 | `tan1X` = $a_3$ | $p_{x1}/p_{z1}$ (Steigung Spur 1 in x) | `mom1[0]/mom1[2]` |
| 4 | `tan1Y` = $a_4$ | $p_{y1}/p_{z1}$ | `mom1[1]/mom1[2]` |
| 5 | `1/mom1` = $a_5$ | $1/\lvert\vec p_1\rvert$ | `1/mom1.Mag()` |
| 6 | `tan2X` = $a_6$ | $p_{x2}/p_{z2}$ | `mom2[0]/mom2[2]` |
| 7 | `tan2Y` = $a_7$ | $p_{y2}/p_{z2}$ | `mom2[1]/mom2[2]` |
| 8 | `1/mom2` = $a_8$ | $1/\lvert\vec p_2\rvert$ | `1/mom2.Mag()` |

`emat` = 9×9 Kovarianzmatrix dieser Parameter (`gMinuit.mnemat`), flach als 81 Zahlen.

`m1, m2` = angenommene Massen der beiden Spuren, aus
`self.PDG.GetParticle(pdgCode).Mass()`. `pdgCode` ist die **Fit‑Hypothese**
der Spur (praktisch immer $\mu^\pm$, PDG $\pm13$; Protonen werden vorher zu
Pionen umgesetzt). Die Masse kommt also **aus der Reco‑Hypothese, nicht aus
Truth**.

Abkürzungen im Code:

$$
A_5 \equiv 1 + a_3^2 + a_4^2,\qquad
A_8 \equiv 1 + a_6^2 + a_7^2 .
$$

---

## 1. Impuls jeder Spur

```python
px1 = a3 / (a5 * sqrt(1 + a3**2 + a4**2))
py1 = a4 / (a5 * sqrt(1 + a3**2 + a4**2))
pz1 = 1  / (a5 * sqrt(1 + a3**2 + a4**2))
```

**Begründung.** Gegeben sind die Steigungen $t_x = a_3 = p_x/p_z$,
$t_y = a_4 = p_y/p_z$ und der Betrag $\lvert\vec p_1\rvert = 1/a_5$.
Aus

$$
\lvert\vec p_1\rvert^2 = p_x^2+p_y^2+p_z^2 = p_z^2\,(t_x^2+t_y^2+1) = p_z^2\,A_5
$$

folgt (mit der Annahme $p_z>0$, für SHiP‑Spuren immer erfüllt)

$$
p_{z1} = \frac{\lvert\vec p_1\rvert}{\sqrt{A_5}} = \frac{1}{a_5\sqrt{A_5}},\qquad
p_{x1} = t_x\,p_{z1} = \frac{a_3}{a_5\sqrt{A_5}},\qquad
p_{y1} = \frac{a_4}{a_5\sqrt{A_5}} .
$$

Analog für Spur 2 mit $a_6,a_7,a_8,A_8$. Diese Umrechnung ist **exakt**.

---

## 2. Gesamtimpuls des Mutterteilchens

```python
Px = px1 + px2
Py = py1 + py2
Pz = pz1 + pz2
```

**Begründung.** Impulserhaltung am Zerfallsvertex: $\vec P = \vec p_1 + \vec p_2$.

---

## 3. Energien der beiden Spuren

```python
E1 = sqrt(px1**2 + py1**2 + pz1**2 + m1**2)
E2 = sqrt(px2**2 + py2**2 + pz2**2 + m2**2)
```

**Begründung.** Relativistische Energie‑Impuls‑Beziehung
$E_i = \sqrt{\lvert\vec p_i\rvert^2 + m_i^2}$.
Wegen $\lvert\vec p_i\rvert = 1/a_5$ bzw. $1/a_8$ gilt äquivalent
$E_1 = \sqrt{1/a_5^2 + m_1^2}$, $E_2 = \sqrt{1/a_8^2 + m_2^2}$ —
$E_i$ hängt **nur** vom jeweiligen $a_5$ / $a_8$ ab, nicht von den Steigungen.

---

## 4. Invariante Masse

```python
M = sqrt(2*E1*E2 + m1**2 + m2**2 - 2*pz1*pz2*(1 + a3*a6 + a4*a7))
```

**Begründung.** Definition über das Quadrat des Gesamt‑4‑Impulses:

$$
M^2 = (E_1+E_2)^2 - \lvert\vec p_1+\vec p_2\rvert^2 .
$$

Ausmultiplizieren und $E_i^2 - \lvert\vec p_i\rvert^2 = m_i^2$ einsetzen:

$$
M^2 = m_1^2 + m_2^2 + 2E_1E_2 - 2\,\vec p_1\!\cdot\!\vec p_2 .
$$

Das Skalarprodukt schreibt sich mit $p_{x i} = a\,p_{z i}$ als

$$
\vec p_1\!\cdot\!\vec p_2
= p_{x1}p_{x2} + p_{y1}p_{y2} + p_{z1}p_{z2}
= p_{z1}p_{z2}\,\big(a_3 a_6 + a_4 a_7 + 1\big) .
$$

Einsetzen ergibt genau die Codezeile. Diese Form ist **exakt** (nur umgestellt,
damit man `pz1`, `pz2` direkt wiederverwenden kann).

```python
P = ROOT.TLorentzVector()
P.SetXYZM(Px, Py, Pz, M)          # speichert (Px, Py, Pz, E) mit E = sqrt(P^2 + M^2)
```

---

## 5. Jacobi‑Matrix `M_AtoP` (4×6)

Fehlerfortpflanzung braucht

$$
J \equiv \frac{\partial (P_x, P_y, P_z, M)}{\partial (a_3, a_4, a_5, a_6, a_7, a_8)} .
$$

Zeilen 0–2 = $\partial \vec P/\partial a_k$, Zeile 3 = $\partial M/\partial a_k$.
Im Code `MM = 2*M`, weil $\partial M/\partial a_k = \dfrac{1}{2M}\,\partial M^2/\partial a_k$.

### Zeilen 0–2 (Impulsableitungen) — exakt

Beispiel $\partial P_x/\partial a_3$ (nur $p_{x1}$ hängt von $a_3$ ab):

$$
\frac{\partial}{\partial a_3}\!\left(\frac{a_3}{a_5\sqrt{A_5}}\right)
= \frac{1}{a_5}\Big(A_5^{-1/2} - a_3^2 A_5^{-3/2}\Big)
= \frac{1 - a_3^2/A_5}{a_5\sqrt{A_5}} .
$$

→ `M_AtoP[0][0] = (1 - a3*a3/A5)/(a5*sqrt(A5))`. Ebenso:

$$
\frac{\partial P_x}{\partial a_4} = \frac{-a_3 a_4/A_5}{a_5\sqrt{A_5}},\qquad
\frac{\partial P_x}{\partial a_5} = \frac{-a_3}{a_5^2\sqrt{A_5}},
$$

$$
\frac{\partial P_z}{\partial a_3} = \frac{-a_3/A_5}{a_5\sqrt{A_5}},\qquad
\frac{\partial P_z}{\partial a_5} = \frac{-1}{a_5^2\sqrt{A_5}} .
$$

Die Spalten 3–5 (Spur 2) sind die exakte Spiegelung mit $a_3\!\to\!a_6$,
$a_4\!\to\!a_7$, $a_5\!\to\!a_8$, $A_5\!\to\!A_8$. Zeile 1 ($P_y$) ist Zeile 0
mit $a_3\leftrightarrow a_4$.

### Zeile 3 (Massenableitungen)

Mit der Abkürzung

```python
a5a8 = a5 * a8 * sqrt(A5) * sqrt(A8)        #  = 1/(pz1*pz2)
```

**Ableitung nach $a_5$ (und $a_8$) — exakt.**
Aus $M^2 = m_1^2+m_2^2+2E_1E_2 - 2p_{z1}p_{z2}(1+a_3a_6+a_4a_7)$ mit
$\partial E_1/\partial a_5 = -1/(a_5^3 E_1)$ und
$\partial p_{z1}/\partial a_5 = -1/(a_5^2\sqrt{A_5})$:

$$
\frac{\partial M^2}{\partial a_5}
= \frac{2(1+a_3a_6+a_4a_7)}{a_5\,(a_5a_8)} - \frac{2E_2}{a_5^3 E_1},
\qquad
\frac{\partial M}{\partial a_5} = \frac{1}{MM}\,\frac{\partial M^2}{\partial a_5}.
$$

→ `M_AtoP[3][2]` (und spiegelbildlich `M_AtoP[3][5]`) stimmen **exakt**.

**Ableitung nach den Steigungen $a_3, a_4, a_6, a_7$ — genäherte Form.**
Der Code setzt z. B.

```python
M_AtoP[3][0] = (-2*a6/a5a8 + 2*a3*E2/(a5*a5*A5*E1)) / MM
```

Der **erste Term** $-2a_6/(a_5a_8)$ ist exakt: er stammt aus
$-2\,p_{z1}p_{z2}\,\partial(a_3a_6)/\partial a_3 = -2 p_{z1}p_{z2} a_6$ und
$p_{z1}p_{z2} = 1/(a_5a_8)$.

Der **zweite Term** ist eine Näherung. Exakt (da $E_1$ **nicht** von $a_3$
abhängt) wäre nur der Beitrag aus $\partial p_{z1}/\partial a_3$:

$$
\left.\frac{\partial M}{\partial a_3}\right|_{\text{exakt}}
= \frac{1}{MM}\left[-\frac{2a_6}{a_5a_8}
+ \frac{2a_3\,p_{z2}(1+a_3a_6+a_4a_7)}{a_5\,A_5^{3/2}}\right].
$$

Der Code ersetzt den zweiten Summanden durch $2a_3 E_2/(a_5^2 A_5 E_1)$.
Beide sind identisch **im masselosen Limes** $E_i\to\lvert\vec p_i\rvert$;
der relative Fehler ist $\mathcal O(m^2/p^2)$ (für $\mu$ bei GeV‑Impulsen
$\lesssim$ ein Promille). Betrifft **nur die Massenunsicherheit** `covP`,
nicht $M$ selbst. Das `# fixme: mass from track reconstruction needed` im
Code bezieht sich auf diese Baustelle.

`M_AtoP[3][3]`, `M_AtoP[3][4]` sind die Spiegelung ($1\leftrightarrow2$).

### Transponierte

```python
for i in range(4):
    for j in range(6):
        MT_AtoP[j][i] = M_AtoP[i][j]      # MT_AtoP = J^T  (6x4)
```

---

## 6. Parameter‑Kovarianz `covA` (6×6)

```python
for i in range(36):
    covA[i//6][i%6] = cov[i//6 + 3 + (i%6 + 3)*9]
```

D. h. $\text{covA}[r][c] = \text{emat}\big[(r+3) + 9\,(c+3)\big]$ für
$r,c \in \{0,\dots,5\}$.

**Begründung.** `emat` ist die volle 9×9‑Fit‑Kovarianz. `getP` braucht nur
den Block der **6 Impulsparameter** $a_3\dots a_8$ (Indizes 3–8); die 3
Positionsparameter (0–2) werden weggelassen, weil $\vec P$ und $M$ nicht von
der Vertexlage abhängen. Der Offset `+3` bzw. `+3` in Zeile/Spalte schneidet
genau diesen Unterblock heraus. (`emat` ist symmetrisch, daher ist die
Reihenfolge Zeile/Spalte im Flachindex egal.)

---

## 7. Kovarianz des 4‑Impulses `covP` (4×4)

```python
tmp  = M_AtoP · covA          # (4x6)
covP = tmp   · MT_AtoP        # (4x4)
```

$$
\mathrm{cov}(P_x,P_y,P_z,M) = J\,\Sigma_{a}\,J^{\mathsf T},
\qquad \Sigma_a = \text{covA}.
$$

**Begründung.** Lineare (Gauß’sche) Fehlerfortpflanzung: für
$y = f(x)$ mit Kovarianz $\Sigma_x$ gilt in erster Ordnung
$\Sigma_y = J\Sigma_x J^{\mathsf T}$ mit $J = \partial f/\partial x$.
Gültig, solange die Fit‑Fehler klein gegen die Nichtlinearität von $f$ sind.

Rückgabe: `P` (`TLorentzVector`) und `covP` (4×4, Reihenfolge $P_x,P_y,P_z,M$).

Im Aufrufer wird `covP` auf die 10 unabhängigen Elemente der oberen
Dreiecksmatrix reduziert und mit `particle.SetCovP(covP)` gespeichert:

```
[00, 01, 02, 03, 11, 12, 13, 22, 23, 33]
```

---

## Zusammenfassung: exakt vs. genähert

| Größe | Status |
|-------|--------|
| $p_{x,y,z}$ pro Spur | exakt (gegebene Parametrisierung) |
| $\vec P = \vec p_1+\vec p_2$ | exakt |
| $E_1, E_2$, $M$ | exakt |
| `M_AtoP` Zeilen 0–2 ($\partial\vec P$) | exakt |
| `M_AtoP[3]` Spalten $a_5,a_8$ | exakt |
| `M_AtoP[3]` Spalten $a_3,a_4,a_6,a_7$ | genähert ($\mathcal O(m^2/p^2)$), nur `covP` betroffen |
| `covP = J covA J^T` | lineare Fehlerfortpflanzung (1. Ordnung) |

**Wichtigster physikalischer Punkt:** $m_1, m_2$ sind die **Reco‑Massen­hypothesen**
der Spuren (i. d. R. Myon für beide), nicht die wahren Teilchenmassen.
`selected_mom.M()` in der Analyse ist genau dieses $M$.
