Modelo BIPV/T no CitySim Pro – Especificação Técnica

Este documento especifica **como o CitySim Pro deve ser alterado** para contabilizar o balanço energético de uma fachada com PV ventilado (**BIPV/T**) quando indicado no ficheiro XML de entrada. fileciteturn1file0

---

## 1. Condição de ativação do modelo BIPV/T

Cada `Wall` no XML passa a poder ter um parâmetro:

```xml
<Wall ... WALLPVTYPE="1" ... />
```

- Se **`WALLPVTYPE="1"`**, essa superfície deve:
  - ser tratada como **fachada BIPV/T ventilada**;
  - usar o sistema de 13 equações descrito neste ficheiro;
  - substituir o balanço energético “normal” dessa superfície.

- Se **`WALLPVTYPE` não existir** ou **`WALLPVTYPE!="1"`**:
  - o programa mantém o **modelo original de balanço de parede** (sem BIPV/T).

Além disso, no XML deve existir:

```xml
<WALLPVTYPE id="1" name="BIPV/T" />
<WALLPV name="1"
        kbp="..."
        ebp="..."
        alfacsw="..."
        taugsw="..."
        alfacir="..."
        taugir="..."
        nuetamp="..."
        keva="..."
        eeva="..."
        kairg="..."
        eairg="..."
        etaeleref="..."
        Tcellref="..."
        algagsw="..."
        algagir="..."
        eair="..."
        epsilonbp="..."
        epsilonos="..."
        epsilongo="..."
        epsilongi="..."
/>
```

Cada atributo mapeia para um parâmetro físico do modelo (ver Secção 4.3).

---

## 2. Vetor de incógnitas do sistema (13 variáveis)

Para cada parede com BIPV/T ativa, em cada passo de tempo, o CitySim deve resolver o vetor:

\[
X =
\begin{bmatrix}
T_{os} &
T_{bp} &
T_{cell} &
T_{gi} &
T_{go} &
h_{os} &
h_{bp} &
T_{air,mean} &
T_{air,out} &
X &
Y &
V_{air} &
H
\end{bmatrix}^T
\]

onde:

1. **T_os** – temperatura da superfície exterior do módulo PV  
2. **T_bp** – temperatura da **backplate** PV  
3. **T_cell** – temperatura das células fotovoltaicas  
4. **T_gi** – temperatura do vidro interior  
5. **T_go** – temperatura do vidro exterior  
6. **h_os** – coeficiente de convecção na superfície exterior PV  
7. **h_bp** – coeficiente de convecção na superfície da backplate  
8. **T_air,mean** – temperatura média do ar na cavidade ventilada  
9. **T_air,out** – temperatura do ar à saída da cavidade (na altura L)  
10. **X** – parâmetro de escoamento (função de integração da equação de energia do ar)  
11. **Y** – parâmetro de escoamento (coeficiente de decaimento)  
12. **V_air** – velocidade média do ar na cavidade  
13. **H** – termo de aquecimento (heating) associado ao balanço de energia do espaço interior

> **IMPORTANTE:** Todas as temperaturas usadas nas equações não lineares (termos de ordem superior, por exemplo T⁴) devem estar em **Kelvin**. fileciteturn1file0

---

## 3. Forma matricial do sistema

O sistema deve ser implementado na forma:

\[
\mathbf{0} = A \cdot X + B + C
\]

onde:

- **A** é uma matriz 13×13 com coeficientes lineares;
- **B** é um vetor de termos conhecidos (fontes, condições de fronteira, termos lineares independentes);
- **C** é um vetor de termos **não lineares** (convecção, radiação, funções de escoamento, etc.).

Na prática, para cada equação i:

\[
F_i(X) = \sum_j A_{ij} X_j + B_i + C_i(X) = 0
\]

---

## 4. Coeficientes do sistema A, B, C

### 4.1 Estrutura da matriz A (entradas não nulas)

A ordem das equações \(F_1 \dots F_{13}\) segue o capítulo original. fileciteturn1file0

#### Equação 1 – Balanço na superfície exterior PV (T_os)

- Incógnita principal: T_os

\[
A_{1,1} = \frac{K_{w1}}{2} \frac{\Delta t}{C_w}
+ \Delta t (k_2 + K_{w1})
- K_{w1}
- K_e (1 - pvratio)
\]

B₁:

\[
B_1 = K_e (1 - pvratio) T_{ext}
+ (Q_{sun,1} + Q_{ir})(1 - pvratio)
+ \frac{K_{w1}}{C_w} T_w^{n-1}
+ \Delta t \, k_2 T_a
+ \Delta t \frac{k_2}{K_i} (Q_{sun,2} W_w + L_r)
\]

C₁:

\[
C_1 = h_{os} \, pvratio \, (T_{air,mean} - T_{os})
+ pvratio \, \sigma \, (\varepsilon_{bp} T_{bp}^4 - \varepsilon_{os} T_{os}^4)
\]

---

#### Equação 2 – Balanço na backplate (T_bp)

\[
A_{2,2} = -\frac{k_{bp}}{e_{bp}}
\quad,\quad
A_{2,3} = \frac{k_{bp}}{e_{bp}}
\]

(sem termos B₂; C₂ inclui convecção e radiação – ver equação 3)

---

#### Equação 3 – Balanço nas células PV (T_cell)

\[
A_{3,2} = \frac{k_{bp}}{e_{bp}} A_{wall}
\]

\[
A_{3,3} =
- (\alpha_{c,sw} \tau_{g,sw}^2 Q_{sun,1}
+ \alpha_{c,ir} \tau_{g,ir}^2 Q_{ir}) \, \mu \eta_{mp}
- \frac{k_{bp}}{e_{bp}} A_{wall}
- \frac{k_{eva}}{e_{eva}} A_{wall}
\]

\[
A_{3,4} = \frac{k_{eva}}{e_{eva}} A_{wall}
\]

B₃:

\[
B_3 =
(\alpha_{c,sw} \tau_{g,sw}^2 Q_{sun,1}
+ \alpha_{c,ir} \tau_{g,ir}^2 Q_{ir})
-
(\alpha_{c,sw} \tau_{g,sw}^2 Q_{sun,1}
+ \alpha_{c,ir} \tau_{g,ir}^2 Q_{ir})
(\eta_{ele,ref} - \mu \eta_{mp} T_{cell,ref})
\]

C₂ (backplate):

\[
C_2 =
h_{bp} (T_{air,mean} - T_{bp})
+ \sigma (\varepsilon_{os} T_{os}^4 - \varepsilon_{bp} T_{bp}^4)
\]

---

#### Equação 4 – Vidro interior (T_gi)

\[
A_{4,3} = \frac{k_{eva}}{e_{eva}}
\quad,\quad
A_{4,4} = -\frac{k_{air,g}}{e_{air,g}} - \frac{k_{eva}}{e_{eva}}
\quad,\quad
A_{4,5} = \frac{k_{air,g}}{e_{air,g}}
\]

C₄:

\[
C_4 = \sigma(\varepsilon_{go} T_{go}^4 - \varepsilon_{gi} T_{gi}^4)
\]

---

#### Equação 5 – Vidro exterior (T_go)

\[
A_{5,5} = -K_e - \frac{k_{air,g}}{e_{air,g}} A_{wall}
\quad,\quad
A_{5,4} = \frac{k_{air,g}}{e_{air,g}} A_{wall}
\]

B₅:

\[
B_5 = K_e T_{ext} + (\alpha_{g,sw} Q_{sun,1} + \alpha_{g,ir} Q_{ir})
\]

C₅:

\[
C_5 = \sigma A_{wall} (\varepsilon_{gi} T_{gi}^4 - \varepsilon_{go} T_{go}^4)
\]

---

#### Equações 6 e 7 – h_os e h_bp (coef. de convecção)

Entradas da matriz:

\[
A_{6,12} =
\begin{cases}
0 & , T_a \le T_{ext} \\
-0.85 \times 1.33 & , T_a > T_{ext}
\end{cases}
\]

\[
A_{7,12} =
\begin{cases}
0 & , T_a \le T_{ext} \\
-0.85 \times 1.33 & , T_a > T_{ext}
\end{cases}
\]

B₆, B₇:

\[
B_6 = B_7 =
\begin{cases}
0 & , T_a \le T_{ext} \\
-0.85 \times 1.959 & , T_a > T_{ext}
\end{cases}
\]

C₆, C₇ (naturais vs mistas, ver Secção 5).

---

#### Equação 8 – Parâmetro X (escoamento)

\[
C_8 =
\begin{cases}
-(T_a - X)\dfrac{1 - \exp(-Y L)}{Y L}, & T_a \le T_{ext}\\[4pt]
-(T_{ext} - X)\dfrac{1 - \exp(-Y L)}{Y L}, & T_a > T_{ext}
\end{cases}
\]

---

#### Equação 9 – Parâmetro Y (escoamento)

\[
C_9 =
\begin{cases}
-(T_a - X) \exp(-Y L), & T_a \le T_{ext}\\[4pt]
-(T_{ext} - X) \exp(-Y L), & T_a > T_{ext}
\end{cases}
\]

---

#### Equação 10 – Temperatura média do ar

\[
C_{10} = -\dfrac{h_{os} T_{os} - h_{bp} T_{bp}}{h_{os} + h_{bp}}
\]

---

#### Equação 11 – Energia do ar na cavidade

\[
C_{11} = -\dfrac{h_{os} + h_{bp}}{e \, V_{air} \, \rho_{air} \, c_{p,air}}
\]

---

#### Equação 12 – Velocidade do ar V_air

\[
C_{12} =
\begin{cases}
- C_D e w 
\sqrt{\dfrac{2 g L |T_{air,mean} - T_a|}{T_a}}, & T_a \le T_{ext}\\[6pt]
- C_D e w 
\sqrt{\dfrac{2 g L |T_{air,mean} - T_{ext}|}{T_{ext}}}, & T_a > T_{ext}
\end{cases}
\]

\[
B_{12} =
\begin{cases}
0, & T_a \le T_{ext}\\[4pt]
- C_v e w v_{wind}, & T_a > T_{ext}
\end{cases}
\]

onde \(v_{wind}\) é a **componente de vento normal à fachada**. fileciteturn1file0

---

#### Equação 13 – Balanço de energia do espaço interior (H)

\[
A_{13,1} = -\dfrac{k_2 \Delta t K_{w1}}{\Delta t(-k_2 - K_{w1}) - C_w}
\]

\[
A_{13,12} =
\begin{cases}
- e w \rho_{air} c_{p,air} T_a, & T_a \le T_{heat} \land T_{air,out} > T_a\\[4pt]
0, & \text{else}
\end{cases}
\]

B₁₃:

\[
B_{13} =
\left(
-\dfrac{C_i}{\Delta t}
- U A_n - k_2 - \dfrac{k_2^2 \Delta t}{2}
\right)
\dfrac{T_a}{\Delta t(-k_2 - K_{w1}) - C_w}
+ \dfrac{C_i}{\Delta t} T_a^{n-1}
\]

\[
+ \big(U A_n T_{ext}
+ \dfrac{k_2}{K_{w1}}(Q_{sun,2} W_w + L_r)
+ Q_{sun,2} W_a + L_c\big)
- k_2
\dfrac{C_w T_w^{n-1} - \Delta t \frac{k_2}{K_i}(Q_{sun,2} W_w + L_c)}
{\Delta t(-k_2 - K_{w1}) - C_w}
\]

C₁₃:

\[
C_{13} =
\begin{cases}
e w \rho_{air} c_{p,air} T_{air,out} V_{air}, & T_a \le T_{heat} \land T_{air,out} > T_a\\[4pt]
0, & \text{else}
\end{cases}
\]

---

## 5. Propriedades do ar e números de Rayleigh

### 5.1 Propriedades termo-físicas do ar

Todas avaliadas à temperatura média da cavidade \(T_{air,mean}\):

\[
g = 9.81\ \text{m/s}^2
\]

\[
\rho_{air}(t) = \dfrac{101325}{287.05 \, T_{air,mean}(t)}
\]

\[
c_{p,air} = 1005 \ \text{J/(kg·K)}
\]

\[
\beta(t) = \dfrac{1}{T_{air,mean}(t)}
\]

\[
k_{air}(t) = 2.873\times 10^{-3} + 7.76\times 10^{-8} T_{air,mean}(t)
\]

\[
\nu_{air}(t) = 3.723\times 10^{-6} + 4.94\times 10^{-8} T_{air,mean}(t)
\]

\[
\alpha_{air}(t) = \dfrac{k_{air}(t)}{\rho_{air}(t) c_{p,air}}
\]

---

### 5.2 Números de Rayleigh

Exterior:

\[
Ra_{os}(t) =
\dfrac{g \beta(t) |T_{os}(t) - T_{air,mean}(t)| e^3}{\alpha_{air}(t) \nu_{air}(t)}
\]

Backplate:

\[
Ra_{bp}(t) =
\dfrac{g \beta(t) |T_{bp}(t) - T_{air,mean}(t)| e^3}{\alpha_{air}(t) \nu_{air}(t)}
\]

---

### 5.3 Coeficientes convectivos – natural vs mista

Convecção natural (Bar-Cohen & Rohsenow):

\[
\overline{Nu} =
\left[
\dfrac{144}{(Ra \cdot e/L)^2}
+ \dfrac{2.87}{(Ra \cdot e/L)^{1/2}}
\right]^{-1/2}
\]

\[
h = \dfrac{\overline{Nu} \, k_{air}}{e}
\]

Convecção mista (com vento), já incluída em C₆, C₇, B₆, B₇ e A₆,₁₂, A₇,₁₂ (coeficientes 0.85 × …).

---

## 6. Método numérico – Newton–Raphson

O sistema é resolvido por Newton–Raphson:

\[
X^{(i+1)} = X^{(i)} - J^{-1}(X^{(i)}) F(X^{(i)})
\]

onde o Jacobiano é calculado numericamente por diferenças finitas:

\[
\frac{\partial F_k}{\partial X_j}
\approx
\frac{F_k(X_1,\dots,X_j+\Delta X_j,\dots,X_n) - F_k(X_1,\dots,X_j,\dots,X_n)}
{\Delta X_j}
\]

---

## 7. Dados já fornecidos pelo CitySim Pro (não recalcular)

Segundo o resumo: fileciteturn1file0

Estes parâmetros já existem no código original e **devem ser reutilizados**:

- Kw1 – condutância da parede  
- Ke – coeficiente de convecção exterior (para superfície *Awall* apenas)  
- k2 – coeficiente de troca com o interior  
- pvratio – fração da parede coberta por PV (do XML)  
- Δt – deve ser fixado em 3600 s (1 h) para este modelo  
- Cw – capacidade térmica da parede  
- Awall – área da Wall surface  
- Qsun1(n), Qsun2(n) – radiações solares já calculadas  
- Qir(n) – radiação infravermelha  
- Text(n) – temperatura exterior (ficheiro climático)  
- Ta(n) – temperatura interior estimada pelo programa  
- Ci – capacidade térmica do ar interior  
- Theat – temperatura de set point de aquecimento  
- Tw(n-1) – temperatura da parede no passo anterior  
- Lc, Lr – ganhos internos convectivos e radiativos  
- UA(n) – coeficiente global perdas  
- Ta(n-1) – temperatura interior passo anterior  
- Ki – coeficiente ligado a trocas com o interior  
- vwind(n) – velocidade do vento (usar componente perpendicular à fachada!)

---

## 8. Parâmetros lidos do XML para o PV

Da entrada `<WALLPV>` (id = 1, BIPV/T), devem vir:

- kbp → `kbp`  
- ebp → `ebp`  
- alfacsw → αc,sw  
- taugsw → τg,sw  
- alfacir → αc,ir  
- taugir → τg,ir  
- nuetamp → μ ηmp  
- keva → k_eva  
- eeva → e_eva  
- kairg → k_air,g  
- eairg → e_air,g  
- etaeleref → ηele,ref  
- Tcellref → Tcell,ref  
- algagsw → αg,sw  
- algagir → αg,ir  
- eair → e (espessura da cavidade)  
- epsilonbp → εbp  
- epsilonos → εos  
- epsilongo → εgo  
- epsilongi → εgi  

---

## 9. Geometria da cavidade

- **L** – deve ser definido como a **altura mínima (em z)** da Wall surface.  
- **w** – definido como:

\[
w = pvratio \dfrac{A_{wall}}{L}
\]

- **Ww** e **Wa** continuam a ser usados como no programa original (não substituir por w).

---

## 10. Heating (H) e Cooling (C)

Este modelo define **H** (heating). Para **cooling (C)**, deve ser usada a **mesma lógica**.

A estrutura das equações mantém-se, mas o sinal de H ou C no balanço de energia global da zona é ajustado conforme a convenção já usada no CitySim Pro.

---

Este ficheiro deve servir de **base direta para implementação em C++**.
