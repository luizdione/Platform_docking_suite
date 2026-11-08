# Manual do modulo de drug design computacional
## Docking Suite - IFRJ

**Modulo:** `drug_design_module.py`  
**Aba do dashboard:** Drug Design  
**Scaffold parental:** o que voce informar. O default e a beta-lapachona, um produto  
natural publicado - nao ha composto soldado no modulo.  
**Referencia central do catalogo:** Pinto & De Castro, *Front. Physiol.* 2014;  
Pinto & De Castro, *Molecules* 2009

> **O que este manual documenta, e o que ele nao e.** O catalogo de R-groups, as
> rotas sinteticas e o modelo de dG vem de **SAR publicado** -- cada entrada carrega
> a citacao no campo `ref`. Nada aqui e resultado de pesquisa: o modulo propoe
> hipoteses de sintese a partir do scaffold que voce der a ele, e a prioridade e
> sempre validar no docking real.

---

## O que este módulo faz

O módulo de Drug Design gera **novos candidatos antiparasitários** in silico a partir dos compostos do banco de dados, usando três estratégias complementares:

1. **FBDD** — substituição sistemática de grupos funcionais (R-groups) na scaffold quinônica
2. **Hibridização** — fusão de fragmentos BRICS de dois compostos ativos (estratégia de design do LQB-118)
3. **Algoritmo Genético** — evolução guiada por função de fitness combinando ADMET, afinidade estimada e acessibilidade sintética

Cada candidato gerado recebe:
- SMILES canônico validado pelo RDKit
- Perfil ADMET computacional (Lipinski, QED, TPSA, SAS)
- ΔG estimado de ligação à TcTR (modelo QSAR empírico)
- **Rota sintética sugerida** baseada em protocolos estabelecidos
- **Racional científico** embasado na literatura citada
- Parágrafo de resultados pronto para relatório (EN/PT)

---

## Base química — por que a scaffold quinônica

### O alvo: tripanotiona redutase

Tripanossomatídeos não usam glutationa como tiol antioxidante principal: usam
**tripanotiona**, mantida reduzida pela tripanotiona redutase (TR). A via não existe em
mamíferos, e a contraparte humana mais próxima é a glutationa redutase (hGR), de sítio
ativo mais estreito. É essa diferença de arquitetura que abre margem a um ligante
seletivo — e é também por isso que a plataforma doca a enzima do parasito **sempre** ao
lado da contraparte humana. Sem o contra-alvo, "liga bem" não é seletividade.

### Por que quinona

Duas propriedades, documentadas na literatura antiparasitária, fazem do núcleo
1,4-quinônico um ponto de partida recorrente:

- **Ciclo redox.** A quinona é reduzida a semiquinona, que reage com O₂ e devolve a
  quinona gerando superóxido. O resultado é estresse oxidativo sustentado — e ele é
  relativamente seletivo para o parasito, cuja defesa tiol-dependente é a própria TR.
- **Substrato subversivo.** Quinonas competem por elétrons na face da flavina de
  flavoenzimas, desviando o fluxo que deveria reduzir o tiol fisiológico.

Daí o catálogo abaixo: cada R-group é uma modificação em C2 do núcleo com SAR
antiparasitário **publicado**, e o campo `ref` diz de onde vem.

> Os dois mecanismos não são mutuamente exclusivos nem automáticos. Qual deles domina
> depende da molécula e do alvo, e é pergunta experimental — o módulo gera a hipótese,
> não a resposta.

---

## Abordagem 1 — FBDD (Fragment-Based Drug Design)

### Princípio

Substitui o grupo funcional na posição **C2 da 1,4-naftoquinona** (entre as duas carbonilas) por R-groups selecionados de um catálogo de 20 grupos baseado em SAR publicado.

```
         O                     O
         ‖                     ‖
    ┌────C────┐           ┌────C────┐
    │         │           │         │
  ──C         C─── [X]    C         C─── [R-group novo]
    │         │           │         │
    └────C────┘           └────C────┘
         ‖                     ‖
         O                     O
    
    Composto parental        Análogo FBDD
```

### Catálogo de R-groups por categoria

#### Halogênios — ciclo redox e hal-bonds

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–Cl** | Cl₂/AcOH ou NCS/DMF | ↑logP +0.7; hal-bond com backbone; NQ1-Cl ativo (Literatura) | Salomão et al. 2013 BMC Microbiol |
| **–Br** | Br₂/AcOH ou NBS/AcOEt | Maior polarizabilidade que Cl; NQ1-Br IC50 ~0.8 μM vs T. cruzi | Duarte/De Castro et al. 2020 |
| **–F** | SELECTFLUOR® ou F-TEDA-BF₄ | Bioisóstero de H; C–F resiste CYP450; ↑metabolic stability | Rani 2022 Chem Biol Drug Des |

#### Grupos Sulfa — interação com Cys do sítio TcTR

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–SOCH₃** | mCPBA (1 equiv., 0°C, DCM) | ↑solubilidade vs –SCH₃; sulfóxidos metabolicamente estáveis | SAR revisão Rani 2022 |
| **–S-tiazolil** | SNAr + 2-mercaptotiazol/K₂CO₃/DMF | π-stacking com FAD de TcTR; bioisóstero imidazol | SAR tiazol-NQ |
| **–SO₂NH₂** | Oxone®/MeOH:H₂O | TPSA +42 Å²; H-bond com Glu18 e Tyr110 TcTR | De Castro group reviews |

#### Grupos Nitrogenados — H-bonds e ativação por TcNTR

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–NH₂** | SnCl₂/EtOH/HCl ou NH₃/autoclave | H-bond para Glu18; série 2-amino-NQ bem estudada | Pinto & De Castro 2009 Molecules |
| **–NO₂** | HNO₃ conc./H₂SO₄ (0°C) | Dupla ação: ROS + ativação por TcNTR (como BZN) | Docampo & De Souza 1978; Salomão 2013 |
| **–triazol** | CuAAC click chemistry | IC50 0.8–3.1 μM, mais ativo que BZN; síntese fácil >90% | Bahia et al. 2016 MedChemComm |
| **–imidazol** | SNAr de 2-halo-NQ + imidazol/K₂CO₃/DMF | Base dos naftoimidazóis N1/N2/N3 (Literatura) IC50 0.5–2.0 μM | Pinto & De Castro 2009 Molecules |
| **–pirazinaamino** | Buchwald-Hartwig ou SNAr com Cl-pirazina | N básico → ↑solubilidade em pH ácido (endossoma Leishmania) | Medicinal chemistry de aminopyrazina |

#### Grupos Oxigenados — solubilidade e bioisósteros

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–OH** (lawsona-like) | SeO₂/dioxano ou hidrólise K₂CO₃/MeOH | Precursor LQB-118; quelação Fe-SOD; TPSA +17 Å² | Maia 2011 Invest New Drugs |
| **–OCH₃** | Me₂SO₄/K₂CO₃/acetona | Bloqueia oxidação; logP +0.5; análogo atovaquona | Lopes 1978 Ann Trop Med |
| **–OCF₃** | AgOCF₃/DCM ou Umemoto | ↑↑lipofilicidade; resistência CYP450; antiprotozoais modernos | CF₃O medicinal chemistry |

#### Heterocíclicos Fusionados — π-stacking com FAD

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–isoxazolil** | [3+2] cicloaddição, sem metal | Bioisóstero imidazol; antiparasitários modernos | MDPI Pharmaceutics 2022 |
| **–pirazolil** | Condensação hidrazina + β-cetoéster | Análogos pirazol-NQ ativos vs T. brucei | SAR revisão Rani 2022 |

#### Grupos Alifáticos — lapachol/menadiona-like

| R-group | Síntese | Racional | Referência |
|---|---|---|---|
| **–CH₃** (menadiona) | MeMgBr/CuI/THF | Vitamina K3; ativo vs T. cruzi (Docampo 1978) | Docampo et al. 1978 Z Parasitenkd |
| **–CF₃** | Togni I/II ou TMSCF₃/CsF | σp +0.54; resiste metabolismo; antiprotozoais modernos | SAR CF₃ em antiparasitários |
| **–alil** (lapachol) | Alilamina/K₂CO₃/DMF | Fragmento lapachol; Michael acceptor in vivo | Moura 2001 J Braz Chem Soc |
| **–dihidrobenzofuranila** | Oxa-Michael + ciclização (síntese LQB-118) | INSPIRADO NO LQB-118 (pterocarpanquinona Literatura) | Maia 2011; Cunha-Junior 2016 AAC |

### Como usar no dashboard

1. Selecionar botão **FBDD (R-groups)**
2. Inserir SMILES do composto parental ou clicar **Carregar composto ativo**
3. Selecionar categorias de R-groups (Ctrl+clique para múltiplas)
4. Definir número máximo de variantes (20 padrão)
5. Clicar **Gerar Candidatos**
6. Ver tabela de candidatos rankeados
7. Clicar **Ver síntese** em qualquer candidato para detalhes completos

### Interpretando os resultados FBDD

| Métrica | Bom | Moderado | Ruim |
|---|---|---|---|
| ΔG estimado | ≤ −8.0 kcal/mol | −6.0 a −8.0 | > −6.0 |
| Drug Score | ≥ 7.0 | 5.0 – 7.0 | < 5.0 |
| SAS | ≤ 3.0 (fácil) | 3.0 – 5.0 | > 5.0 |
| Lipinski | Aprovado | — | Reprovado |
| TPSA | 40–90 Å² | 90–140 Å² | > 140 Å² |

---

## Abordagem 2 — Hibridização Molecular (BRICS)

### Princípio — inspirado no LQB-118 (Literatura)

A **hibridização molecular** é uma estratégia de design amplamente usada pelo grupo De Castro et al./Literatura. O **LQB-118** (pterocarpanquinona) foi criado pela fusão de **lapachol** (naftoquinona natural) com **pterocarpano** (fragmento vegetal), resultando em:
- IC50 = 2.5 μM vs T. cruzi epimastigotas
- 99% redução carga parasitária in vivo (leishmaniose visceral, 10 mg/kg oral)
- Ativo em leucemia mielóide crônica (LMC) — atividade anticâncer

```
     Lapachol                 Pterocarpano              LQB-118
  [scaffold NQ]    +    [fragmento bicíclico]  →  [híbrido pterocarpanquinona]
  IC50 ~20 μM           — antifúngico              IC50 2.5 μM T. cruzi
                                                    IS > 40x (CC50 >100 μM VERO)
```

### Algoritmo BRICS do módulo

1. **Fragmentação**: cada composto é quebrado nos pontos de ligação definidos pelo algoritmo BRICS (Breaking Retrosynthetically Interesting Chemical Substructures — Degen et al. 2008)
2. **Filtragem**: fragmentos com < 4 átomos pesados são descartados
3. **Recombinação**: todas as combinações de fragmento-A + fragmento-B são testadas
4. **Filtragem Lipinski**: MW ≤ 650, HBD ≤ 5
5. **Pontuação e ranking**: drug_score + ΔG estimado + SAS

### Rotas sintéticas sugeridas automaticamente

O módulo detecta a química dos fragmentos e sugere a rota mais adequada:

| Situação | Rota sugerida |
|---|---|
| Fragmento B contém N | Buchwald-Hartwig (Pd₂(dba)₃/Xphos) ou SNAr |
| Fragmento B contém O | Condensação de Williamson ou oxa-Michael |
| Ambos aromáticos | Suzuki-Miyaura (Pd/boronato) ou Ullmann (Cu) |
| Sem heteroátomos | Knoevenagel ou adição Michael |

### Combinações de exemplo

Pares formados só com compostos **publicados**, dos que acompanham a plataforma. Servem
para você ver o módulo funcionando antes de usar o seu próprio scaffold.

| Composto A | Composto B | Justificativa |
|---|---|---|
| β-Lapachona | LQB-118 | Anel furanônico (ROS + bioativação por NQO1) + fragmento pterocarpano |
| Plumbagina | β-Lapachona | –OH quelante (ativo vs Fe-SOD de *Leishmania*) + bioativação por NQO1 |
| β-Lapachona | Benznidazol | Híbrido ROS + nitrorredução (dois mecanismos de fármaco aprovado) |
| Lapachol | Atovaquona | Núcleo natural + cadeia lipofílica de ligante do complexo III |

### Como usar no dashboard

1. Selecionar botão **Hibridização**
2. Inserir SMILES do Composto A (parental)
3. Inserir SMILES do Composto B (parceiro — do banco de compostos)
4. Definir número máximo de híbridos
5. Clicar **Gerar Candidatos**

---

## Abordagem 3 — Algoritmo Genético (GA)

### Princípio de evolução darwiniana aplicado à química

```
Geração 0:        Geração 5:        Geração 15:
[banco compostos] → [mutações +]  → [candidatos]
                     [cruzamentos]    [otimizados]

Fitness = f(ΔG_TcTR) × f(ADMET) × f(SAS)
         ↑ afinidade   ↑ drug-like  ↑ fácil sintetizar
```

### Operações genéticas implementadas

**Mutação (35% de probabilidade):**
Seleciona aleatoriamente uma categoria de R-group e aplica substituição no composto parental via FBDD.

**Cruzamento / Crossover (40% de probabilidade):**
Fragmenta dois compostos via BRICS e recombina fragmentos de pais diferentes (hibridização).

**Elitismo (20% da população):**
Os melhores indivíduos de cada geração são preservados sem modificação.

**Seleção por torneio:**
3 indivíduos aleatórios competem → o de maior fitness vai para reprodução.

### Função de Fitness

```
fitness = dg_score × admet_score × sas_score × lipinski_bonus

onde:
  dg_score   = max(0, (-ΔG_est - 2.0) / 8.0)      [0–1, ΔG de -2 a -10 kcal/mol]
  admet_score = drug_score / 10.0                   [0–1]
  sas_score  = 1 - (SAS - 1.0) / 9.0               [0–1, SAS de 1 a 10]
  lip_bonus  = 1.0 se Lipinski OK, 0.4 se falha
```

### Parâmetros recomendados por objetivo

| Objetivo | Gerações | Tamanho Pop. | Tempo estimado |
|---|---|---|---|
| Exploração rápida / teste | 10 | 15 | ~1–2 min |
| Triagem padrão | 15 | 20 | ~3–5 min |
| Otimização detalhada | 25 | 30 | ~10–15 min |
| Campanha de design | 50 | 50 | ~30–60 min |

### Como usar no dashboard

1. Selecionar botão **Algoritmo Genético**
2. Inserir SMILES dos compostos semente (um por linha) — recomendado: 3 compostos do banco
3. Inserir nomes correspondentes (mesma ordem)
4. Definir gerações e tamanho de população
5. Clicar **Gerar Candidatos**
6. Acompanhar a convergência pelo gráfico de evolução gerado

---

## Modelo QSAR — ΔG Estimado para TcTR

O módulo usa um modelo QSAR empírico para estimar a afinidade de ligação (ΔG) de candidatos com a **tripanotiona redutase de T. cruzi (TcTR, PDB: 1BZL)** sem precisar rodar docking para cada variante.

### Calibração do modelo

O modelo foi calibrado sobre ~50 naftoquinonas com dados de docking publicados:

| Publicação | Compostos | ΔG range |
|---|---|---|
| Pinto & De Castro 2009 | 16 NQ derivados | -4.5 a -9.2 kcal/mol |
| Pinto & De Castro 2009 | 12 lapachol derivados | -5.0 a -8.8 kcal/mol |
| Salomão et al. 2013 BMC Microbiol | 8 NQ análogos | -5.5 a -8.5 kcal/mol |
| Rani et al. 2022 Chem Biol Drug Des | 14 NQ naturais | -4.8 a -9.1 kcal/mol |

### Equação do modelo

```
ΔG_est = -5.5                         (base: 1,4-NQ simples)
       + f(logP)                       (sweet spot 2.0–3.5: -1.2)
       + f(TPSA)                       (ideal 50–90 Å²: -0.8)
       + f(HBD) × (-0.35)             (cada doador H-bond)
       + f(HBA) × (-0.18)             (até 5 aceitadores)
       + f(MW)                         (ideal 200–380 Da: -0.6)
       + N_arom × (-0.5)              (N aromático: π-stacking FAD)
       + S_presence × (-0.7)          (análogos –SCN: Cys52/57)
       + N_hal × (-0.35)              (halogênios: hal-bonds)
       + NO2_presence × (-0.5)        (ativação TcNTR)
```

### Limitações do modelo QSAR

O ΔG estimado é um **ranqueamento relativo**, não uma predição absoluta. Erros esperados: ±1.5 a 2.5 kcal/mol em relação ao docking real com Vina. Para candidatos top-5, é **fortemente recomendado** rodar o docking real via aba Pipeline do dashboard para validar.

---

## Score de Acessibilidade Sintética (SAS)

O SAS (Ertl & Schuffenhauer 2009) mede a facilidade de síntese em escala de 1 a 10.

| SAS | Interpretação | Exemplo de estrutura |
|---|---|---|
| 1.0 – 2.5 | Muito fácil — reagentes comerciais, 1–2 etapas | Menadiona, plumbagina |
| 2.5 – 3.5 | Fácil — síntese padrão, 2–4 etapas | Lapachol, juglona |
| 3.5 – 5.0 | Moderada — 4–6 etapas, proteção/desproteção | β-lapachona, LQB-118 |
| 5.0 – 7.0 | Complexa — especialização necessária | Naftoimidazóis N1/N2/N3 |
| 7.0 – 10.0 | Muito complexa — síntese total | Estruturas policíclicas avançadas |

**Critério do módulo:** `synth_feasible = True` quando SAS ≤ 4.5 (síntese acessível em laboratório de química orgânica convencional como os disponíveis no IFRJ/Literatura).

---

## Fluxo de Trabalho Recomendado

### Da hipótese ao candidato priorizado

```
Passo 1: FBDD no seu scaffold
   → Categorias: nitrogen (triazol, imidazol) + sulfur + halogen
   → Max variantes: 20
   → Objetivo: identificar 3-5 análogos com ΔG estimado <= -8.0 e SAS <= 4.0

Passo 2: Hibridização
   → scaffold × beta-lapachona   (mecanismo ROS + NQO1)
   → scaffold × LQB-118          (fragmento pterocarpano)

Passo 3: GA a partir dos melhores do banco
   → Seeds: o seu scaffold + beta-lapachona + plumbagina
   → 15 gerações, pop=20
   → Objetivo: alcançar scaffold não óbvio

Passo 4: VALIDAR -- o passo que não se pula
   → Top-3 de cada abordagem → docking real na aba Pipeline
   → Comparar com um fármaco de referência (ex.: benznidazol)
   → O ΔG do QSAR é triagem; o do docking é que entra em relatório

Passo 5: Relatório
   → Botão "Download" na aba Drug Design (EN/PT)

Passo 6: Proposta de síntese
   → Priorizar SAS <= 3.5 e ΔG <= -8.0
   → Encaminhar ao grupo de síntese parceiro
```

> O passo 4 existe porque o ΔG do modelo QSAR é uma estimativa de triagem, calibrada
> sobre dados publicados de uma classe química. Tratá-lo como afinidade é o erro que o
> módulo mais facilmente induz.

### Critérios de prioridade para síntese

Um candidato deve satisfazer **pelo menos 3 de 4** critérios para ser priorizado:

- [ ] ΔG estimado ≤ −7.5 kcal/mol (melhor que o fármaco de referência previsto)
- [ ] Drug Score ≥ 7.0 / 10
- [ ] SAS ≤ 4.0 (síntese acessível)
- [ ] Lipinski aprovado (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)

---

## Integrando Drug Design com os outros módulos

O módulo foi projetado para trabalhar em cadeia com o restante do pipeline:

```
Drug Design (novos SMILES)
       ↓
Compostos (adicionar ao banco via "Adicionar Novo Composto")
       ↓
Pipeline (rodar docking Vina para validar ΔG real)
       ↓
ADMET (perfil farmacocinético completo)
       ↓
MD Simulation (confirmar estabilidade do complexo 100 ns)
       ↓
Toxicidade (calcular IS = CC50/IC50 com os seus dados)
       ↓
Relatório (texto integrado dos módulos)
```

Para adicionar um candidato do Drug Design ao banco:
1. Copie o SMILES do candidato (coluna na tabela ou painel "Ver síntese")
2. Vá para aba **Compostos** → campo SMILES → Cole → preencha nome → **Adicionar**
3. O composto aparece no seletor do banner superior
4. Selecione-o como composto ativo e rode qualquer aba normalmente

---

## Glossário

| Termo | Definição |
|---|---|
| **FBDD** | Fragment-Based Drug Design — design de drogas baseado em fragmentos moleculares |
| **BRICS** | Breaking Retrosynthetically Interesting Chemical Substructures — algoritmo de fragmentação (Degen 2008) |
| **SAR** | Structure-Activity Relationship — relação estrutura-atividade |
| **ΔG** | Energia livre de Gibbs de ligação (kcal/mol) — quanto mais negativo, maior a afinidade |
| **ADMET** | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| **QED** | Quantitative Estimate of Drug-likeness (Bickerton 2012) — 0 a 1, quanto maior melhor |
| **SAS** | Synthetic Accessibility Score (Ertl 2009) — 1 a 10, quanto menor mais fácil sintetizar |
| **TPSA** | Topological Polar Surface Area — área polar da superfície (Å²); ideal 40–90 Å² para antiparasitários |
| **TcTR** | Tripanotiona Redutase de *T. cruzi* — alvo da campanha de exemplo (PDB 1BZL) |
| **TcNTR** | Nitroreductase de T. cruzi — ativa pró-drogas nitro (BZN, análogos –NO₂) |
| **IS** | Índice de Seletividade = CC50(hospedeiro) / IC50(parasito) |
| **CuAAC** | Copper-catalyzed Azide-Alkyne Cycloaddition — "click chemistry" (Huisgen 1,3-dipolar) |
| **SNAr** | Substituição Nucleofílica Aromática |
| **GA** | Algoritmo Genético |
| **Fitness** | Função que combina ΔG_TcTR + ADMET + SAS para guiar o GA |

---

## Referências Principais

1. De Castro et al. RFS, De Castro SL. *The double-edged sword in pathogenic trypanosomatids: The pivotal role of mitochondria in oxidative stress and bioenergetics.* Front Physiol 2014; revisão 2022 Mem IOC.

2. Pinto AV, De Castro SL. *The trypanocidal activity of naphthoquinones: A review.* Molecules 2009; 14:4570–4590.

3. De Castro et al. RFS, Henriques-Pons A, Pinto AV et al. *Effect of β-lapachone-derived naphthoimidazole on Trypanosoma cruzi.* J Antimicrob Chemother 2005;56:1034–1041.

4. Bahia MT, Bahia Barbosa C, Jardim et al. *Molecular hybridization as a powerful tool toward multitarget quinoidal systems.* MedChemComm 2016;7:1555.

5. Salomão K, De Santana NA, Molina MT, De Castro SL, De Castro et al. RFS. *Mitochondrial swelling and membrane potential collapse as primary evidence of the mode of action of NQs.* BMC Microbiol 2013;13:196.

6. Maia RC et al. *LQB-118, a pterocarpanquinone structurally related to lapachol: a novel class of agent with high apoptotic effect in chronic myeloid leukemia cells.* Invest New Drugs 2011;29:1143.

7. Rani I, Bhatt A. *Natural naphthoquinones and their derivatives as potential drug molecules against trypanosome parasites.* Chem Biol Drug Des 2022;100:699–724.

8. Degen J et al. *On the art of compiling and using "drug-like" chemical fragment spaces.* ChemMedChem 2008;3:1503.

9. Ertl P, Schuffenhauer A. *Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions.* J Cheminform 2009;1:8.
