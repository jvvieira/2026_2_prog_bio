# Spec — Etapa 3: visualização e interpretação dos dados

**Data:** 2026-11-04  
**Projeto:** Análise de associações genéticas de fibromialgia (GWAS Catalog, `MONDO_0005546`)  
**Etapa do PDF:** Etapa 3 — Visualização e interpretação dos dados

## Objetivo

Produzir quatro gráficos diferentes a partir do DataFrame tratado na Etapa 2, cada um com:

- pergunta ou objetivo da visualização;
- justificativa do tipo de gráfico;
- código que gera o gráfico;
- descrição dos principais padrões observados;
- interpretação no contexto da fibromialgia, incluindo limitações.

## Dados de entrada

O DataFrame resultante da Etapa 2 tem 215 linhas e 23 colunas. As colunas usadas nesta etapa são:

- `locations`;
- `mapped_genes`;
- `or_per_copy_num` ou `or_value`;
- `p_value`;
- `reported_trait`.

## Gráfico 1 — Prevalência de localização no genoma

**Pergunta:** como as associações estão distribuídas entre os cromossomos?

**Comportamento:**

- Extrair o cromossomo de cada valor de `locations` no formato `cromossomo:posição`.
- Contar cada associação uma única vez por cromossomo, mesmo que a lista de localizações tenha entradas repetidas.
- Desconsiderar associações sem localização e informar quantas foram excluídas.
- Ordenar os cromossomos em ordem natural (`1` a `22`, depois `X` e `Y`, quando existirem).
- Gerar um gráfico de barras vertical com:
  - título “Distribuição das associações por cromossomo”;
  - eixo X “Cromossomo”;
  - eixo Y “Número de associações”;
  - valor numérico sobre cada barra;
  - nota informando quantas associações tinham localização e quantas não tinham.

**Justificativa:** o cromossomo é uma variável categórica e o gráfico de barras permite comparar frequências entre categorias.

**Resultado esperado no dataset atual:** 148 associações com localização, 67 sem localização, 19 cromossomos representados. O cromossomo 6 deve concentrar o maior número de associações.

## Gráfico 2 — Prevalência de genes e recorrência entre associações

**Pergunta:** quais genes aparecem com mais frequência e quais aparecem em mais de uma associação?

**Comportamento:**

- Expandir a coluna `mapped_genes`, separando os genes por vírgula quando necessário.
- Remover genes vazios.
- Contar cada gene uma única vez por associação, para não inflar a contagem quando o mesmo gene aparece mais de uma vez na mesma linha.
- Identificar os 15 genes mais frequentes.
- Gerar um gráfico de barras horizontal com:
  - título “Genes mais frequentes nas associações”;
  - eixo X “Número de associações”;
  - eixo Y “Gene”;
  - valor numérico ao lado de cada barra;
  - cor diferente para genes que aparecem em mais de uma associação.
- Exibir também um resumo com:
  - número total de ocorrências de genes;
  - número de genes únicos;
  - número de genes que aparecem em mais de uma associação;
  - lista dos genes recorrentes e suas contagens.

**Justificativa:** a contagem por gene é uma frequência de uma variável categórica; o gráfico de barras horizontais facilita a leitura de nomes longos.

**Resultado esperado no dataset atual:** 199 ocorrências, 131 genes únicos, 50 genes em mais de uma associação. Os genes mais frequentes devem ter 4 ocorrências.

## Gráfico 3 — Distribuição dos odds ratios

**Pergunta:** como os odds ratios das associações se distribuem em relação a OR = 1?

**Comportamento:**

- Usar apenas as associações com OR numérico disponível.
- Gerar um histograma com 15 bins.
- Adicionar uma linha vertical de referência em OR = 1.
- Colorir as barras conforme o OR seja menor que 1 ou maior que 1.
- Exibir título “Distribuição dos odds ratios” e eixos “Odds ratio” e “Frequência”.
- Informar no console:
  - quantidade de ORs disponíveis;
  - quantidade de ORs < 1;
  - quantidade de ORs > 1;
  - mediana.

**Justificativa:** o OR é uma variável numérica contínua; o histograma mostra a forma da distribuição e a posição em relação ao efeito nulo.

**Resultado esperado no dataset atual:** 50 ORs disponíveis; 28 abaixo de 1, 22 acima de 1; mediana aproximadamente 0,96.

## Gráfico 4 — Dispersão entre odds ratio e significância

**Pergunta:** existe relação entre a magnitude do odds ratio e a significância estatística (`-log10(p)`)?

**Comportamento:**

- Usar apenas as associações com OR e p-valor disponíveis.
- Calcular `-log10(p_value)`.
- Gerar um gráfico de dispersão com OR no eixo X e `-log10(p)` no eixo Y.
- Colorir os pontos por `reported_trait` (“Fibromyalgia” ou “Fibromyalgia (MTAG)”).
- Adicionar:
  - linha vertical em OR = 1;
  - linha horizontal em `-log10(5e-8)` (limiar de significância genômica);
  - legenda para os traits;
  - título “Relação entre odds ratio e significância”;
  - eixos “Odds ratio” e “-log10(p)”.
- Calcular e exibir a correlação de Pearson e de Spearman entre OR e `-log10(p)`.
- Não interpretar correlação como causalidade.

**Justificativa:** a dispersão é adequada para investigar a relação entre duas variáveis numéricas.

**Resultado esperado no dataset atual:** 50 pontos. A correlação deve ser calculada e interpretada com cautela, pois o dataset contém apenas associações significativas e mistura resultados de estudos diferentes.

## Gráfico 5 — Forest plot dos odds ratios

**Pergunta:** quais variantes/genes apresentam efeito protetor ou de risco, e qual a incerteza desses efeitos?

**Comportamento:**

- Usar apenas associações com OR e IC95% completos (`ci_lower` e `ci_upper`).
- Extrair de `snp_allele` o `rs_id` e o `effect_allele`.
- Usar `mapped_genes` para montar o rótulo de cada linha no formato `rsID (alelo) — gene`.
- Ordenar as linhas pelo OR.
- Gerar um forest plot com:
  - OR no eixo X, em escala logarítmica;
  - uma linha por variante no eixo Y;
  - ponto no OR e barra horizontal representando o IC95%;
  - linha vertical de referência em OR = 1;
  - cores diferentes para OR < 1 e OR > 1;
  - título “Forest plot dos odds ratios com IC95%”;
  - eixos “Odds ratio (escala log)” e “Variante (alelo de efeito) — gene mapeado”.
- Informar quantas associações com OR foram omitidas por não terem IC95% completo.

**Justificativa:** o forest plot mostra simultaneamente a direção do efeito, a magnitude e a incerteza, sendo mais informativo para levantar hipóteses de alvos do que uma contagem de alelos.

**Resultado esperado no dataset atual:** 27 associações com OR e IC95% completo; 23 associações com OR foram omitidas por falta de IC.

## Saída

- Salvar cada gráfico como PNG na pasta `figuras/`, com 300 dpi e layout ajustado.
- Nomes dos arquivos:
  - `01_localizacoes_por_cromossomo.png`;
  - `02_genes_mais_frequentes.png`;
  - `03_distribuicao_odds_ratio.png`;
  - `04_dispersao_or_pvalor.png`;
  - `05_forest_plot_or.png`.
- Exibir cada gráfico ao final da construção.
- Imprimir, para cada gráfico, um bloco curto de interpretação com os principais padrões observados.

## Fora de escopo

- Testes estatísticos inferenciais além das correlações do Gráfico 4.
- Modelagem, regressão ou ajuste de múltiplas comparações.
- Alterações nas Etapas 1 e 2.
- Gráficos interativos.

## Critérios de aceitação

- O script executa sem erros a partir do `data.json`.
- As cinco figuras são geradas e salvas em `figuras/`.
- Associações sem localização, sem gene ou sem OR são excluídas apenas dos gráficos que dependem desses campos, com a quantidade informada no console.
- Nenhum gráfico usa `efo_traits`, `bg_efo_traits`, `pvalue_description`, `multi_snp_haplotype` ou `snp_interaction`, pois essas colunas foram removidas na Etapa 2.
