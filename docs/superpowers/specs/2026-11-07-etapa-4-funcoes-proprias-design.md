# Spec — Etapa 4: duas funções próprias e aplicação ao dataset

**Data:** 2026-11-07  
**Projeto:** Análise de associações genéticas de fibromialgia (GWAS Catalog, `MONDO_0005546`)  
**Etapa do PDF:** Etapa 4 — Desenvolvimento e aplicação de funções

## Objetivo

Implementar e demonstrar **duas funções próprias** aplicadas efetivamente ao dataset tratado na Etapa 2. As funções devem ter finalidade clara, parâmetros documentados, resultados reproduzíveis e uso prático no contexto da fibromialgia.

A entrega desta etapa inclui:

1. o código das duas funções;
2. uma demonstração de aplicação ao dataset real;
3. documentação com parâmetros, resultados e interpretação;
4. testes automatizados das funções.

## Função 1 — `resumo_descritivo(df, coluna)`

### Finalidade

Produzir um resumo descritivo padronizado de qualquer coluna do dataset, tratando corretamente os códigos de ausência usados pelo GWAS Catalog. A função evita repetir manualmente a mesma lógica de contagem, estatísticas e valores ausentes em cada análise.

### Entradas

- `df` (`pandas.DataFrame`): dataset a analisar.
- `coluna` (`str`): nome da coluna a resumir.

### Comportamento

- Validar que `df` é um `DataFrame` e que `coluna` existe; caso contrário, levantar `ValueError` com mensagem clara.
- Considerar ausentes: `NaN`, `None`, string vazia ou só com espaços, `"-"`, `"NR"`, lista vazia e dicionário vazio.
- Calcular, para todas as colunas:
  - total de linhas;
  - quantidade e percentual de valores presentes;
  - quantidade e percentual de valores ausentes;
  - quantidade de valores únicos entre os presentes.
- Para colunas numéricas, calcular também:
  - média, mediana, desvio-padrão, mínimo, primeiro quartil (Q1), terceiro quartil (Q3) e máximo.
- Para colunas categóricas, de texto ou listas, calcular também:
  - valor mais frequente (moda), sua contagem e seu percentual;
  - tabela completa de contagens por valor presente.
- Preservar os tipos originais; não converter nem imputar valores.

### Saída

Um dicionário com as chaves:

- `coluna` (`str`);
- `tipo` (`str`: `"numérica"`, `"categórica"` ou `"lista"`);
- `n_total` (`int`);
- `n_presentes` (`int`);
- `n_ausentes` (`int`);
- `percentual_ausentes` (`float`, de 0 a 100);
- `n_unicos` (`int`);
- `estatisticas` (`dict` ou `None`): presentes apenas para colunas numéricas;
- `contagens` (`pandas.Series` ou `None`): presentes apenas para colunas categóricas/listas.

### Aplicação prevista no dataset

- `or_per_copy_num`: 50 valores presentes, 165 ausentes (76,7%), mediana 0,9615, mínimo 0,04 e máximo 2,779.
- `p_value`: 215 valores presentes, 48 únicos, mediana 5e-9 e máximo 5e-8.
- `risk_frequency`: 49 valores numéricos presentes e 166 `"NR"` (77,2% ausentes).
- `mapped_genes`: 137 listas não vazias e 78 listas vazias (36,3% ausentes).

## Função 2 — `classificar_efeito_or(df, coluna_or="or_per_copy_num", limiar=1.0)`

### Finalidade

Classificar cada associação genética conforme a direção do odds ratio em relação ao efeito nulo (`OR = 1`), criando uma nova variável categórica reutilizável para contagens, proporções e comparações. A função atende ao exemplo do PDF de “classificar em categorias com base em critérios definidos”.

### Entradas

- `df` (`pandas.DataFrame`): dataset a classificar.
- `coluna_or` (`str`, padrão `"or_per_copy_num"`): coluna com o odds ratio.
- `limiar` (`float`, padrão `1.0`): valor de referência para classificar o efeito.

### Comportamento

- Validar que `df` é um `DataFrame`, que `coluna_or` existe e que `limiar` é numérico e finito.
- Converter a coluna de OR para numérico com `pd.to_numeric(..., errors="coerce")`, sem alterar a coluna original.
- Criar uma cópia do DataFrame com a nova coluna `efeito_genetico`, contendo:
  - `"Protetor (OR < limiar)"` quando `OR < limiar`;
  - `"Risco (OR > limiar)"` quando `OR > limiar`;
  - `"Neutro (OR = limiar)"` quando `OR == limiar`;
  - `"Sem OR"` quando o valor é ausente ou não numérico.
- Não remover linhas e não alterar o DataFrame recebido.
- Não inferir OR ausente.

### Saída

Um novo `pandas.DataFrame` com todas as colunas originais mais `efeito_genetico`.

### Aplicação prevista no dataset

- 165 associações (76,7%) ficam em `"Sem OR"`.
- 28 associações (13,0% do total; 56,0% das 50 com OR) ficam em `"Protetor (OR < 1)"`.
- 22 associações (10,2% do total; 44,0% das 50 com OR) ficam em `"Risco (OR > 1)"`.
- Nenhuma associação tem `OR = 1`.

## Demonstração — `etapa_4(df)`

A etapa deve:

1. aplicar `resumo_descritivo` a `or_per_copy_num`, `p_value`, `risk_frequency` e `mapped_genes`;
2. aplicar `classificar_efeito_or` ao DataFrame tratado;
3. exibir a contagem e o percentual de cada categoria de `efeito_genetico`;
4. imprimir uma interpretação curta, incluindo limitações;
5. retornar um dicionário com os resultados, para reutilização.

A demonstração não deve alterar as Etapas 1, 2 e 3 nem os gráficos já gerados.

## Interpretação esperada

- A maior parte das associações (165/215) não possui OR disponível; qualquer conclusão sobre direção de efeito fica restrita às 50 com OR.
- Entre as 50 com OR, 28 têm OR < 1 e 22 têm OR > 1.
- O OR < 1 deve ser descrito como direção do efeito em relação ao alelo reportado, **não** como prova de proteção clínica, pois a codificação do alelo pode variar entre estudos.
- 207 das 215 associações têm `p_value < 5e-8`; o dataset é composto majoritariamente por associações já significativas em nível genômico, o que limita comparações de significância.
- A ausência de OR, IC95% e genes mapeados é frequente e deve ser explicitada como limitação.

## Documentação

Criar `docs/etapa-4-funcoes.md` contendo:

- finalidade de cada função;
- explicação dos parâmetros e do valor de retorno;
- trecho de código de cada função;
- chamadas de exemplo e resultados reais no dataset;
- interpretação no contexto da fibromialgia;
- limitações.

Atualizar/criar `README.md` com:

- descrição resumida do projeto;
- instruções para executar `python3 trabalho-final.py`;
- descrição das duas funções da Etapa 4;
- indicação de que os gráficos são salvos em `figuras/`.

## Testes

Criar testes automatizados com `pytest` cobrindo:

- resumo numérico com valores presentes e ausentes (`"NR"`, `"-"`, `""`, `NaN`);
- resumo categórico com contagens e moda;
- erro para coluna inexistente;
- classificação de OR abaixo, acima, igual e ausente;
- preservação do DataFrame original;
- funcionamento do parâmetro `limiar`.

## Fora de escopo

- Imputação de valores ausentes.
- Alterações nas Etapas 1, 2 e 3.
- Novos gráficos.
- Modelagem estatística ou inferência causal.
- Consolidação do notebook final da disciplina (será feita em etapa posterior).

## Critérios de aceitação

- `python3 trabalho-final.py` executa sem erros, mantém os cinco gráficos e imprime a saída da Etapa 4.
- As duas funções são aplicadas ao dataset real e produzem os resultados esperados.
- A documentação contém código, parâmetros, resultados e interpretação.
- Os testes passam.
- O `README.md` explica como executar o projeto e onde estão as funções.
