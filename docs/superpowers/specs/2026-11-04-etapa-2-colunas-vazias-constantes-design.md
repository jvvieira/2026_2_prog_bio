# Spec — Etapa 2: remoção de colunas vazias ou constantes

**Data:** 2026-11-04  
**Projeto:** Análise de associações genéticas de fibromialgia (GWAS Catalog, `MONDO_0005546`)  
**Etapa do PDF:** Etapa 2 — Avaliação geral, manipulação e normalização dos dados

## Objetivo

Identificar, no DataFrame carregado na Etapa 1, quais colunas estão **completamente vazias** ou **não apresentam variação** (todos os valores são iguais) e retornar uma cópia do DataFrame sem essas colunas.

## Definições de comportamento

Uma coluna é considerada **completamente vazia** quando todos os seus valores são ausentes. São tratados como ausentes:

- `NaN` / `None`;
- string vazia ou composta apenas por espaços;
- lista vazia;
- dicionário vazio;
- os códigos `-` e `NR`, usados pelo GWAS Catalog para “não disponível” e “não reportado”.

Uma coluna é considerada **constante** quando, após ignorar os valores ausentes:

- existe pelo menos um valor não ausente; e
- todos os valores não ausentes são iguais entre si.

Listas e dicionários são comparados pelo seu conteúdo, não pela identidade do objeto. A ordem dos elementos de uma lista é considerada parte do valor.

## Resultado esperado

A Etapa 2 deve:

1. Exibir um relatório no console com:
   - a dimensão do DataFrame antes da remoção;
   - a lista de colunas completamente vazias;
   - a lista de colunas constantes;
   - para cada coluna removida, o motivo (`vazia` ou `constante`) e um exemplo de valor;
   - a dimensão do DataFrame depois da remoção.
2. Retornar um **novo DataFrame**, sem as colunas identificadas, preservando o DataFrame original.
3. Não remover colunas que tenham variação real, mesmo que tenham muitos valores ausentes.

## Critérios de aceitação

- No dataset atual, devem ser identificadas como completamente vazias: `bg_efo_traits` e `pvalue_description`.
- No dataset atual, devem ser identificadas como constantes: `multi_snp_haplotype`, `snp_interaction` e `efo_traits`.
- O DataFrame retornado deve ter **23 colunas** (28 originais menos 5 removidas) e **215 linhas**.
- As colunas `locations`, `mapped_genes`, `_mapped_genes`, `risk_frequency`, `beta`, `range`, `or_value`, `ci_lower` e `ci_upper` **não** devem ser removidas, pois variam.
- O DataFrame original deve permanecer com as 28 colunas.

## Fora de escopo

- Imputação de valores ausentes.
- Normalização Min-Max ou padronização Z-score.
- Conversão de tipos das colunas restantes.
- Análise de duplicatas.
- Gráficos.

Esses itens podem ser tratados em outras partes da Etapa 2 ou em etapas posteriores.

## Observações

- A Etapa 1 adiciona a coluna derivada `_mapped_genes` (genes separados por vírgula), por isso o DataFrame recebido pela Etapa 2 tem 28 colunas, e não 27.
- `efo_traits` é constante porque todo o dataset foi filtrado para fibromialgia (`MONDO_0005546`). A remoção é consequência da regra, não uma indicação de que o traço seja irrelevante.
- `bg_efo_traits` é uma lista vazia em todas as linhas.
- `pvalue_description` é uma string vazia em todas as linhas.
- `multi_snp_haplotype` e `snp_interaction` são `False` em todas as linhas.
