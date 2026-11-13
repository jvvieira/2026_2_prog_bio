# Dicionário de dados — GWAS Catalog API v2

## Integrantes do Grupo

**João Victor Rossetti Vieira** (RA 202630554)

**Viviane** 

##
**Fonte:** [GWAS Catalog REST API v2](https://www.ebi.ac.uk/gwas/rest/api/v2/docs)  
**Traço consultado:** fibromialgia (`MONDO_0005546`)  
**Arquivos de referência:** `data.json` (API)   
**Data da consulta:** 2026-11-04

## Contexto do dataset

- **215 associações**
- **10 estudos** (`accession_id`)
- **3 publicações** (`pubmed_id`: `40899795`, `42521668`, `42521817`)
- **50** associações com OR
- **71** associações com beta
- **94** associações sem OR nem beta
- **98** associações com IC95%
- **166** associações com RAF = `NR`
- **78** associações sem gene mapeado
- **67** associações sem localização genômica

## Dicionário da API

| Campo na API | Tipo | Descrição oficial (OpenAPI) | Equivalente no TSV | Presente no `data.json`? |
|---|---|---|---|---|
| `association_id` | integer (int64) | ID primário da tabela de associações | — | Sempre; 215 valores únicos |
| `risk_frequency` | string | Frequência do alelo de risco/efeito associado ao SNP mais forte, em controles | `riskFrequency` | Sempre (`NR` ou número como string) |
| `pvalue_description` | string | Informação descrevendo o contexto do p-valor (ex.: “females, smokers”) | `pValueAnnotation` | Sempre; string vazia neste dataset |
| `pvalue_mantissa` | integer (int32) | Mantissa do p-valor reportado para o SNP de risco mais forte | — (derivável de `pValue`) | Sempre |
| `pvalue_exponent` | integer (int32) | Expoente do p-valor reportado para o SNP de risco mais forte | — (derivável de `pValue`) | Sempre |
| `multi_snp_haplotype` | boolean | Indica se a associação é de um haplótipo de múltiplos SNPs | — | Sempre; todos `false` |
| `snp_interaction` | boolean | Indica se a associação é de interação SNP-SNP | — | Sempre; todos `false` |
| `snp_type` | string | Indica se o SNP já foi reportado anteriormente (ex.: `known`) | — | Não apareceu no `data.json` |
| `range` | string | Intervalo de confiança de 95% | `ci` | Sempre (`-` quando ausente) |
| `description` | string | Comentário adicional relacionado ao valor de beta ou OR (ex.: `Discovery`) | — | Não apareceu no `data.json` |
| `or_per_copy_num` | number (float) | Valor numérico do odds ratio | `orValue` (versão numérica) | 50 registros |
| `or_value` | string | Odds ratio | `orValue` | 50 registros |
| `beta_num` | number (float) | Valor numérico do coeficiente beta | — (derivável de `beta`) | Não apareceu no `data.json` |
| `beta_unit` | string | Unidade do coeficiente beta | — (derivável de `beta`) | Não apareceu no `data.json` |
| `beta_direction` | string | Direção do coeficiente beta (`increase`/`decrease`) | — (derivável de `beta`) | Não apareceu no `data.json` |
| `beta` | string | Texto concatenado com número, direção e unidade do beta | `beta` | 71 registros |
| `last_mapping_date` | string (date-time) | Última data em que a associação foi mapeada para o Ensembl | — | Não apareceu no `data.json` |
| `last_update_date` | string (date-time) | Última data em que a associação foi atualizada | — | Não apareceu no `data.json` |
| `p_value` | number (double) | P-valor reportado para a associação | `pValue` | Sempre |
| `efo_traits` | array de objetos | Trait(s) EFO da associação; cada item tem `efo_id` e `efo_trait` | `efoTraits` (apenas o rótulo) | Sempre; 1 item por registro |
| `reported_trait` | array de strings | Trait(s) reportado(s) pelo estudo | `traitName` | Sempre; 1 item por registro |
| `accession_id` | string | Accession ID do estudo no GWAS Catalog (`GCST...`) | `accessionId` | Sempre; 10 estudos |
| `locations` | array de strings | Localização(ões) genômica(s) do SNP, no formato `cromossomo:posição` | `locations` | Sempre; lista vazia quando não mapeado |
| `mapped_genes` | array de strings | Gene(s) sobreposto(s) ao variant; se intergênico, os genes 5′ e 3′ mais próximos | `mappedGenes` | Sempre; lista vazia quando não mapeado |
| `snp_effect_allele` | array de strings | SNP(s) e alelo de risco mais fortemente associado(s) no locus | `riskAllele` | Sempre; 1 item por registro |
| `snp_allele` | array de objetos | SNP(s) e alelo de risco de forma estruturada; cada item tem `rs_id` e `effect_allele` | — (derivável de `riskAllele`) | Sempre; 1 item por registro |
| `bg_efo_traits` | array de objetos | Trait(s) de background compartilhado(s) por todos os indivíduos do GWAS | `bgTraits` (apenas o rótulo) | Sempre; lista vazia neste dataset |
| `pubmed_id` | string | PMID da publicação | `pubmedId` | Sempre; 3 publicações |
| `first_author` | string | Primeiro autor da publicação | `author` | Sempre; 3 autores |
| `ci_lower` | number (float) | Limite inferior do IC95% | `ci` (parte inferior) | 98 registros |
| `ci_upper` | number (float) | Limite superior do IC95% | `ci` (parte superior) | 98 registros |
| `links` | array de objetos | Links HATEOAS do recurso | — | Não como `links`; na prática vem como `_links` |

### Subcampos de `_links`

| Campo na API | Descrição | Equivalente no TSV |
|---|---|---|
| `_links.self.href` | URL do registro da associação na API | — |
| `_links.loci.href` | URL do endpoint de loci da associação | — |
| `_links.snp.href` | URL do SNP no GWAS Catalog | — |

## Observações importantes

1. **`-`, `NR`, `""` e `[]` são códigos de ausência diferentes:**
   - `-` = não disponível/não aplicável;
   - `NR` = *not reported* (não reportado no artigo);
   - `""` = texto vazio (ex.: `pvalue_description`);
   - `[]` = lista vazia (ex.: `mapped_genes`, `locations`, `bg_efo_traits`).
2. **A API usa listas** em `efo_traits`, `reported_trait`, `locations`, `mapped_genes`, `snp_effect_allele`, `snp_allele` e `bg_efo_traits`. O TSV achata esses valores em strings.
3. **`or_value` é string; `or_per_copy_num` é o mesmo OR em formato numérico.** Use `or_per_copy_num` para cálculos.
4. **`beta` é texto** com número, direção e unidade (ex.: `0.035 unit increase`). Os campos `beta_num`, `beta_unit` e `beta_direction` estão documentados no OpenAPI, mas não vieram no `data.json`.
5. **`range` é o IC95% em texto; `ci_lower` e `ci_upper` são os limites numéricos.** Use os campos numéricos para cálculos.
6. **`p_value` já é numérico.** `pvalue_mantissa` e `pvalue_exponent` permitem reconstruir o valor exibido (ex.: `4` e `-8` → `4 × 10⁻⁸`).
7. **`association_id` é a chave única da associação** (215 valores únicos). `accession_id` identifica o estudo, não o SNP.
8. **`efo_traits` acrescenta o `efo_id`** (`MONDO_0005546`) ao rótulo `fibromyalgia`, informação que o TSV não traz.
9. **`mapped_genes` não é evidência de causalidade**; é mapeamento posicional. A ordem dos genes pode diferir da ordem do TSV.
10. **Há uma localização duplicada na API** para `rs9279546-CT`: `["6:32264119", "6:32264119"]`. O TSV mostra apenas uma ocorrência.
11. **P-valores são arredondados para 1 dígito significativo** pelo GWAS Catalog.
12. **OR < 1 é possível** em associações adicionadas após 2021; o Catalog não inverte mais os ORs para ficarem sempre > 1.
13. **`reported_trait = "Fibromyalgia (MTAG)"`** indica resultado de *Multi-Trait Analysis of GWAS*, não uma coorte independente.
14. **A consulta por `efo_trait=fibromyalgia` funciona**, mas para o termo exato é mais preciso usar `efo_id=MONDO_0005546`. A consulta `efo_trait=MONDO_0005546` retorna 0 resultados.

## Links úteis

- **Swagger UI / API Reference:** https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference
- **Swagger UI direto:** https://www.ebi.ac.uk/gwas/rest/api/v2/swagger-ui/index.html
- **Especificação OpenAPI (YAML):** https://www.ebi.ac.uk/gwas/rest/api/v2/rest-api-doc.yaml
- **Schema `associations`:** seção `components.schemas.associations` do YAML
- **Documentação dos cabeçalhos do TSV:** https://www.ebi.ac.uk/gwas/docs/fileheaders
- **Página do traço:** https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005546
- **Metadados da API:** https://www.ebi.ac.uk/gwas/rest/api/v2/metadata
