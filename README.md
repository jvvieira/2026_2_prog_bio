# Trabalho final — Associações genéticas da fibromialgia (GWAS Catalog)

Análise exploratória das associações genéticas de **fibromialgia** (`MONDO_0005546`)
publicadas no [GWAS Catalog](https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005546),
usando a API REST v2 como fonte de dados.

O script principal do projeto é o [`trabalho-final.py`](trabalho-final.py), que baixa
(quando habilitado), carrega, normaliza, limpa e visualiza os dados. A descrição
detalhada de cada campo devolvido pela API está em
[`dicionario_api_gwas.md`](dicionario_api_gwas.md).

## Arquivos do projeto

| Arquivo | Descrição |
|---|---|
| `trabalho-final.py` | Script principal: carga, normalização, limpeza e geração dos gráficos. |
| `data.json` | Resposta bruta da API do GWAS Catalog (215 associações). |
| `dicionario_api_gwas.md` | Dicionário de dados da API v2, com tipos, descrições e observações sobre valores ausentes. |
| `figuras/` | Gráficos gerados pelo script (PNG, 300 dpi). |
| `projeto-final.pdf` | Enunciado do trabalho. |
| `docs/` | Specs de desenvolvimento das etapas do projeto. |

## Requisitos

- Python 3.10 ou superior
- [`pandas`](https://pandas.pydata.org/)
- [`matplotlib`](https://matplotlib.org/)

Exemplo de preparação do ambiente:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pandas matplotlib
```

## Como executar

Execute o script a partir da raiz do projeto, pois ele abre `data.json` pelo
caminho relativo e grava os gráficos em `figuras/`:

```bash
python3 trabalho-final.py
```

O script imprime no terminal as dimensões e colunas do DataFrame, o relatório de
colunas removidas, a interpretação de cada gráfico e grava as figuras em
`figuras/`.

### Download dos dados

O bloco de download via API está comentado no início do arquivo (linhas 24–27)
para não sobrecarregar o serviço. Para baixar novamente o `data.json`, basta
descomentar o trecho:

```python
with urllib.request.urlopen(
    "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?size=300&efo_id=MONDO_0005546"
) as url:
    data = json.loads(url.read().decode())
    with open("data.json", "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
```

## O que o script faz

O código é organizado em três etapas sequenciais, executadas em `__main__`.

### Etapa 1 — Carga e normalização (`etapa_1`)

1. Lê `data.json` com `json.load`.
2. Converte a lista `_embedded.associations` em um DataFrame com
   `pd.json_normalize`.
3. Cria três colunas derivadas para facilitar a análise:
   - `_mapped_genes` — genes mapeados unidos em uma única string separada por vírgulas;
   - `_allele` — `rs_id` do primeiro SNP da lista `snp_allele`;
   - `_cromossome` — cromossomo extraído de `locations` no formato `cromossomo:posição`.
4. Imprime as dimensões e os nomes das colunas do DataFrame.

### Etapa 2 — Remoção de colunas sem variação (`etapa_2`)

1. Identifica colunas **completamente vazias** e colunas **constantes** com
   `identificar_colunas_sem_variancia`.
2. Trata como ausentes: `NaN`, `None`, string vazia, apenas espaços, `"-"`, `"NR"`,
   listas/dicionários vazios.
3. Exibe um relatório com as colunas encontradas e um exemplo de valor de cada uma.
4. Retorna uma **cópia** do DataFrame sem essas colunas, sem alterar o original.

A remoção é necessária porque a base é composta por associações genéticas já
publicadas: não é possível inferir os valores ausentes sem inventar informação.

### Etapa 3 — Visualização e interpretação (`etapa_3`)

Gera três gráficos e imprime uma interpretação curta de cada um:

| Função | Figura gerada | Pergunta |
|---|---|---|
| `grafico_cromossomos` | `figuras/01_cromossomos.png` | Como as associações se distribuem entre os cromossomos? |
| `grafico_genes` | `figuras/02_genes.png` | Quais genes aparecem com mais frequência nas associações? |
| `grafico_beta_cromossomo` | `figuras/03_beta_por_cromossomo.png` | Qual é o beta mediano por cromossomo? |

Funções auxiliares usadas nessa etapa:

- `_salvar_grafico(fig, nome_arquivo)` — cria `figuras/` se necessário, salva em
  300 dpi e fecha a figura;
- `_ordenar_cromossomos(cromossomos)` — ordena `1`–`22`, `X`, `Y`, `MT` e deixa
  valores inesperados no final;
- `_beta_com_sinal(valor)` — converte textos como `"0.06 unit increase"` em
  `0.06` e `"0.06 unit decrease"` em `-0.06`;
- `_valor_ausente` e `_valor_canonico` — normalizam a detecção e a comparação de
  valores ausentes e de listas/dicionários.

O gráfico de beta usa apenas as associações cujo campo `beta` está no formato
textual com número, direção (`increase`/`decrease`) e unidade. Os demais valores
(por exemplo, `"-"` ou `"NR"`) são descartados dessa figura.

## Painel de funções

| Função | Entrada | Saída |
|---|---|---|
| `etapa_1()` | — | `DataFrame` normalizado. |
| `etapa_2(df)` | `DataFrame` da Etapa 1 | Cópia sem colunas vazias/constantes. |
| `etapa_3(df)` | `DataFrame` da Etapa 2 | Gera os três gráficos e imprime interpretações. |
| `identificar_colunas_sem_variancia(df)` | `DataFrame` | Tupla `(colunas_vazias, colunas_constantes)`. |
| `grafico_cromossomos(df)` | `DataFrame` | Salva `01_cromossomos.png`. |
| `grafico_genes(df)` | `DataFrame` | Salva `02_genes.png`. |
| `grafico_beta_cromossomo(df)` | `DataFrame` | Salva `03_beta_por_cromossomo.png`. |

## Referência dos dados: `dicionario_api_gwas.md`

O arquivo [`dicionario_api_gwas.md`](dicionario_api_gwas.md) documenta a
estrutura da resposta da API v2 do GWAS Catalog, incluindo:

- contexto do dataset (215 associações, 10 estudos, 3 publicações);
- tabela com todos os campos da API, tipo, descrição oficial e presença no `data.json`;
- subcampos de `_links`;
- observações importantes, como:
  - `"-"`, `"NR"`, `""` e `[]` são **códigos de ausência diferentes**;
  - a API devolve listas em `mapped_genes`, `locations`, `snp_allele`, entre outros;
  - `or_per_copy_num` é a versão numérica de `or_value`;
  - `beta` é texto e precisa ser convertido para número;
  - `ci_lower`/`ci_upper` são os limites numéricos do IC95%;
  - `association_id` é a chave única da associação.

Consulte esse dicionário antes de alterar o script ou de interpretar colunas
derivadas como `_mapped_genes`, `_allele` e `_cromossome`.

## Limitações

- Não há imputação de valores ausentes: colunas sem variação são removidas.
- A base contém apenas associações já reportadas/publicadas, em grande parte
  significativas em nível genômico, o que limita comparações de significância.
- OR e beta estão disponíveis apenas para parte das associações.
- A direção do efeito (OR < 1 ou beta negativo) depende do alelo reportado e
  **não** deve ser interpretada como proteção ou risco clínico.
- Os gráficos agregam estudos diferentes; o script não realiza uma meta-análise.

## Referências

- [GWAS Catalog — fibromialgia (MONDO_0005546)](https://www.ebi.ac.uk/gwas/efotraits/MONDO_0005546)
- [GWAS Catalog REST API v2 — documentação](https://www.ebi.ac.uk/gwas/rest/api/v2/docs)
- [GWAS Catalog REST API v2 — referência](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference)
- [Especificação OpenAPI (YAML)](https://www.ebi.ac.uk/gwas/rest/api/v2/rest-api-doc.yaml)
- [Documentação dos cabeçalhos do TSV](https://www.ebi.ac.uk/gwas/docs/fileheaders)
- Dicionário de dados do projeto: [`dicionario_api_gwas.md`](dicionario_api_gwas.md)
- Enunciado do trabalho: [`projeto-final.pdf`](projeto-final.pdf)
