import math


ROTULOS_GRAFICO_3 = {
    "titulo": "Média do beta com sinal por rs_id",
    "xlabel": "Beta médio (increase > 0, decrease < 0) — escala symlog",
    "ylabel": "rs_id (SNP)",
}


def rs_id_de_snp_allele(valores):
    """Extrai o rs_id do primeiro item de snp_allele."""
    if not isinstance(valores, list) or not valores:
        return None

    primeiro = valores[0]
    if not isinstance(primeiro, dict):
        return None

    return primeiro.get("rs_id")


def resumo_beta_por_rsid(pares):
    """Resume (snp_allele, beta) como (rs_id, media, n) ordenado por |media|."""
    return resumo_beta_por_allele(
        (rs_id_de_snp_allele(snp_allele), beta) for snp_allele, beta in pares
    )


def _e_nan(valor):
    try:
        return math.isnan(valor)
    except TypeError:
        return False


def resumo_beta_por_allele(pares):
    """Agrupa pares (allele, beta) e devolve (allele, media, n) por |media|."""
    grupos = {}

    for allele, beta in pares:
        if allele is None or beta is None or _e_nan(beta):
            continue
        grupos.setdefault(str(allele), []).append(beta)

    resumo = [
        (allele, sum(valores) / len(valores), len(valores))
        for allele, valores in grupos.items()
    ]

    return sorted(resumo, key=lambda linha: abs(linha[1]), reverse=True)
