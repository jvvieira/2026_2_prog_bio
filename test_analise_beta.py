import unittest

from analise_beta import (
    ROTULOS_GRAFICO_3,
    resumo_beta_por_allele,
    resumo_beta_por_rsid,
    rs_id_de_snp_allele,
)


class RsIdDeSnpAlleleTest(unittest.TestCase):
    def test_extrai_rs_id_do_primeiro_snp(self):
        snp_allele = [
            {"rs_id": "rs71080556", "effect_allele": "G"},
            {"rs_id": "rs12345678", "effect_allele": "A"},
        ]

        self.assertEqual(rs_id_de_snp_allele(snp_allele), "rs71080556")

    def test_sem_snp_allele_devolve_none(self):
        self.assertIsNone(rs_id_de_snp_allele([]))
        self.assertIsNone(rs_id_de_snp_allele(None))


class ResumoBetaPorRsIdTest(unittest.TestCase):
    def test_agrupa_pelo_rs_id_de_snp_allele(self):
        pares = [
            ([{"rs_id": "rs1", "effect_allele": "T"}], 0.2),
            ([{"rs_id": "rs1", "effect_allele": "T"}], -0.1),
            ([{"rs_id": "rs2", "effect_allele": "A"}], 0.3),
        ]

        self.assertEqual(
            resumo_beta_por_rsid(pares),
            [("rs2", 0.3, 1), ("rs1", 0.05, 2)],
        )

    def test_descarta_snp_sem_rs_id(self):
        pares = [([], 0.9), ([{"effect_allele": "A"}], 0.8), ([{"rs_id": "rs1"}], 0.1)]

        self.assertEqual(resumo_beta_por_rsid(pares), [("rs1", 0.1, 1)])


class RotulosGrafico3Test(unittest.TestCase):
    def test_rotulos_referem_se_ao_rs_id(self):
        self.assertEqual(ROTULOS_GRAFICO_3["ylabel"], "rs_id (SNP)")
        self.assertIn("rs_id", ROTULOS_GRAFICO_3["titulo"])


class ResumoBetaPorAlleleTest(unittest.TestCase):
    def test_agrupa_por_allele_e_calcula_media_com_sinal(self):
        pares = [("rs1", 0.2), ("rs1", -0.1), ("rs2", 0.3)]

        self.assertEqual(
            resumo_beta_por_allele(pares),
            [("rs2", 0.3, 1), ("rs1", 0.05, 2)],
        )

    def test_ordena_por_magnitude_absoluta_da_media(self):
        pares = [("menor", 0.1), ("maior", -0.5)]

        self.assertEqual(
            resumo_beta_por_allele(pares),
            [("maior", -0.5, 1), ("menor", 0.1, 1)],
        )

    def test_ignora_allele_ausente_e_beta_nao_parseavel(self):
        pares = [(None, 0.9), ("rs1", None), ("rs2", 0.1)]

        self.assertEqual(resumo_beta_por_allele(pares), [("rs2", 0.1, 1)])

    def test_ignora_beta_nan_do_dataframe(self):
        pares = [("rs1", float("nan")), ("rs2", 0.1)]

        self.assertEqual(resumo_beta_por_allele(pares), [("rs2", 0.1, 1)])

    def test_media_pode_cancelar_quando_sinais_sao_opostos(self):
        pares = [("rs1", 0.2), ("rs1", -0.2)]

        self.assertEqual(resumo_beta_por_allele(pares), [("rs1", 0.0, 2)])


if __name__ == "__main__":
    unittest.main()
