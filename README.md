# NUPDECs Florestais no Estado do Rio de Janeiro: priorização por análise estatística espacial de focos de calor

Código e dados do Trabalho de Conclusão de Curso **"Implantação de Núcleos Comunitários de Proteção e Defesa Civil Florestais no Estado do Rio de Janeiro: priorização por meio da análise estatística espacial de focos de calor"**, de Bruno Polycarpo Palmerim Dias, apresentado à Universidade Federal Rural do Rio de Janeiro (UFRRJ), sob orientação do Prof. Manoel Pereira, em 2026.

Os quatro scripts são os mesmos publicados nos Apêndices A a D do trabalho. Os métodos seguem Rodrigues et al. (2025), *Estatística espacial aplicada à vigilância em saúde*, com os procedimentos do pacote R spatstat reproduzidos em Python.

## Ordem de execução

Os scripts devem ser executados nesta ordem, a partir da pasta raiz do repositório, porque cada um usa arquivos gerados pelo anterior.

| Ordem | Script | Apêndice | O que faz | Gera |
|---|---|---|---|---|
| 1 | `preparo_ucs.py` | A | Lê as camadas de UCs do INEA, corrige geometrias, seleciona as UCs de proteção integral e recorta pelo limite do estado; delimita as REDECs; mede o porte das APAs e a sobreposição com as UCs de proteção integral | `pi_final.pkl`, `redec_rj.gpkg` |
| 2 | `analise_padroes_pontuais.py` | B | Teste de aleatoriedade espacial por quadrantes, funções G e K com envelopes e mapa de calor (Kernel) | `fig_1` a `fig_5` |
| 3 | `hierarquizacao_ucs.py` | C | Critério da menor UC, zonas de 5 km e 10 km, taxas brutas e suavizadas (Bayes empírico global), classes de prioridade (Jenks) e sensibilidade | `ucs_resultado.pkl`, `hierarquizacao_ucs.csv` |
| 4 | `mapa_prioridades.py` | D | Mapa de prioridades, planilha com as 229 UCs e contagem por classe e REDEC | `fig_6`, `hierarquizacao_ucs_protecao_integral.xlsx` |

```bash
pip install -r requirements.txt
python preparo_ucs.py
python analise_padroes_pontuais.py
python hierarquizacao_ucs.py
python mapa_prioridades.py
```

Testado com Python 3.13 e as versões listadas em `requirements.txt`. A execução completa leva alguns minutos.

## Resultados esperados

Com os dados deste repositório, os scripts reproduzem os números do trabalho, entre eles:

- 6.652 focos de calor analisados (2016 a 2025);
- 229 UCs de proteção integral, com 4.764 km² (10,9% do estado);
- 287 focos no interior das UCs (4,3%) e 3.076 no entorno de 0 a 5 km (46,2%);
- 12 UCs de prioridade alta, 54 de prioridade média e 163 de prioridade baixa (entorno de 5 km).

## Dados

| Pasta ou arquivo | Conteúdo | Fonte | Obtido em |
|---|---|---|---|
| `dados/` | Focos de calor do satélite de referência AQUA_M-T, estado do Rio de Janeiro, 2016 a 2025, um arquivo por ano | INPE, Programa Queimadas, BDQueimadas | 15/09/2026 |
| `ucs_shp/est/` | Limites das UCs estaduais (2023) | INEA, Geoportal | 09/10/2026 |
| `ucs_shp/fed/` | Limites das UCs federais (novembro de 2020) | INEA, Geoportal | 09/10/2026 |
| `ucs_shp/mun/` | Limites das UCs municipais (2024, ICMS Ecológico 2025) | INEA, Geoportal | 09/10/2026 |
| `ibge/` | Limite do estado do Rio de Janeiro | IBGE, Malha Municipal Digital 2025 | 14/10/2026 |
| `ibge_mun/` | Limites dos 92 municípios | IBGE, Malha Municipal Digital 2025 | 14/10/2026 |
| `redec_municipios.csv` | Municípios de cada Coordenadoria Regional de Defesa Civil (REDEC) | SEDEC-RJ, site oficial | 14/10/2026 |

Observações sobre os dados:

- As coordenadas dos focos de 2016 a 2020 e de parte de 2021 têm três casas decimais (precisão da ordem de 100 m); as demais têm cinco.
- Um foco de calor indica uma anomalia térmica em um pixel; não equivale a um incêndio nem mede área queimada.
- A camada federal é a mais antiga das três (2020).
- Os dados pertencem às respectivas fontes e seguem as condições de uso de cada uma. Para trabalhos novos, recomenda-se obtê-los diretamente nas fontes.

## Referências principais

INSTITUTO NACIONAL DE PESQUISAS ESPACIAIS (INPE). **BDQueimadas**. São José dos Campos: INPE. Disponível em: https://terrabrasilis.dpi.inpe.br/queimadas/bdqueimadas/.

INSTITUTO ESTADUAL DO AMBIENTE (INEA). **Geoportal INEA**. Rio de Janeiro: INEA. Disponível em: https://geoportal.inea.rj.gov.br/.

INSTITUTO BRASILEIRO DE GEOGRAFIA E ESTATÍSTICA (IBGE). **Malha Municipal Digital do Brasil: versão 2025**. Rio de Janeiro: IBGE, 2025.

RODRIGUES, D. O.; TASSINARI, W. S.; MORAIS, J. H. A.; FREITAS, L. P.; CRUZ, O. G. **Estatística espacial aplicada à vigilância em saúde**. Florianópolis: UFSC; Instituto Todos pela Saúde; ABRASCO, 2025.

## Uso de inteligência artificial

O código foi desenvolvido com apoio da ferramenta de inteligência artificial Claude (Anthropic), conforme declarado na seção 3.3.3 do trabalho. Os procedimentos, os resultados e as interpretações foram verificados pelo autor.

## Licença

O código está sob a licença MIT (arquivo `LICENSE`). A licença não se aplica aos dados, que seguem as condições das fontes.
