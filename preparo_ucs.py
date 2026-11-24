# -*- coding: utf-8 -*-
"""Preparo das UCs: leitura das camadas do INEA, correção de geometrias, classificação pelas categorias de
proteção integral do SNUC (Lei 9.985/2000, art. 8º), exclusões e recorte pelo limite terrestre do estado (IBGE, MMD 2025).
Também delimita as REDECs (SEDEC-RJ) pela união dos municípios que as compõem e mede o porte das APAs
e a sua sobreposição com as unidades de proteção integral (seção 1.2).
Entradas: ucs_shp/ (camadas do INEA); ibge/RJ_UF_2025.shp; ibge_mun/RJ_Municipios_2025.shp; redec_municipios.csv
(relação de municípios por REDEC, conforme a SEDEC-RJ). Saídas: pi_final.pkl; redec_rj.gpkg."""
import warnings; warnings.filterwarnings('ignore')
import unicodedata, pandas as pd, geopandas as gpd, shapely
from shapely.validation import make_valid
CRS = 31983
W = make_valid(gpd.read_file('ibge/RJ_UF_2025.shp').to_crs(CRS).geometry.iloc[0])
PI_SNUC = {'ESEC','EEE','REBIO','PARQUE','PARNA','MONA','REVIS'}
def carregar(pasta, arq, esfera, campo_cat, campo_uf=None):
    g = gpd.read_file(f'ucs_shp/{pasta}/{arq}.shp').to_crs(CRS)
    invalidas = int((~g.is_valid).sum())
    g['geometry'] = g.geometry.apply(make_valid).buffer(0)
    out = gpd.GeoDataFrame({'nome': g['nome'], 'esfera': esfera, 'cat': g[campo_cat].astype(str),
                            'uf': g[campo_uf].astype(str) if campo_uf else 'RJ'}, geometry=g.geometry, crs=CRS)
    print(f'{esfera}: {len(g)} feições, {invalidas} geometrias inválidas corrigidas')
    return out
uc = pd.concat([carregar('est','gpl_ucs_estaduais_limites_inearj_me_2023','Estadual','categoria'),
                carregar('fed','gpl_ucs_federais_novembro_2020','Federal','sigla','uf'),
                carregar('mun','gpl_ucs_municipais_2024_icms_2025','Municipal','categoria')], ignore_index=True)
uc.loc[(uc.esfera=='Municipal') & (uc.cat=='PI'), 'cat'] = 'PARQUE'      # PNM Atalaia: categoria registrada como "PI"
pi = uc[uc.cat.isin(PI_SNUC)].copy()                                      # exclui uso sustentável e a RESEC Juatinga
norm = lambda s: unicodedata.normalize('NFKD', str(s)).encode('ascii','ignore').decode().upper()
pi = pi[~pi.nome.map(norm).str.contains(r'\bMARINHA\b')]
pi['area_total_km2'] = pi.area/1e6
pi['geom_rj'] = [make_valid(g.intersection(W)) for g in pi.geometry]
pi['area_terra_rj_km2'] = [g.area/1e6 for g in pi.geom_rj]
pi = pi[pi.area_terra_rj_km2 > 0].reset_index(drop=True)
print(f'UCs de proteção integral: {len(pi)}'); print(pd.crosstab(pi.cat, pi.esfera, margins=True).to_string())
u = shapely.union_all(pi.geom_rj.values)
print(f'Área terrestre coberta: {u.area/1e6:,.0f} km² = {u.area/W.area*100:.1f}% do estado')
pd.DataFrame(pi.drop(columns='geometry')).to_pickle('pi_final.pkl')

# REDECs: união dos municípios da Malha Municipal Digital 2025 (IBGE) segundo a relação da SEDEC-RJ
mun = gpd.read_file('ibge_mun/RJ_Municipios_2025.shp').to_crs(CRS)
rel = pd.read_csv('redec_municipios.csv')
mun['k'] = mun.NM_MUN.map(norm); rel['k'] = rel.municipio.map(norm)
mun = mun.merge(rel[['k', 'redec']], on='k', how='left')
assert mun.redec.notna().all(), 'município sem REDEC'
mun['geometry'] = mun.geometry.apply(make_valid)
mun.dissolve('redec').reset_index()[['redec', 'geometry']].to_file('redec_rj.gpkg', driver='GPKG')
print(f'REDECs: {mun.redec.nunique()} regionais, {len(mun)} municípios')

# Porte das APAs e sobreposição com as unidades de proteção integral (seção 1.2)
us = uc[~uc.cat.isin(PI_SNUC)].copy()
us['geom_rj'] = [make_valid(g.intersection(W)) for g in us.geometry]
us['area_km2'] = [g.area/1e6 for g in us.geom_rj]
apa = us[(us.cat == 'APA') & (us.area_km2 > 0.001)]
print(f'APAs no estado: {len(apa)}; mediana {apa.area_km2.median():.1f} km²; {(apa.area_km2 >= 100).sum()} com 100 km² ou mais')
U = make_valid(shapely.union_all(apa.geom_rj.values))
f = [g.intersection(U).area/g.area for g in pi.geom_rj]
print(f'UCs de proteção integral com metade ou mais da área em APA: {sum(x >= 0.5 for x in f)} de {len(pi)}')
grandes = apa[apa.area_km2 >= 100]
n = sum(g.intersection(u).area/1e6 > 0.01 for g in grandes.geom_rj)
print(f'APAs com 100 km² ou mais que contêm UC de proteção integral: {n} de {len(grandes)}')
