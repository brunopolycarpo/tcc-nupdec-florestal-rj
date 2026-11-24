# -*- coding: utf-8 -*-
"""
Entradas: dados/focos_*.csv; ibge/RJ_UF_2025.shp; pi_final.pkl (Apêndice A). Saídas: ucs_resultado.pkl; hierarquizacao_ucs.csv.
Hierarquização das UCs de proteção integral do RJ segundo a incidência de focos de calor (2016 a 2025).
Etapas: (1) critério da menor UC para sobreposições; (2) zona = UC + entorno (5 km e 10 km), só em terra;
(3) taxa de focos por 100 km² por ano; (4) suavização Bayesiana empírica global (esda);
(5) classes de prioridade por quebras de Jenks (mapclassify.FisherJenks, k = 3); (6) sensibilidade 5 km x 10 km.
"""
import warnings; warnings.filterwarnings('ignore')
import glob, numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.validation import make_valid
from shapely import contains_xy
from esda.smoothing import Empirical_Bayes
import mapclassify

CRS = 31983; T = 10
W = make_valid(gpd.read_file('ibge/RJ_UF_2025.shp').to_crs(CRS).geometry.iloc[0])   # limite do estado, IBGE MMD 2025
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('dados/focos_*.csv'))], ignore_index=True)
p = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326).to_crs(CRS)
p = p[contains_xy(W, p.geometry.x, p.geometry.y)].reset_index(drop=True)
X, Y = p.geometry.x.values, p.geometry.y.values
print(f'Focos analisados: {len(p)}')

pi = pd.read_pickle('pi_final.pkl').reset_index(drop=True)
cat_map = {'EEE':'Estação ecológica','ESEC':'Estação ecológica','REBIO':'Reserva biológica','PARQUE':'Parque',
           'PARNA':'Parque','MONA':'Monumento natural','REVIS':'Refúgio de vida silvestre'}
g = gpd.GeoDataFrame({'nome': pi.nome, 'esfera': pi.esfera, 'categoria': pi.cat.map(cat_map)},
                     geometry=gpd.GeoSeries(pi.geom_rj.values, crs=CRS), crs=CRS)
g['area_original_km2'] = g.area/1e6

# (1) critério da menor UC: o trecho sobreposto fica com a UC de menor área
acc = None; nova = {}
for i in g.sort_values('area_original_km2').index:
    geom = g.geometry[i]
    nova[i] = geom if acc is None else make_valid(geom.difference(acc))
    acc = geom if acc is None else make_valid(shapely.union(acc, geom))
g['geometry'] = gpd.GeoSeries(pd.Series(nova), crs=CRS)
g['area_uc_km2'] = g.area/1e6
g['reatribuida'] = (g.area_original_km2 - g.area_uc_km2) > 0.01
print(f'UCs com parte da área reatribuída: {int(g.reatribuida.sum())} | UCs: {len(g)}')

# focos dentro das UCs (sem dupla contagem, pois os polígonos já não se sobrepõem)
allu = shapely.union_all(g.geometry.values)
g['focos_dentro_uc'] = [int(contains_xy(geom, X, Y).sum()) for geom in g.geometry]

def hierarquizar(dist_m):
    zonas = [make_valid(geom.buffer(dist_m).intersection(W)) for geom in g.geometry]
    a = np.array([z.area/1e6 for z in zonas])
    n = np.array([int(contains_xy(z, X, Y).sum()) for z in zonas])
    taxa = n/(a*T)*100
    eb = Empirical_Bayes(n.astype(float), a*T)          # exposição = área x anos
    taxa_eb = np.asarray(eb.r).ravel()*100
    fj = mapclassify.FisherJenks(taxa_eb, k=3)
    classe = np.array(['Baixa','Média','Alta'])[fj.yb]
    return zonas, a, n, taxa, taxa_eb, classe, fj.bins

res = {}
for d in (5000, 10000):
    res[d] = hierarquizar(d)
    z, a, n, t, te, c, bins = res[d]
    print(f'\nEntorno {d//1000} km | taxa global = {n.sum()/(a.sum()*T)*100:.2f} focos/100 km²/ano (zonas sobrepostas contam mais de uma vez)')
    print(f'  limites superiores das classes (Jenks): {np.round(bins,2).tolist()}')
    print('  UCs por classe:', pd.Series(c).value_counts().reindex(['Alta','Média','Baixa']).to_dict())

z5, a5, n5, t5, te5, c5, b5 = res[5000]
z10, a10, n10, t10, te10, c10, b10 = res[10000]
g['area_zona5_km2']=a5; g['focos_zona5']=n5; g['taxa_bruta5']=t5; g['taxa_eb5']=te5; g['prioridade5']=c5
g['area_zona10_km2']=a10; g['focos_zona10']=n10; g['taxa_bruta10']=t10; g['taxa_eb10']=te10; g['prioridade10']=c10
alta5, alta10 = set(g.index[g.prioridade5=='Alta']), set(g.index[g.prioridade10=='Alta'])
print(f'\nSENSIBILIDADE | prioridade alta: {len(alta5)} UCs com 5 km, {len(alta10)} com 10 km, {len(alta5 & alta10)} em comum '
      f'({len(alta5 & alta10)/len(alta5)*100:.0f}% das de 5 km mantidas)')
mud = pd.crosstab(g.prioridade5, g.prioridade10).reindex(index=['Alta','Média','Baixa'], columns=['Alta','Média','Baixa']).fillna(0).astype(int)
print('Tabela cruzada (linhas: 5 km; colunas: 10 km):'); print(mud.to_string())

# focos dentro das UCs x entorno (descritivo)
anel = make_valid(shapely.union_all(z5).difference(allu))
nd, na = int(contains_xy(allu, X, Y).sum()), int(contains_xy(anel, X, Y).sum())
print(f'\nFocos dentro das UCs: {nd} ({nd/len(p)*100:.1f}%) | densidade {nd/(allu.area/1e6)/T*100:.2f} focos/100 km²/ano')
print(f'Focos no entorno de 0 a 5 km (fora das UCs): {na} ({na/len(p)*100:.1f}%) | densidade {na/(anel.area/1e6)/T*100:.2f}')
print(f'Estado: densidade {len(p)/(W.area/1e6)/T*100:.2f}')

g.to_pickle('ucs_resultado.pkl')
cols = ['nome','esfera','categoria','area_original_km2','area_uc_km2','focos_dentro_uc','area_zona5_km2','focos_zona5',
        'taxa_bruta5','taxa_eb5','prioridade5','focos_zona10','taxa_eb10','prioridade10']
out = g.sort_values('taxa_eb5', ascending=False)[cols].round(3)
out.insert(0, 'posicao', range(1, len(out)+1))
out.to_csv('hierarquizacao_ucs.csv', index=False, encoding='utf-8-sig')
print('\nTOP 20 (entorno de 5 km, taxa suavizada):')
print(out.head(20)[['posicao','nome','esfera','categoria','area_uc_km2','focos_zona5','taxa_bruta5','taxa_eb5','prioridade5','prioridade10']].to_string(index=False))
