# -*- coding: utf-8 -*-
"""Mapa de prioridades (Figura 8) e tabela completa da hierarquização.
Entradas: ucs_resultado.pkl (Apêndice C); ibge/RJ_UF_2025.shp; redec_rj.gpkg (Apêndice A).
Saídas: fig_6_prioridade_ucs.png; hierarquizacao_ucs_protecao_integral.xlsx; número de UCs por classe e REDEC."""
import warnings; warnings.filterwarnings('ignore')
import pandas as pd, geopandas as gpd, shapely, numpy as np
from shapely.validation import make_valid
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pyproj import Transformer
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
ink, ink2 = '#1f1f1f', '#5c5c5c'
TR = Transformer.from_crs(4326, 31983, always_xy=True)
def moldura(ax, fs=8, escala=None):
    """Coordenadas geográficas na moldura (meridianos e paralelos de 1 em 1 grau) e, se pedida, escala gráfica e norte."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); xt, xl, yt, yl = [], [], [], []
    t = np.linspace(-26, -18, 400)
    for lon in range(-46, -39):
        x, y = TR.transform(np.full_like(t, lon), t); ax.plot(x, y, color='#d0d0d0', lw=0.4, zorder=0.5)
        xb = np.interp(y0, y, x)
        if x0 < xb < x1: xt.append(xb); xl.append(f'{-lon}°O')
    t = np.linspace(-48, -38, 400)
    for lat in range(-24, -19):
        x, y = TR.transform(t, np.full_like(t, lat)); ax.plot(x, y, color='#d0d0d0', lw=0.4, zorder=0.5)
        yb = np.interp(x0, x, y)
        if y0 < yb < y1: yt.append(yb); yl.append(f'{-lat}°S')
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_axis_on(); ax.set_xlabel(''); ax.set_ylabel('')
    ax.set_xticks(xt, xl, fontsize=fs, color=ink2); ax.set_yticks(yt, yl, fontsize=fs, color=ink2, rotation=90, va='center')
    ax.tick_params(length=3, color=ink2)
    for s in ax.spines.values(): s.set_color(ink2); s.set_linewidth(0.6)
    if escala:
        sx, sy = x0 + 0.04*(x1-x0), y0 + 0.05*(y1-y0)
        ax.plot([sx, sx+escala*1000], [sy, sy], color=ink, lw=2.5, solid_capstyle='butt', zorder=7)
        ax.text(sx+escala*500, sy+0.015*(y1-y0), f'{escala} km', ha='center', fontsize=fs, color=ink, zorder=7)
        ax.annotate('N', xy=(x1-0.04*(x1-x0), y1-0.04*(y1-y0)), xytext=(x1-0.04*(x1-x0), y1-0.14*(y1-y0)), ha='center',
                    fontsize=fs+2, color=ink, arrowprops=dict(arrowstyle='-|>', color=ink, lw=1.2), zorder=7)
cores = {'Alta':'#8f3210','Média':'#e8743b','Baixa':'#f6c3a3'}
g = pd.read_pickle('ucs_resultado.pkl'); g = gpd.GeoDataFrame(g, geometry='geometry', crs=31983)
W = make_valid(gpd.read_file('ibge/RJ_UF_2025.shp').to_crs(31983).geometry.iloc[0]); minx, miny, maxx, maxy = W.bounds  # limite do estado, IBGE MMD 2025
g = g.sort_values('taxa_eb5', ascending=False).reset_index(drop=True); g['posicao'] = np.arange(1, len(g)+1)
fig, ax = plt.subplots(figsize=(12, 7.6))
gpd.GeoSeries([W]).plot(ax=ax, color='#f4f3f0', edgecolor=ink2, lw=0.8)
import matplotlib.patheffects as pe
DESL = {'Baixada Fluminense': (-14000, 13000), 'Metropolitana': (14000, -9000)}   # deslocamento dos rótulos (m) para evitar sobreposição
redec = gpd.read_file('redec_rj.gpkg').to_crs(31983)          # REDECs: municípios da MMD 2025 (IBGE) agrupados conforme a SEDEC-RJ
redec.boundary.plot(ax=ax, color='#7a7a7a', lw=0.8, zorder=2)
for _, rr in redec.iterrows():
    pt = rr.geometry.representative_point(); dx, dy = DESL.get(rr.redec, (0, 0))
    ax.text(pt.x + dx, pt.y + dy, rr.redec, ha='center', va='center', fontsize=8, color=ink2, style='italic',
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white')], zorder=5)
for cl in ['Baixa','Média','Alta']:
    s = g[g.prioridade5==cl]
    s.plot(ax=ax, color=cores[cl], edgecolor='white', lw=0.2)
    rp = s.geometry.representative_point()
    ax.scatter(rp.x, rp.y, s={'Baixa':10,'Média':22,'Alta':60}[cl], color=cores[cl], edgecolor='white', linewidth=0.6, zorder=3)
alta = g[g.prioridade5=='Alta']
for _, r in alta.iterrows():
    pt = r.geometry.representative_point()
    off = {9:(-14,-12), 5:(5,3), 1:(-12,2)}.get(int(r.posicao),(4,4))
    ax.annotate(str(r.posicao), (pt.x, pt.y), xytext=off, textcoords='offset points', fontsize=8.5, color=ink, weight='bold', zorder=4)
leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor=cores[c], markeredgecolor='white', markersize=m, label=f'{c} ({(g.prioridade5==c).sum()} UCs)')
       for c, m in [('Alta',10),('Média',8),('Baixa',6)]] + [Line2D([0],[0], color='#7a7a7a', lw=1.0, label='Limite de REDEC (SEDEC-RJ)')]
ax.legend(handles=leg, title='Prioridade', loc='lower right', frameon=False, fontsize=9, title_fontsize=9.5, labelcolor=ink2)
ax.set_xlim(minx-5000, maxx+5000); ax.set_ylim(miny-5000, maxy+5000); moldura(ax)
ax.set_title('Prioridade das unidades de conservação de proteção integral para a implantação de NUPDECs Florestais', fontsize=13, color=ink, loc='left', pad=20)
ax.text(0, 1.0, 'Taxa de focos de calor (2016 a 2025) na UC e no entorno de 5 km, suavizada pelo método Bayesiano empírico; classes por quebras de Jenks. '
        'Números: posição das UCs de prioridade alta.', transform=ax.transAxes, fontsize=9, color=ink2, va='bottom')
sx0, sy0 = minx+190000, miny+9000
ax.plot([sx0, sx0+50000], [sy0, sy0], color=ink, lw=2.5, solid_capstyle='butt'); ax.text(sx0+25000, sy0+3500, '50 km', ha='center', fontsize=9, color=ink)
ax.annotate('N', xy=(maxx-8000, maxy-4000), xytext=(maxx-8000, maxy-22000), ha='center', fontsize=11, color=ink, arrowprops=dict(arrowstyle='-|>', color=ink, lw=1.5))
ax.text(1, -0.07, 'Fonte dos dados: BDQueimadas/INPE (AQUA_M-T); unidades de conservação: INEA;\nlimite estadual e municípios: IBGE (Malha Municipal Digital 2025); REDECs: SEDEC-RJ. Projeção: SIRGAS 2000 / UTM 23S.',
        transform=ax.transAxes, ha='right', va='top', fontsize=8, color=ink2)
plt.savefig('fig_6_prioridade_ucs.png', dpi=250, bbox_inches='tight', facecolor='white'); plt.close()

# tabela Excel
t = g[['posicao','nome','esfera','categoria','area_uc_km2','focos_dentro_uc','area_zona5_km2','focos_zona5','taxa_bruta5','taxa_eb5','prioridade5','taxa_eb10','prioridade10','reatribuida']].copy()
t.columns = ['Posição','Unidade de conservação','Esfera','Categoria','Área da UC (km²)','Focos dentro da UC','Área da zona 5 km (km²)',
             'Focos na zona 5 km','Taxa bruta 5 km (focos/100 km²/ano)','Taxa suavizada 5 km (focos/100 km²/ano)','Prioridade 5 km',
             'Taxa suavizada 10 km (focos/100 km²/ano)','Prioridade 10 km','Área reatribuída pelo critério da menor UC']
t['Unidade de conservação'] = t['Unidade de conservação'].str.title().str.replace(' Da ',' da ').str.replace(' De ',' de ').str.replace(' Do ',' do ').str.replace(' Das ',' das ').str.replace(' Dos ',' dos ').str.replace(' E ',' e ')
t['Área reatribuída pelo critério da menor UC'] = t['Área reatribuída pelo critério da menor UC'].map({True:'Sim', False:'Não'})
for c in t.columns[4:]:
    if t[c].dtype.kind=='f': t[c]=t[c].round(2)
with pd.ExcelWriter('hierarquizacao_ucs_protecao_integral.xlsx', engine='openpyxl') as xw:
    t.to_excel(xw, index=False, sheet_name='Hierarquização')
    ws = xw.sheets['Hierarquização']
    from openpyxl.styles import Font, Alignment, PatternFill
    for cell in ws[1]: cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical='center'); cell.fill = PatternFill('solid', fgColor='E7E6E6')
    widths = [8,58,11,24,11,11,12,11,14,14,11,14,11,14]
    for i, w in enumerate(widths): ws.column_dimensions[chr(65+i)].width = w
    ws.row_dimensions[1].height = 60; ws.freeze_panes = 'C2'
print(t.groupby(['Prioridade 5 km','Esfera']).size().unstack(fill_value=0).reindex(['Alta','Média','Baixa']).to_string())

# número de unidades por classe de prioridade em cada REDEC (cada UC associada à REDEC que contém a maior parte da sua área)
ov = gpd.overlay(gpd.GeoDataFrame(g[['geometry']].reset_index(), crs=31983), redec, how='intersection', keep_geom_type=True)
ov['a'] = ov.area
g['redec'] = ov.sort_values('a').groupby('index').tail(1).set_index('index').redec
tr = pd.crosstab(g.redec, g.prioridade5).reindex(columns=['Alta', 'Média', 'Baixa']).fillna(0).astype(int)
tr['Total'] = tr.sum(axis=1)
print(tr.sort_values(['Alta', 'Média'], ascending=False).to_string())
