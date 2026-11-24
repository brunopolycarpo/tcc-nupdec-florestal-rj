# -*- coding: utf-8 -*-
"""
Análise de padrão pontual e mapa de calor dos focos de calor no estado do Rio de Janeiro (2016 a 2025).

Procedimentos reproduzidos de "Estatística espacial aplicada à vigilância em saúde"
(RODRIGUES et al., 2025), Módulo 2 (Padrões pontuais), parte prática, convertidos de R para Python:
  1. Teste de Completa Aleatoriedade Espacial (CSR) por contagem em quadrantes, grade nx = 5, ny = 6,
     teste qui-quadrado (equivalente a quadratcount() e quadrat.test() do spatstat).
  2. Função G com envelope de 9 simulações de CSR (equivalente a envelope(fun = Gest, nsim = 9)).
  3. Função K de Ripley com envelope de 9 simulações de CSR (equivalente a envelope(fun = Kest, nsim = 9)).
  4. Estimativa de densidade de Kernel gaussiana com correção de borda de Diggle
     (equivalente a density(ppp, sigma = ..., diggle = TRUE)), com comparação de larguras de banda
     e cálculo da largura de banda "ótima" pela regra MASS::bandwidth.nrd.
A janela de análise é o polígono do estado (equivalente a Window(ppp) <- as.owin(...)).

Entradas: pasta dados/ com focos_2016.csv a focos_2025.csv (BDQueimadas/INPE, AQUA_M-T, RJ);
  ibge/RJ_UF_2025.shp (limite estadual, IBGE MMD 2025); pi_final.pkl (UCs de proteção integral, Apêndice A).
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import box
from shapely.validation import make_valid
from shapely import contains_xy
from scipy.spatial import cKDTree
from scipy.stats import chi2, iqr
from scipy.ndimage import gaussian_filter
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pyproj import Transformer

CRS_M = 31983                      # SIRGAS 2000 / UTM 23S (distâncias em metros)
SIGMA = 4000                       # largura de banda adotada (m)
SIGMAS_COMPARACAO = [2000, 4000, 10000]   # comparação de larguras (como sigma = 100, 200 e 500 no livro)
NSIM = 9                           # número de simulações dos envelopes (como no livro)
rng = np.random.default_rng(2025)
ink, ink2 = '#1f1f1f', '#5c5c5c'
cmap = LinearSegmentedColormap.from_list('laranja', ['#fdf1ea','#f9c9ad','#f29a6c','#eb6834','#c44a1b','#8f3210','#5a1d07'])
saida = []
def log(s=''): print(s); saida.append(s)
TR = Transformer.from_crs(4326, CRS_M, always_xy=True)
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

# ------------------------------------------------------------------ dados e janela
import glob
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('dados/focos_*.csv'))], ignore_index=True)
ucpi = gpd.GeoSeries(pd.read_pickle('pi_final.pkl').geom_rj.values, crs=CRS_M)   # contornos das UCs de proteção integral
W = make_valid(gpd.read_file('ibge/RJ_UF_2025.shp').to_crs(CRS_M).geometry.iloc[0])   # janela de análise: limite do estado (IBGE, MMD 2025)
pts = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326).to_crs(CRS_M)
fora = ~contains_xy(W, pts.geometry.x.values, pts.geometry.y.values)
log(f'Focos fora do limite terrestre (excluídos): {int(fora.sum())}')
pts = pts[~fora].reset_index(drop=True)
x, y = pts.geometry.x.values, pts.geometry.y.values
xy = np.c_[x, y]; n = len(xy); area = W.area; lam = n/area
minx, miny, maxx, maxy = W.bounds
log(f'Focos: {n} | área da janela: {area/1e6:,.0f} km² | intensidade média: {lam*1e8:.2f} focos/100 km²')

def csr(k=None):
    """Uma realização de CSR (Poisson homogêneo) na janela W com a mesma intensidade dos dados."""
    k = rng.poisson(lam*area) if k is None else k
    out = np.empty((0, 2))
    while len(out) < k:
        c = np.c_[rng.uniform(minx, maxx, 3*k+10), rng.uniform(miny, maxy, 3*k+10)]
        out = np.vstack([out, c[contains_xy(W, c[:, 0], c[:, 1])]])
    return out[:k]

# ------------------------------------------------------------------ 1. quadrantes
log('\n1. TESTE DE CSR POR CONTAGEM EM QUADRANTES (grade nx = 5, ny = 6)')
nx, ny = 5, 6
xs, ys = np.linspace(minx, maxx, nx+1), np.linspace(miny, maxy, ny+1)
tiles, a, obs = [], [], []
for i in range(nx):
    for j in range(ny):
        t = box(xs[i], ys[j], xs[i+1], ys[j+1]).intersection(W)
        if t.area <= 0: continue
        tiles.append(t); a.append(t.area)
        obs.append(int(contains_xy(t, x, y).sum()))
a, obs = np.array(a), np.array(obs, float)
esp = n*a/a.sum()                                   # esperado sob CSR, proporcional à área do quadrante na janela
X2 = float(((obs-esp)**2/esp).sum()); gl = len(obs)-1; p = 2*min(chi2.sf(X2, gl), chi2.cdf(X2, gl))  # bilateral
log(f'Quadrantes com área na janela: {len(obs)}')
log(f'X2 = {X2:.2f} | gl = {gl} | p-valor = {p:.2e} (hipótese alternativa bilateral)')
log('Como p < 0,05, rejeita-se a hipótese nula de completa aleatoriedade espacial (CSR).' if p < 0.05
    else 'Como p >= 0,05, não se rejeita a hipótese nula de CSR.')
fig, ax = plt.subplots(figsize=(8, 6))
gpd.GeoSeries(tiles).boundary.plot(ax=ax, color='#b8b8b8', lw=0.6)
gpd.GeoSeries([W]).boundary.plot(ax=ax, color=ink2, lw=0.8)
ax.scatter(x, y, s=4, color='#c44a1b', alpha=0.6, linewidths=0)
for t, o in zip(tiles, obs):
    c = t.representative_point(); ax.text(c.x, c.y, int(o), ha='center', va='center', fontsize=9, color=ink, weight='bold')
ax.set_title(f'Contagem em quadrantes (nx = 5, ny = 6)\nX² = {X2:.1f}; gl = {gl}; p = {p:.1e}', fontsize=11, color=ink)
moldura(ax, escala=50); plt.tight_layout(); plt.savefig('fig_1_quadrantes.png', dpi=200, facecolor='white'); plt.close()

# ------------------------------------------------------------------ 2 e 3. funções G e K com envelope
# Correções de borda (as mesmas oferecidas pelo spatstat nas funções Gest e Kest):
#   Função G: amostra reduzida ("rs"), usando só os focos a mais de r do limite da janela;
#   Função K: translação ("translate"), com peso |W| / |W ∩ (W + (xj - xi))| para cada par.
from scipy.signal import fftconvolve
_rs = 500.0
_xe = np.arange(minx, maxx+_rs, _rs); _ye = np.arange(miny, maxy+_rs, _rs)
_XX, _YY = np.meshgrid((_xe[:-1]+_xe[1:])/2, (_ye[:-1]+_ye[1:])/2)
_mask = contains_xy(W, _XX.ravel(), _YY.ravel()).reshape(_XX.shape).astype(float)
_gW = fftconvolve(_mask, _mask[::-1, ::-1], mode='full')*_rs*_rs     # covariância geométrica |W ∩ W_h|
_cy, _cx = _mask.shape[0]-1, _mask.shape[1]-1
_bnd = W.boundary
def dist_borda(P):
    return shapely.distance(_bnd, shapely.points(P))

def G(P, r, b=None):
    b = dist_borda(P) if b is None else b
    d = cKDTree(P).query(P, k=2)[0][:, 1]
    out = []
    for ri in r:
        ok = b >= ri
        out.append(((d <= ri) & ok).sum()/ok.sum() if ok.sum() else np.nan)
    return np.array(out)

def K(P, r):
    m = len(P); pares = cKDTree(P).query_pairs(r.max(), output_type='ndarray')
    v = P[pares[:, 1]]-P[pares[:, 0]]; dij = np.hypot(v[:, 0], v[:, 1])
    iy_ = np.clip(np.round(v[:, 1]/_rs).astype(int)+_cy, 0, _gW.shape[0]-1)
    ix_ = np.clip(np.round(v[:, 0]/_rs).astype(int)+_cx, 0, _gW.shape[1]-1)
    e = area/np.maximum(_gW[iy_, ix_], 1.0)                              # peso de translação
    return np.array([2*e[dij <= ri].sum() for ri in r])*area/(m*(m-1))

rG, rK = np.linspace(0, 5000, 51), np.linspace(0, 60000, 61)
sims = [csr() for _ in range(NSIM)]
g_obs, k_obs = G(xy, rG), K(xy, rK)
g_sim = np.array([G(s, rG) for s in sims]); k_sim = np.array([K(s, rK) for s in sims])
g_teo = 1-np.exp(-lam*np.pi*rG**2); k_teo = np.pi*rK**2      # valores teóricos sob CSR (Poisson)
for nome, r, ob, sm, te, arq, esc, un in [
        ('Função G', rG, g_obs, g_sim, g_teo, 'fig_2_funcao_G.png', 1, 'G(r)'),
        ('Função K', rK, k_obs, k_sim, k_teo, 'fig_3_funcao_K.png', 1e6, 'K(r) (km²)')]:
    lo, hi = sm.min(0), sm.max(0)                           # envelope: mínimo e máximo das simulações
    fora = (ob[1:] > hi[1:]).sum()
    log(f'\n{"2" if nome=="Função G" else "3"}. {nome.upper()} | envelope de {NSIM} simulações de CSR')
    log(f'Curva observada acima do envelope em {fora} de {len(r)-1} distâncias avaliadas '
        f'(0 a {r[-1]/1000:.0f} km). Leitura gráfica: curva acima do envelope indica padrão agrupado.')
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.fill_between(r/1000, lo/esc, hi/esc, color='#d9d9d9', label=f'Envelope de {NSIM} simulações de CSR')
    ax.plot(r/1000, te/esc, '--', color=ink2, lw=1.2, label='Teórica (CSR / Poisson)')
    ax.plot(r/1000, ob/esc, color='#c44a1b', lw=2, label='Observada (focos de 2016 a 2025)')
    ax.set_xlabel('Distância r (km)'); ax.set_ylabel(un); ax.legend(frameon=False, fontsize=9)
    ax.set_title(f'{nome}: focos de calor, RJ, 2016 a 2025', fontsize=11, color=ink)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    plt.tight_layout(); plt.savefig(arq, dpi=200, facecolor='white'); plt.close()

# ------------------------------------------------------------------ 4. kernel
nrd = lambda v: 4*1.06*min(v.std(ddof=1), iqr(v)/1.34)*n**(-0.2)     # MASS::bandwidth.nrd
bw_x, bw_y = nrd(x), nrd(y)
log(f'\n4. ESTIMATIVA DE DENSIDADE DE KERNEL')
log(f'Largura de banda "ótima" MASS::bandwidth.nrd: x = {bw_x/1000:.1f} km, y = {bw_y/1000:.1f} km '
    f'(desvio-padrão do kernel = 1/4 desses valores: {bw_x/4000:.1f} e {bw_y/4000:.1f} km)')
res, pad = 500, 5000
xe = np.arange(minx-pad, maxx+pad+res, res); ye = np.arange(miny-pad, maxy+pad+res, res)
xc, yc = (xe[:-1]+xe[1:])/2, (ye[:-1]+ye[1:])/2
XX, YY = np.meshgrid(xc, yc)
inW = contains_xy(W, XX.ravel(), YY.ravel()).reshape(XX.shape).astype(float)
ix = np.clip(np.searchsorted(xe, x)-1, 0, len(xc)-1); iy = np.clip(np.searchsorted(ye, y)-1, 0, len(yc)-1)

def kernel(sx, sy=None):
    """Kernel gaussiano com correção de borda de Diggle (density(..., diggle = TRUE)). Saída: focos/100 km²."""
    sy = sx if sy is None else sy
    s = (sy/res, sx/res)
    e = gaussian_filter(inW, s, mode='constant')                # massa do kernel dentro da janela
    w = 1/np.maximum(e[iy, ix], 1e-6)                           # peso de cada foco: 1/e(x_i)
    cnt = np.zeros_like(inW); np.add.at(cnt, (iy, ix), w)
    d = gaussian_filter(cnt, s, mode='constant')/(res*res)*1e8
    d[inW == 0] = np.nan
    return d

redec = gpd.read_file('redec_rj.gpkg').to_crs(CRS_M)          # REDECs: municípios da MMD 2025 (IBGE) agrupados conforme a SEDEC-RJ
import matplotlib.patheffects as pe
DESL = {'Baixada Fluminense': (-14000, 13000), 'Metropolitana': (14000, -9000)}   # deslocamento dos rótulos (m) para evitar sobreposição
from shapely.ops import unary_union, linemerge
REDEC_INT = gpd.GeoSeries([linemerge(unary_union(list(redec.boundary)).difference(unary_union(list(redec.geometry)).boundary.buffer(100)))], crs=CRS_M)
def contorno(ax, rotulos=False, escala=50, fs=8):
    ucpi.boundary.plot(ax=ax, color='#1b5e57', lw=0.45)
    if rotulos:   # mapa de calor final: divisas entre REDECs em traço escuro espesso sobre faixa branca
        REDEC_INT.plot(ax=ax, color='white', lw=4.2, zorder=5, capstyle='round', joinstyle='round')
        REDEC_INT.plot(ax=ax, color='#1f1f1f', lw=1.9, zorder=5.1, capstyle='round', joinstyle='round')
    else:
        redec.boundary.plot(ax=ax, color='#7a7a7a', lw=0.8)
    gpd.GeoSeries([W]).boundary.plot(ax=ax, color=ink2, lw=0.8)
    if rotulos:
        for _, r in redec.iterrows():
            pt = r.geometry.representative_point(); dx, dy = DESL.get(r.redec, (0, 0))
            ax.text(pt.x + dx, pt.y + dy, r.redec, ha='center', va='center', fontsize=10, color=ink, fontweight='bold',
                    path_effects=[pe.withStroke(linewidth=3.5, foreground='white')], zorder=6)
    ax.set_xlim(minx-pad, maxx+pad); ax.set_ylim(miny-pad, maxy+pad); moldura(ax, fs=fs, escala=escala)

# comparação de larguras de banda (como a figura sigma = 100, 200 e 500 do livro) + largura MASS
paineis = [(f'sigma = {s/1000:.0f} km', s, s) for s in SIGMAS_COMPARACAO] + \
          [(f'MASS::bandwidth.nrd\nsigma x = {bw_x/4000:.1f} km; sigma y = {bw_y/4000:.1f} km', bw_x/4, bw_y/4)]
fig, axs = plt.subplots(2, 2, figsize=(14, 9.5))
for ax, (tit, sx, sy) in zip(axs.ravel(), paineis):
    d = kernel(sx, sy)
    im = ax.imshow(d, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap=cmap, vmin=0)
    ax.scatter(x, y, s=1.5, c=ink, alpha=0.4, linewidths=0); contorno(ax, fs=7)
    ax.set_title(f'{tit}\nmáximo = {np.nanmax(d):.1f} focos/100 km²', fontsize=11, color=ink)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=8, colors=ink2)
fig.suptitle('Comparação de larguras de banda da estimativa de densidade de Kernel (focos de calor, RJ, 2016 a 2025). Cada painel tem escala própria.',
             fontsize=12, color=ink)
plt.tight_layout(); plt.savefig('fig_4_kernel_larguras_de_banda.png', dpi=200, facecolor='white'); plt.close()

# mapa de calor final
d = kernel(SIGMA)
fig, ax = plt.subplots(figsize=(12, 7.6))
im = ax.imshow(d, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap=cmap, vmin=0)
contorno(ax, rotulos=True, escala=None)
ax.scatter(x, y, s=2.5, c=ink, alpha=0.45, linewidths=0, label=f'Focos de calor, 2016 a 2025 (n = {n})')
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, shrink=0.8)
cb.set_label('Densidade de focos (focos/100 km²)', color=ink); cb.outline.set_visible(False); cb.ax.tick_params(colors=ink2)
ax.set_title('Mapa de calor dos focos de calor no estado do Rio de Janeiro, 2016 a 2025', fontsize=14, color=ink, loc='left', pad=20)
ax.text(0, 1.0, f'Estimativa de densidade de Kernel gaussiana, sigma = {SIGMA/1000:.0f} km, com correção de borda de Diggle',
        transform=ax.transAxes, fontsize=9.5, color=ink2, va='bottom')
from matplotlib.lines import Line2D
hs, ls = ax.get_legend_handles_labels()
hs.append(Line2D([0],[0], color='#1b5e57', lw=1.2)); ls.append('Limite de UC de proteção integral')
hs.append(Line2D([0],[0], color='#1f1f1f', lw=1.9, path_effects=[pe.withStroke(linewidth=4.2, foreground='white')])); ls.append('Divisa entre REDECs (SEDEC-RJ)')
ax.legend(hs, ls, loc='lower right', frameon=False, fontsize=9, markerscale=4, labelcolor=ink2)
sx0, sy0 = minx+190000, miny+9000
ax.plot([sx0, sx0+50000], [sy0, sy0], color=ink, lw=2.5, solid_capstyle='butt')
ax.text(sx0+25000, sy0+3500, '50 km', ha='center', fontsize=9, color=ink)
ax.annotate('N', xy=(maxx-8000, maxy-4000), xytext=(maxx-8000, maxy-22000), ha='center', fontsize=11, color=ink,
            arrowprops=dict(arrowstyle='-|>', color=ink, lw=1.5))
ax.text(1, -0.07, 'Fonte dos dados: BDQueimadas/INPE (satélite de referência AQUA_M-T); unidades de conservação: INEA;\nlimite estadual e municípios: IBGE (Malha Municipal Digital 2025); REDECs: SEDEC-RJ. '
        'Projeção: SIRGAS 2000 / UTM 23S.', transform=ax.transAxes, ha='right', va='top', fontsize=8, color=ink2)
plt.savefig('fig_5_mapa_calor_focos_RJ_2016_2025.png', dpi=250, bbox_inches='tight', facecolor='white'); plt.close()
log(f'Mapa de calor final: sigma = {SIGMA/1000:.0f} km | densidade máxima = {np.nanmax(d):.1f} focos/100 km²')
open('resultados.txt', 'w', encoding='utf-8').write('\n'.join(saida))
