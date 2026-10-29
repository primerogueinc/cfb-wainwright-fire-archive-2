import json,os,numpy as np,pandas as pd,geopandas as gpd,rasterio
from rasterio.windows import from_bounds
from scipy.stats import norm
from pyproj import Transformer
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LightSource
import matplotlib.patches as mpatches
plt.rcParams['font.family']='DejaVu Sans'
os.environ['CURL_CA_BUNDLE']='/root/.ccr/ca-bundle.crt'
tr=Transformer.from_crs(4326,32612,always_xy=True); px,py=tr.transform(-110.9005,52.7695)
B=np.load('best.npz'); sh=tuple(B['shape']); rx,ry=B['rx'],B['ry']
ext=(rx.min()-30+2520,rx.max()+30,ry.min()-30,ry.max()+30)
from scipy.ndimage import median_filter
CROP=int(2520/60)
def cls(L,M=140):
    L=median_filter(L.reshape(sh),3)[:,CROP:]
    p=norm.cdf((M-L)/8.0); c=np.zeros_like(p); c[p>=.1]=1; c[p>=.5]=2; c[p>=.9]=3; return c
cmap=ListedColormap(['#b2182b','#f4a582','#92c5de','#2166ac'])
labels=['Unlikely (<10%)','Marginal (10–50%)','Probable (50–90%)','Likely (≥90%)']
# basemap S2 Oct 19
r=json.load(open('../v6/s2_lens/search.json')); f=[f for f in r['features'] if f['id']=='S2B_12UWD_20261019_0_L2A'][0]
with rasterio.open(f['assets']['visual']['href']) as s:
    w=from_bounds(ext[0],ext[2],ext[1],ext[3],s.transform); rgb=s.read(window=w,boundless=True).transpose(1,2,0)
burn=gpd.read_file('/home/claude/cfb-wainwright-fire-archive-2/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson').to_crs(32612)
t=pd.read_csv('towers.csv'); t['x'],t['y']=tr.transform(t.lon.values,t.lat.values)
wg=pd.read_csv('../radio/wigle_points.csv'); wg=wg[wg.kind=='cell']; wg['x'],wg['y']=tr.transform(wg.lon.values,wg.lat.values)
pts={'Pairs range\n(best guess)':(-110.9005,52.7695),'Range 2':(-110.90089,52.772)}
def deco(ax,legend=True,small=False):
    burn.boundary.plot(ax=ax,color='white',lw=1.4 if not small else 1); burn.boundary.plot(ax=ax,color='#7a1020',lw=.6)
    tt=t[(t.x>ext[0])&(t.x<ext[1])&(t.y>ext[2])&(t.y<ext[3])]
    ax.scatter(tt.x,tt.y,marker='^',c='#1e3452',s=90 if not small else 50,edgecolor='white',zorder=5)
    lab={'WAINWRIGHT BASE':'Rogers "Wainwright Base" tower','922 3 AV':'TELUS and Rogers towers (town)'}
    for _,r in tt.iterrows():
        if r['name'] in lab: ax.annotate(lab[r['name']],(r.x,r.y),xytext=(-8,-4) if r['name']=='WAINWRIGHT BASE' else (8,-12),ha='right' if r['name']=='WAINWRIGHT BASE' else 'left',textcoords='offset points',fontsize=7.5 if not small else 6,color='white',weight='bold')
    ax.scatter(wg.x,wg.y,marker='o',facecolor='none',edgecolor='white',s=28,lw=1.1,zorder=5)
    for nm,(lo,la) in pts.items():
        x,y=tr.transform(lo,la); ax.plot(x,y,'X' if 'Pairs' in nm else 'o',color='k',mec='white',ms=11 if 'Pairs' in nm else 6,zorder=6)
        if not small: ax.annotate(nm,(x,y),xytext=(10,-22) if 'Pairs' in nm else (8,4),textcoords='offset points',fontsize=8,color='white',weight='bold')
    ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
    sb=2000; ax.plot([ext[0]+700,ext[0]+700+sb],[ext[2]+700]*2,color='white',lw=3); ax.text(ext[0]+700,ext[2]+900,'2 km',color='white',fontsize=8,weight='bold')
# main map
fig,ax=plt.subplots(figsize=(10,9),dpi=170)
ax.imshow(rgb,extent=ext); ax.imshow(cls(B['Any_network']),extent=ext,cmap=cmap,vmin=-.5,vmax=3.5,alpha=.45,interpolation='nearest')
deco(ax)
h=[mpatches.Patch(color=c,label=l,alpha=.7) for c,l in zip(cmap.colors[::-1],labels[::-1])]
h+=[plt.Line2D([],[],color='#7a1020',lw=1.5,label='Oct. 17 burn outline'),plt.Line2D([],[],marker='^',ls='',color='#1e3452',ms=9,label='Licensed cell tower'),plt.Line2D([],[],marker='o',ls='',mfc='none',mec='#333',ms=6,label='WiGLE cell sighting (2014–2020)')]
h[:0]=[]
ax.legend(handles=h,loc='upper left',fontsize=7.8,framealpha=.92,title='Chance of usable phone signal\n(any network, low-band LTE)',title_fontsize=8)
fig.text(.01,.005,'Modelled, not measured. Copernicus 30 m terrain; ISED tower sites and heights; 700 MHz; handset outdoors at 1.5 m; path-loss threshold 140 dB, 8 dB spread.\nBase image: Sentinel-2, Oct. 19, 2026 (Copernicus). Burn outline and model: Prime Rogue Inc.',fontsize=6.6,color='#333')
bx0,by0,bx1,by1=burn.total_bounds
ax.add_patch(mpatches.Rectangle((bx0-600,by0-600),bx1-bx0+1200,by1-by0+1200,fill=False,ec='white',lw=1.6,ls='--',zorder=7))
ax.text(bx0-600,by1+700,'Zoom (next map)',color='white',fontsize=8,weight='bold')
plt.savefig('fig_cell_any.png',bbox_inches='tight')
# carriers + sensitivity panels over hillshade
dem=rasterio.open('dem_cop30.tif')
fig,axs=plt.subplots(2,2,figsize=(12,11),dpi=150)
panels=[('Rogers network (threshold 140 dB)',B['Rogers'],140),('TELUS network, Bell likely shares it (140 dB)',B['TELUS'],140),('Any network, cautious (135 dB)',B['Any_network'],135),('Any network, optimistic (145 dB)',B['Any_network'],145)]
for ax,(ttl,L,M) in zip(axs.ravel(),panels):
    ax.imshow(rgb,extent=ext); ax.imshow(cls(L,M),extent=ext,cmap=cmap,vmin=-.5,vmax=3.5,alpha=.5,interpolation='nearest'); deco(ax,small=True); ax.set_title(ttl,fontsize=10,loc='left')
fig.legend(handles=h[:4],loc='lower center',ncol=4,fontsize=9,frameon=False)
plt.subplots_adjust(bottom=.05,wspace=.03,hspace=.12)
plt.savefig('fig_cell_panels.png',bbox_inches='tight'); print('ok')

# ---------- zoom on the burn ----------
bx0,by0,bx1,by1=burn.total_bounds; pad=600
zx=(bx0-pad,bx1+pad,by0-pad,by1+pad)
# add locator box to main map
fig0=plt.figure()
gray=np.dot(rgb[...,:3].astype(float),[0.299,0.587,0.114]); gray=(gray-np.percentile(gray,2))/(np.percentile(gray,98)-np.percentile(gray,2)); gray=np.clip(gray,0,1)
C=cls(B['Any_network'])
fig,axs=plt.subplots(1,2,figsize=(14,7.2),dpi=170,gridspec_kw=dict(width_ratios=[1,1]))
for ax,(ttl,L) in zip(axs,[('Any network (best case for a 911 call)',B['Any_network']),('TELUS network (Bell likely shares it)',B['TELUS'])]):
    ax.imshow(gray,extent=ext,cmap='gray',vmin=0,vmax=1.4)
    ax.imshow(cls(L),extent=ext,cmap=cmap,vmin=-.5,vmax=3.5,alpha=.78,interpolation='nearest')
    burn.boundary.plot(ax=ax,color='black',lw=2.2); burn.boundary.plot(ax=ax,color='white',lw=.9)
    for nm,(lo,la),off in [('Pairs range\n(best guess)',(-110.9005,52.7695),(12,-26)),('Range 15',(-110.89946,52.77034),(12,2)),('Range 2',(-110.90089,52.772),(-12,8))]:
        x,y=tr.transform(lo,la); ax.plot(x,y,'X' if 'Pairs' in nm else 'o',color='white' ,mec='black',mew=1.3,ms=12 if 'Pairs' in nm else 7,zorder=6)
        ax.annotate(nm,(x,y),xytext=off,textcoords='offset points',fontsize=8.5,weight='bold',color='black',ha='right' if off[0]<0 else 'left',bbox=dict(fc='white',ec='none',alpha=.85,pad=1.5),zorder=7)
    x,y=tr.transform(-110.8311,52.7060); ax.annotate('',xy=(zx[1]-300,zx[2]+300),xytext=(zx[1]-1100,zx[2]+1300),arrowprops=dict(arrowstyle='->',lw=2,color='black'))
    ax.text(zx[1]-1150,zx[2]+1350,'to Rogers "Wainwright Base"\ntower, 8.5 km SSE',fontsize=7.5,ha='right',bbox=dict(fc='white',ec='none',alpha=.85,pad=1.5))
    ax.set_xlim(zx[0],zx[1]); ax.set_ylim(zx[2],zx[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel(''); ax.set_aspect('equal')
    ax.plot([zx[0]+300,zx[0]+1300],[zx[2]+300]*2,color='black',lw=3); ax.text(zx[0]+300,zx[2]+420,'1 km',fontsize=8,weight='bold')
    ax.set_title(ttl,fontsize=11,loc='left',weight='bold')
hz=[mpatches.Patch(color=c,label=l,alpha=.85) for c,l in zip(cmap.colors[::-1],labels[::-1])]+[plt.Line2D([],[],color='black',lw=2.2,label='Oct. 17 burn outline')]
fig.legend(handles=hz,loc='lower center',bbox_to_anchor=(0.5,0.035),ncol=5,fontsize=9.5,frameon=False,title='Modelled chance of usable phone signal (700 MHz, phone outdoors)',title_fontsize=9.5)
plt.subplots_adjust(bottom=.17,wspace=.04)
fig.text(.01,.005,'Model, not measurement. Threshold 140 dB, 8 dB spread. Base: Sentinel-2, Oct. 19, 2026 (greyscale). Burn outline and model: Prime Rogue Inc.',fontsize=7,color='#333')
plt.savefig('fig_cell_zoom.png',bbox_inches='tight')
# locator box on main map
import matplotlib.image as mpimg
