import json,os,numpy as np,geopandas as gpd,rasterio
from rasterio.windows import from_bounds
from scipy.stats import norm
from scipy.ndimage import median_filter
from pyproj import Transformer
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap; import matplotlib.patches as mpatches
plt.rcParams['font.family']='DejaVu Sans'; os.environ['CURL_CA_BUNDLE']='/root/.ccr/ca-bundle.crt'
tr=Transformer.from_crs(4326,32612,always_xy=True); px,py=tr.transform(-110.9005,52.7695)
Z=np.load('radio.npz'); sh=tuple(Z['shape']); rx,ry=Z['rx'],Z['ry']
CROP=int(4260/60); ext=(rx.min()-30+4260,rx.max()+30,ry.min()-30,ry.max()+30)
r=json.load(open('../v6/s2_lens/search.json')); f=[f for f in r['features'] if f['id']=='S2B_12UWD_20261019_0_L2A'][0]
with rasterio.open(f['assets']['visual']['href']) as s:
    w=from_bounds(ext[0],ext[2],ext[1],ext[3],s.transform); rgb=s.read(window=w,boundless=True).transpose(1,2,0).astype(float)
gray=rgb[...,:3]@[0.299,0.587,0.114]; gray=np.clip((gray-np.percentile(gray[gray>0],2))/(np.percentile(gray,98)-np.percentile(gray[gray>0],2)),0,1)
burn=gpd.read_file('/home/claude/cfb-wainwright-fire-archive-2/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson').to_crs(32612)
cmap=ListedColormap(['#b2182b','#f4a582','#92c5de','#2166ac']); labels=['Unlikely (<10%)','Marginal (10–50%)','Probable (50–90%)','Likely (≥90%)']
def cls(L,M):
    L=median_filter(L.reshape(sh),3)[:,CROP:]; p=norm.cdf((M-L)/8); c=np.zeros_like(p); c[p>=.1]=1; c[p>=.5]=2; c[p>=.9]=3; return c
gar=Z['gar']; tf=Z['townfire']
panels=[('hh_hh','A. Handheld to handheld from the pairs range\n(low VHF, 5 W)','pr'),('veh_hh','B. Vehicle radio at the pairs range to handhelds\n(low VHF, 50 W, 3 m antenna)','pr'),
        ('area_Garrison hub (20 m mast, 50 W, 3 dBi)','C. Handhelds talking to a base station\nassumed in the garrison area','gar'),('town_fire','D. Town of Wainwright fire base station to\nfirefighters\' portables (169 MHz)','tf')]
fig,axs=plt.subplots(2,2,figsize=(13,13),dpi=150)
for ax,(k,ttl,src) in zip(axs.ravel(),panels):
    ax.imshow(gray,extent=ext,cmap='gray',vmin=0,vmax=1.35)
    ax.imshow(cls(Z[k],float(Z[k+'_M'])),extent=ext,cmap=cmap,vmin=-.5,vmax=3.5,alpha=.72,interpolation='nearest')
    burn.boundary.plot(ax=ax,color='black',lw=1.6); burn.boundary.plot(ax=ax,color='white',lw=.6)
    ax.plot(px,py,'X',color='white',mec='black',mew=1.3,ms=12,zorder=6)
    if src=='gar': ax.plot(*gar,'s',color='white',mec='black',ms=9,zorder=6); ax.annotate('assumed hub',gar,xytext=(8,4),textcoords='offset points',fontsize=8,weight='bold',bbox=dict(fc='white',ec='none',alpha=.8,pad=1))
    if src=='tf': ax.plot(*tf,'^',color='white',mec='black',ms=10,zorder=6); ax.annotate('Wainwright fire\nradio site',tf,xytext=(8,4),textcoords='offset points',fontsize=8,weight='bold',bbox=dict(fc='white',ec='none',alpha=.8,pad=1))
    ax.annotate('pairs range',(px,py),xytext=(10,-14),textcoords='offset points',fontsize=8,weight='bold',bbox=dict(fc='white',ec='none',alpha=.8,pad=1))
    ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3]); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
    ax.plot([ext[0]+800,ext[0]+2800],[ext[2]+800]*2,color='black',lw=3); ax.text(ext[0]+800,ext[2]+1100,'2 km',fontsize=8,weight='bold')
    ax.set_title(ttl,fontsize=10.5,loc='left',weight='bold')
h=[mpatches.Patch(color=c,label=l,alpha=.85) for c,l in zip(cmap.colors[::-1],labels[::-1])]+[plt.Line2D([],[],color='black',lw=1.8,label='Oct. 17 burn outline')]
fig.legend(handles=h,loc='lower center',bbox_to_anchor=(.5,.025),ncol=5,fontsize=10,frameon=False,title='Modelled chance of a usable radio link',title_fontsize=10)
plt.subplots_adjust(bottom=.09,wspace=.04,hspace=.1)
fig.text(.01,.005,'Model, not measurement. Generic radios, not DND\'s actual equipment or frequencies. Terrain: Copernicus 30 m; base: Sentinel-2, Oct. 19, 2026 (greyscale). Prime Rogue Inc.',fontsize=7.5,color='#333')
plt.savefig('fig_radio_panels.png',bbox_inches='tight'); print('ok')
