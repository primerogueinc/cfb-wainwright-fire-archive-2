import numpy as np, pandas as pd, geopandas as gpd, json
from scipy.stats import norm
from pyproj import Transformer
from shapely.geometry import Point
z=np.load('losses.npz'); shape=tuple(z['shape']); rx,ry=z['rx'],z['ry']
keys=[k for k in z.files if '|' in k]
best={c:np.min([z[k] for k in keys if k.startswith(c+'|')],axis=0) for c in ('Rogers','TELUS')}
best['Any network']=np.minimum(best['Rogers'],best['TELUS'])
SIG=8.0
def P(L,M=140): return norm.cdf((M-L)/SIG)
tr=Transformer.from_crs(4326,32612,always_xy=True)
pts={'Pairs range (best guess)':(-110.9005,52.7695),'Range 15':(-110.89946,52.77034),'Range 2':(-110.90089,52.772),'Road head':(-110.89965,52.77045)}
def at(arr,lon,lat):
    x,y=tr.transform(lon,lat); i=np.argmin((rx-x)**2+(ry-y)**2); return arr[i]
out={}
for nm,(lo,la) in pts.items():
    out[nm]={c:{'loss_dB':round(float(at(best[c],lo,la)),1),'P140':round(float(P(at(best[c],lo,la))),2),'P135':round(float(P(at(best[c],lo,la),135)),2),'P145':round(float(P(at(best[c],lo,la),145)),2)} for c in best}
# which tower serves pairs range
x,y=tr.transform(-110.9005,52.7695); i=np.argmin((rx-x)**2+(ry-y)**2)
out['per_tower_at_pairs_range']={k:round(float(z[k][i]),1) for k in keys}
# burn area classes
burn=gpd.read_file('/home/claude/cfb-wainwright-fire-archive-2/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson').to_crs(32612).geometry.union_all()
from shapely import contains_xy
inb=contains_xy(burn,rx,ry)
def classes(p): return {'likely(>=90%)':round(float((p>=.9).mean()*100),1),'probable(50-90%)':round(float(((p>=.5)&(p<.9)).mean()*100),1),'marginal(10-50%)':round(float(((p>=.1)&(p<.5)).mean()*100),1),'unlikely(<10%)':round(float((p<.1).mean()*100),1)}
out['burn_area_classes_%']={c:classes(P(best[c][inb])) for c in best}
out['grid_classes_%']={c:classes(P(best[c])) for c in best}
# WiGLE validation (TELUS sightings)
w=pd.read_csv('../radio/wigle_points.csv'); w=w[w.kind=='cell']
val=[]
for _,r in w.iterrows():
    car='TELUS' if 'TELUS' in r.operator.upper() else 'Any network'
    x,y=tr.transform(r.lon,r.lat)
    if not(rx.min()<=x<=rx.max() and ry.min()<=y<=ry.max()): continue
    val.append((r.operator,r.tech,round(float(P(at(best['TELUS'],r.lon,r.lat))),2)))
out['wigle_validation_TELUS_P140']=val
np.savez_compressed('best.npz',rx=rx,ry=ry,shape=shape,**{k.replace(' ','_'):v for k,v in best.items()})
json.dump(out,open('results.json','w'),indent=1); print(json.dumps(out,indent=1))
