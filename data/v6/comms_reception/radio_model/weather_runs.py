import numpy as np, json, pandas as pd
src_code=open('radio_model.py').read().split('# ---- radio assumptions')[0]
src_code=src_code.replace("R=6371e3*4/3","K=4/3\nR=6371e3*K").replace("zs=zs+d1*d2/(2*R)","zs=zs+d1*d2/(2*6371e3*K)")
exec(src_code)
import geopandas as gpd
from shapely import contains_xy
burn=gpd.read_file('/home/claude/cfb-wainwright-fire-archive-2/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson').to_crs(32612).geometry.union_all()
inb=contains_xy(burn,rx,ry); d=np.hypot(rx-px,ry-py)
i_pr=np.argmin((rx-px)**2+(ry-py)**2); m300=np.hypot(rx-px,ry-py)<=300
SENS=-113
gar=tr.transform(-110.885,52.825); tf=tr.transform(-110.8467,52.8458)
t=pd.read_csv('towers.csv'); t['x'],t['y']=tr.transform(t.lon.values,t.lat.values)
scen={'Baseline (standard air, before the fire)':dict(K=4/3,fire=0),
      'Fire-time air (k = 1.2, Oct. 17 6 p.m. sounding)':dict(K=1.2,fire=0),
      'Fire-time air + active fire (+6 dB inside burn)':dict(K=1.2,fire=6),
      'Morning air (inversion, k = 4, for comparison)':dict(K=4.0,fire=0)}
RESULTS={}
for nm,c in scen.items():
    K=c['K']; extra=np.where(inb,c['fire'],0.0)
    r={}
    def summ(L,M,label):
        L=L+extra; p=P(L,M); ok=p>=.5
        r[label]=dict(pairs_median300=round(float(np.median(p[m300])),2),burn_ge50_pct=round(float(ok[inb].mean()*100),1),reach_p90_km=round(float(np.percentile(d[ok],90)/1000),1) if ok.any() else 0)
    summ(loss(px,py,1.5,1.5,50),140,'Handheld to handheld (low VHF)')
    summ(loss(px,py,3.0,1.5,50),155,'Vehicle to handheld (low VHF)')
    summ(loss(gar[0],gar[1],20,1.5,50),148,'Handheld to garrison base station (assumed)')
    summ(loss(tf[0],tf[1],30,1.5,169),157,'Wainwright fire base to portables (169 MHz)')
    for car in ('Rogers','TELUS'):
        Ls=[loss(rw.x,rw.y,0.9*rw.struct_h,1.5,700,clutter=8.0) for _,rw in t[t.carrier==car].iterrows()]
        summ(np.min(Ls,axis=0),140,f'Phone, {car} (700 MHz)')
    RESULTS[nm]=r; print(nm); [print('   ',k,v) for k,v in r.items()]
json.dump(RESULTS,open('weather_runs.json','w'),indent=1)
