import numpy as np, json
exec(open('radio_model.py').read().split('# ---- radio assumptions')[0])
import geopandas as gpd
from shapely import contains_xy
burn=gpd.read_file('/home/claude/cfb-wainwright-fire-archive-2/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson').to_crs(32612).geometry.union_all()
inb=contains_xy(burn,rx,ry); d=np.hypot(rx-px,ry-py)
SENS=-113
pts={'Pairs range (best guess)':(-110.9005,52.7695),'Range 2':(-110.90089,52.772),'Range 15':(-110.89946,52.77034),'West basin (inside burn)':(-110.935,52.765),'South lakes edge (inside burn)':(-110.905,52.748)}
ptxy={k:tr.transform(*v) for k,v in pts.items()}
def stats(L,M):
    p=P(L,M); ok=p>=.5
    cls=lambda q:{'likely':round(float((q>=.9).mean()*100),1),'probable':round(float(((q>=.5)&(q<.9)).mean()*100),1),'marginal':round(float(((q>=.1)&(q<.5)).mean()*100),1),'unlikely':round(float((q<.1).mean()*100),1)}
    reach=[round(float(np.percentile(d[ok],q)/1000),1) for q in (50,90)] if ok.any() else [0,0]
    pp={}
    for k,(x,y) in ptxy.items():
        i=np.argmin((rx-x)**2+(ry-y)**2); m=np.hypot(rx-x,ry-y)<=300
        pp[k]=dict(point=round(float(p[i]),2),median300=round(float(np.median(p[m])),2))
    return dict(MAPL=M,burn=cls(p[inb]),reach_km_p50_p90_of_ok=reach,points=pp)
out={}
base=dict(hh_hh=(50,1.5,1.5,37-5-5-SENS),veh_hh=(50,3.0,1.5,47+0-5-SENS),veh_veh=(50,3.0,3.0,47-SENS))
for k,(F,ht,hr,M) in base.items():
    L=loss(px,py,ht,hr,F)
    out[k]={'base':stats(L,M),'M-6dB':stats(L,M-6),'M+6dB':stats(L,M+6)}
    if k=='hh_hh':
        L2=loss(px,py,ht,hr,150); out[k]['VHF-high 150 MHz']=stats(L2,M)
        L3=loss(px,py,ht,hr,400); out[k]['UHF 400 MHz']=stats(L3,M)
# garrison hub variants (handheld talk-in, MAPL 148)
gar=tr.transform(-110.885,52.825)
for nm,(gx_,gy_,mast) in {'hub_20m':(*gar,20),'hub_10m':(*gar,10),'hub_40m':(*gar,40),'hub_2km_closer':(gar[0],gar[1]-2000,20),'hub_2km_farther':(gar[0],gar[1]+2000,20)}.items():
    L=loss(gx_,gy_,mast,1.5,50); out[nm]=stats(L,148)
fx_,fy_=tr.transform(-110.8467,52.8458)
out['town_fire_169']=stats(loss(fx_,fy_,30,1.5,169),157)
out['town_fire_169_15m']=stats(loss(fx_,fy_,15,1.5,169),157)
json.dump(out,open('radio_sensitivity.json','w'),indent=1)
for k,v in out.items():
    if 'base' in v: print(k,{kk:(vv['burn'],vv['reach_km_p50_p90_of_ok'],vv['points']['Pairs range (best guess)']) for kk,vv in v.items()})
    else: print(k,v['burn'],v['reach_km_p50_p90_of_ok'],{p:v['points'][p]['median300'] for p in v['points']})
