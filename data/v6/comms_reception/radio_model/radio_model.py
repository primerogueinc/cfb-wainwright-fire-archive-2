import numpy as np, pandas as pd, rasterio, json
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from pyproj import Transformer
from scipy.stats import norm
from scipy.ndimage import map_coordinates
tr=Transformer.from_crs(4326,32612,always_xy=True)
px,py=tr.transform(-110.9005,52.7695)
src=rasterio.open('dem_cop30.tif')
x0,x1,y0,y1=px-26e3,px+26e3,py-18e3,py+16e3; res=30; W=int((x1-x0)/res); H=int((y1-y0)/res); T=from_origin(x0,y1,res,res)
dem=np.zeros((H,W),'float32'); reproject(src.read(1),dem,src_transform=src.transform,src_crs=src.crs,dst_transform=T,dst_crs='EPSG:32612',resampling=Resampling.bilinear)
def z_at(x,y): return map_coordinates(dem,[(y1-np.asarray(y))/res-0.5,(np.asarray(x)-x0)/res-0.5],order=1,mode='nearest')
g=60; gx0,gx1,gy0,gy1=px-12e3,px+12e3,py-11e3,py+12e3
xs=np.arange(gx0+g/2,gx1,g); ys=np.arange(gy1-g/2,gy0,-g); GX,GY=np.meshgrid(xs,ys); rx=GX.ravel(); ry=GY.ravel()
R=6371e3*4/3
def J(v):
    out=np.zeros_like(v); m=v>-0.78; vv=v[m]; out[m]=6.9+20*np.log10(np.sqrt((vv-0.1)**2+1)+vv-0.1); return out
def loss(tx,ty,ht,hr,F,qx=None,qy=None,NS=360,clutter=6.0):
    qx=rx if qx is None else np.atleast_1d(qx); qy=ry if qy is None else np.atleast_1d(qy)
    lam=3e8/(F*1e6); zt=z_at([tx],[ty])[0]+ht; zr=z_at(qx,qy)+hr
    d=np.maximum(np.hypot(qx-tx,qy-ty),30); s=np.linspace(0,1,NS)[1:-1]
    sx=tx+np.outer(qx-tx,s); sy=ty+np.outer(qy-ty,s); zs=z_at(sx.ravel(),sy.ravel()).reshape(sx.shape)
    d1=np.outer(d,s); d2=np.outer(d,1-s); zs=zs+d1*d2/(2*R)
    h=zs-(zt+(zr-zt)[:,None]*s[None,:]); v=h*np.sqrt(2*d[:,None]/(lam*d1*d2))
    i=np.argmax(v,axis=1); vmax=v[np.arange(len(d)),i]; L=J(vmax)
    idx=np.arange(v.shape[1])[None,:]
    left=np.where(idx<i[:,None],v,-9).max(axis=1); right=np.where(idx>i[:,None],v,-9).max(axis=1)
    L=L+np.where(vmax>-0.78,0.5*(J(left)+J(right)),0)
    fspl=32.45+20*np.log10(F)+20*np.log10(d/1000); pe=40*np.log10(d)-20*np.log10(max(ht,1))-20*np.log10(hr)
    return np.maximum(fspl,pe)+L+clutter
P=lambda L,M,sig=8: norm.cdf((M-L)/sig)
# ---- radio assumptions (generic, not DND's actual equipment) ----
SENS=-113  # dBm, typical narrowband FM / tactical receiver threshold
cases={
 'hh_hh':dict(F=50,ht=1.5,hr=1.5,ptx=37,gt=-5,gr=-5,desc='Handheld to handheld, low VHF (50 MHz), 5 W'),
 'veh_hh':dict(F=50,ht=3.0,hr=1.5,ptx=47,gt=0,gr=-5,desc='Vehicle radio (50 W, 3 m mast) to handheld, low VHF'),
 'veh_veh':dict(F=50,ht=3.0,hr=3.0,ptx=47,gt=0,gr=0,desc='Vehicle to vehicle, low VHF, 50 W'),
}
res_={}
for k,c in cases.items():
    L=loss(px,py,c['ht'],c['hr'],c['F']); M=c['ptx']+c['gt']+c['gr']-SENS
    res_[k]=dict(M=M,L=L); print(k,'MAPL',M)
# assumed Range Control / garrison hub (approximate, garrison area) and highest-ground relay from terrain
gar=tr.transform(-110.885,52.825)
# highest ground within 15 km of pairs range, inside training area roughly south/west of garrison
gx,gy=np.meshgrid(np.arange(px-15e3,px+15e3,90),np.arange(py-15e3,py+10e3,90)); gz=z_at(gx.ravel(),gy.ravel())
mask=np.hypot(gx.ravel()-px,gy.ravel()-py)<=15e3; j=np.argmax(np.where(mask,gz,-1)); hill=(gx.ravel()[j],gy.ravel()[j]); hz=gz[j]
hill_ll=Transformer.from_crs(32612,4326,always_xy=True).transform(*hill)
print('highest ground',hill_ll,hz,'dist km',np.hypot(hill[0]-px,hill[1]-py)/1000)
links={}
for nm,(sx,sy),sh,sg in [('Garrison hub (20 m mast, 50 W, 3 dBi)',gar,20,3),('Hilltop relay (10 m mast, 50 W, 3 dBi)',hill,10,3)]:
    for k,c in [('handheld (5 W)',dict(ptx=37,g=-5,h=1.5)),('vehicle (50 W, 3 m)',dict(ptx=47,g=0,h=3.0))]:
        L=float(loss(sx,sy,sh,c['h'],50,qx=px,qy=py)[0])
        M=min(c['ptx']+c['g']+sg-SENS, 47+sg+c['g']-SENS)  # talk-in limited by the weaker side
        links[f'{nm} <-> {k}']=dict(loss_dB=round(L,1),MAPL=M,P=round(float(P(L,M)),2),dist_km=round(float(np.hypot(sx-px,sy-py)/1000),1))
    res_['area_'+nm]=dict(M=37-5+sg-SENS,L=loss(sx,sy,sh,1.5,50))
# civilian: Town of Wainwright fire (VHF 169 MHz) reservoir site, assume 30 m antenna, portable 5 W
fx_,fy_=tr.transform(-110.8467,52.8458)
res_['town_fire']=dict(M=37+6-2-(-116),L=loss(fx_,fy_,30,1.5,169))  # base 6 dBi, portable -2 dBi, -116 dBm
lf=float(loss(fx_,fy_,30,1.5,169,qx=px,qy=py)[0]); links['Town fire base (169 MHz, 30 m) <-> portable at pairs range']=dict(loss_dB=round(lf,1),MAPL=157,P=round(float(P(lf,157)),2),dist_km=round(float(np.hypot(fx_-px,fy_-py)/1000),1))
np.savez_compressed('radio.npz',rx=rx,ry=ry,shape=GX.shape,**{k:v['L'] for k,v in res_.items()},**{k+'_M':np.array(v['M']) for k,v in res_.items()},gar=np.array(gar),hill=np.array(hill),townfire=np.array([fx_,fy_]))
summ={}
for k,v in res_.items():
    p=P(v['L'],v['M']); d=np.hypot(rx-px,ry-py)
    summ[k]=dict(MAPL=v['M'],pct_area_ge50=round(float((p>=.5).mean()*100),1),median_reach_km_ge50=round(float(np.percentile(d[p>=.5],90)/1000),1) if (p>=.5).any() else 0)
out=dict(links=links,area=summ,hill=dict(lonlat=[round(hill_ll[0],4),round(hill_ll[1],4)],elev_m=round(float(hz),0)))
json.dump(out,open('radio_results.json','w'),indent=1); print(json.dumps(out,indent=1))
