import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from pyproj import Transformer
from scipy.stats import norm
from scipy.ndimage import map_coordinates
tr=Transformer.from_crs(4326,32612,always_xy=True)
px,py=tr.transform(-110.9005,52.7695)
# reproject DEM to UTM 30 m covering towers
src=rasterio.open('dem_cop30.tif')
x0,x1,y0,y1=px-26e3,px+26e3,py-18e3,py+16e3
res=30; W=int((x1-x0)/res); H=int((y1-y0)/res); T=from_origin(x0,y1,res,res)
dem=np.zeros((H,W),'float32')
reproject(src.read(1),dem,src_transform=src.transform,src_crs=src.crs,dst_transform=T,dst_crs='EPSG:32612',resampling=Resampling.bilinear)
def z_at(x,y):
    c=(x-x0)/res-0.5; r=(y1-y)/res-0.5
    return map_coordinates(dem,[r,c],order=1,mode='nearest')
# receiver grid 60 m around site
gx0,gx1,gy0,gy1=px-9e3,px+9e3,py-8e3,py+8e3; g=60
xs=np.arange(gx0+g/2,gx1,g); ys=np.arange(gy1-g/2,gy0,-g); GX,GY=np.meshgrid(xs,ys)
rx=GX.ravel(); ry=GY.ravel(); zr=z_at(rx,ry)+1.5
F=700.0; lam=3e8/(F*1e6); K=4/3; R=6371e3*K
def J(v):
    out=np.zeros_like(v); m=v>-0.78; vv=v[m]; out[m]=6.9+20*np.log10(np.sqrt((vv-0.1)**2+1)+vv-0.1); return out
def loss_from(tx,ty,ht_agl,NS=320):
    zt=z_at(np.array([tx]),np.array([ty]))[0]+ht_agl
    d=np.hypot(rx-tx,ry-ty); d=np.maximum(d,30)
    s=np.linspace(0,1,NS)[1:-1]          # interior samples
    sx=tx+np.outer(rx-tx,s); sy=ty+np.outer(ry-ty,s)
    zs=z_at(sx.ravel(),sy.ravel()).reshape(sx.shape)
    d1=np.outer(d,s); d2=np.outer(d,1-s)
    zs=zs+d1*d2/(2*R)                    # earth bulge
    los=zt+(zr-zt)[:,None]*s[None,:]
    h=zs-los
    v=h*np.sqrt(2*d[:,None]/(lam*d1*d2))
    i=np.argmax(v,axis=1); vmax=v[np.arange(len(d)),i]
    L=J(vmax)
    # Deygout secondary edges (approx: max v on each side relative to main edge, half weight)
    idx=np.arange(v.shape[1])[None,:]
    left=np.where(idx<i[:,None],v,-9).max(axis=1); right=np.where(idx>i[:,None],v,-9).max(axis=1)
    L=L+np.where(vmax>-0.78,0.5*(J(np.maximum(left,-9))+J(np.maximum(right,-9))),0)
    fspl=32.45+20*np.log10(F)+20*np.log10(d/1000)
    pe=40*np.log10(d)-20*np.log10(max(ht_agl,1))-20*np.log10(1.5)
    return np.maximum(fspl,pe)+L+8.0   # +8 dB clutter/body allowance
t=pd.read_csv('towers.csv'); t['x'],t['y']=tr.transform(t.lon.values,t.lat.values)
Ls={}
for _,r in t.iterrows():
    Ls[(r.carrier,r['name'])]=loss_from(r.x,r.y,0.9*r.struct_h); print(r['name'],r.carrier,'min loss',Ls[(r.carrier,r['name'])].min().round(1))
np.savez_compressed('losses.npz',rx=rx,ry=ry,shape=GX.shape,**{f'{c}|{n}':v for (c,n),v in Ls.items()})
