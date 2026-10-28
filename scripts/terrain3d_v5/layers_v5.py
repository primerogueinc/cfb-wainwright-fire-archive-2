"""Layers for the Rubicon Ridge (working model v5) 3D videos.

Starts from the v4 scene layers (part 1, scripts/terrain3d/prep.py: terrain, canopy, imagery, scar, severity) and replaces
  * the places: Range 2, Range 11, Range 13 and Range 15 at their 2024 Range Standing Orders grids
    (part 1, data/reference/range_standing_orders_2024/rso_points.csv); the road head as used by v5 (Range 15's entrance,
    16 m from the Range 15 barrier). The v4 OpenStreetMap positions are kept under 'camera_ref' only, so the v4 camera
    shots (rendered around them) can be reused unchanged;
  * the fire: the weighted median arrival time of v5's 1,440 weighted scar-limited runs, per 40 m cell (cells reached in
    under half the weight are left unburned), re-run with scripts/rubicon_ridge/model_v5 and checked against the published
    ensemble_runs.csv (every run's arrivals at Range 2, Range 15, the road head and v4's point within 0.005 min).
    Before the front this is the weighted median of v4's four morning cases (start 4:00-6:30 a.m.) as carried into v5.

Run: python3 layers_v5.py V4_LAYERDIR V5_MEDIAN_NPZ OUTDIR
"""
import sys, os, json, shutil
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from pyproj import Transformer
from PIL import Image

L4, MED, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
for f in ('height.npy', 'canopy.npy', 'canopy_valid.npy', 'pre_tc.png', 'dnbr.npy', 'scar.png', 'severity.npy'):
    shutil.copy(os.path.join(L4, f), OUT)
P = json.load(open(os.path.join(L4, 'points.json')))
G = P['grid']; RES = G['res']; W, H = G['w'], G['h']
TG = from_origin(G['x0'], G['y0'] + H * RES, RES, RES)
tr = Transformer.from_crs('EPSG:4326', G['crs'], always_xy=True)
h = np.load(os.path.join(L4, 'height.npy'))


def pt(lon, lat, source):
    x, y = tr.transform(lon, lat); x -= G['x0']; y -= G['y0']
    r, c = int(H - 1 - y / RES), int(x / RES)
    return dict(lat=lat, lon=lon, x=x, y=y, z=float(h[r, c]), source=source)


RSO = 'Range Standing Orders 2024 (CFB Wainwright), grid reference converted'
P['camera_ref'] = {'Range 2': P['pov']['Range 2'], 'Roadhead': P['pov']['Roadhead'], 'Range 11': P['extra']['Range 11']}
P['pov'] = {
    'Range 13': pt(-110.9615, 52.75964, RSO + ' (WD 02598 45534)'),
    'Range 2': pt(-110.90089, 52.772, RSO + ' (WD 06686 46913)'),
    'Roadhead': pt(-110.89965, 52.77045, 'Road head, Range 15 entrance (RSO Range 15 barrier, WD 0677 4673, is 16 m away)'),
    'Garrison': P['pov']['Garrison'],
}
P['extra']['Range 11'] = dict(pt(-110.96129, 52.76796, RSO + ' (WD 02612 46459)'), note='AGL range (RSO para 51i)')
P['extra']['Range 15'] = pt(-110.89946, 52.77034, RSO + ' (WD 06783 46728)')

# ---------------------------------------------------------------- fire
z = np.load(MED)
T40 = from_origin(float(z['minx']), float(z['maxy']), float(z['res']), float(z['res']))
med = z['median'].astype(np.float32)                          # minutes after 3:00 p.m. MDT
fin = np.isfinite(med).astype(np.float32)
filled = np.where(np.isfinite(med), med, 400).astype(np.float32)


def to_grid(a, rs):
    out = np.full((H, W), np.nan, np.float32)
    reproject(a, out, src_transform=T40, src_crs=str(z['crs']), dst_transform=TG, dst_crs=G['crs'], dst_nodata=np.nan, resampling=rs)
    return out


A = to_grid(filled, Resampling.bilinear); M = to_grid(fin, Resampling.bilinear)
arr = np.where((M >= 0.5) & (A < 399), A + 15 * 60, np.nan).astype(np.float32)   # minutes after midnight MDT
scar5 = np.asarray(Image.open(os.path.join(L4, 'scar.png'))) > 127
scar = scar5.reshape(H, 2, W, 2).mean(axis=(1, 3)) >= 0.5
unk = scar & ~np.isfinite(arr)
np.save(os.path.join(OUT, 'arrival_min.npy'), arr); np.save(os.path.join(OUT, 'unknown.npy'), unk)
T_FRONT = 15 * 60 + 27.9
json.dump({'T_SHIFT': T_FRONT, 'source': 'Rubicon Ridge working model v5: weighted median of the 1,440 weighted scar-limited runs',
           'front_at_fire': '15:27.9 (moving line toward 114.5 deg, 77 km/h)'}, open(os.path.join(OUT, 'fire_meta.json'), 'w'))
json.dump(P, open(os.path.join(OUT, 'points.json'), 'w'), indent=1)
km2 = lambda m: m.sum() * RES * RES / 1e6
print(f'burned (median) {km2(np.isfinite(arr)):.2f} km2; scar {km2(scar):.2f}; scar not reached {km2(unk):.2f}')
print(f'earliest cell {np.nanmin(arr) // 60:.0f}:{np.nanmin(arr) % 60:04.1f}')
for k, v in list(P['pov'].items()) + [('Range 11', P['extra']['Range 11']), ('Range 15', P['extra']['Range 15'])]:
    r, c = int(H - 1 - v['y'] / RES), int(v['x'] / RES); a = arr[r, c]
    print(f'  {k}: ' + (f'fire arrives {int(a // 60)}:{a % 60:04.1f}' if np.isfinite(a) else 'not reached'))
