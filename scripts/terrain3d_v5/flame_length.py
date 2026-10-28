"""Flame length per 10 m cell from the v5 arrival field (Byram 1959): L = 0.0775 * I^0.46 m, I = H * w * r (kW/m).
r: local spread rate from the arrival-time gradient (r = 1 / |grad T|, smoothed over ~30 m); H = 18,000 kJ/kg and
w = 0.35 kg/m2 consumed, the values v5 uses for its GOES fire-power check (config.json frp_check).
For mid-October cured grass this gives about 0.5-2 m on slow flanks and backing fire and 5-8 m at the head of the
post-front run. Before the front the rate is capped at 2 km/h (see below). Run: python3 flame_length.py LAYERDIR  (writes flame_len.npy, ros_mpm.npy)"""
import sys, os, numpy as np
from scipy.ndimage import gaussian_filter
Ld = sys.argv[1]; RES = 10.0
a = np.load(os.path.join(Ld, 'arrival_min.npy')).astype(np.float64)
fin = np.isfinite(a)
af = np.where(fin, a, np.nan)
# fill holes for the gradient: nearest-ish by smoothing a filled copy
fill = np.where(fin, a, np.nanmax(a) + 60)
sm = gaussian_filter(fill, 1.5)
gy, gx = np.gradient(sm, RES)                   # min per m
g = np.hypot(gx, gy)
ros = np.clip(1.0 / np.maximum(g, 1e-6), 0.3, 250.0)  # m/min
ros = np.where(fin, ros, np.nan)
# before the front the median field has steps that are not spread (the 12:17/2:49 detection outline is filled at one time;
# v4's morning cases join): cap the pre-front rate at 2 km/h, about the fastest of v4's morning cases in the afternoon wind
pre_m = fin & (a < 15 * 60 + 27.9)
plateau = pre_m & (ros > 2000 / 60)            # no spread signal here (outline filled at one time): no flames drawn
ros = np.where(pre_m, np.minimum(ros, 2000 / 60), ros)
I = 18000.0 * 0.35 * ros / 60.0                 # kW/m
L = np.clip(0.0775 * I ** 0.46, 0.3, 12.0)
L[plateau] = np.nan
np.save(os.path.join(Ld, 'flame_len.npy'), L.astype(np.float32)); np.save(os.path.join(Ld, 'ros_mpm.npy'), ros.astype(np.float32))
pre = fin & (a < 15 * 60 + 27.9); post = fin & ~pre
for n, m in (('before the front', pre), ('after the front', post)):
    q = np.nanpercentile(L[m], [10, 50, 90]); r = np.nanpercentile(ros[m], [10, 50, 90]) * 0.06
    print(f'{n}: flame length p10/50/90 {q.round(1)} m; spread {r.round(2)} km/h')
