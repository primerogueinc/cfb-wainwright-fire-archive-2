import glob,re,csv,os
from datetime import datetime,timedelta,timezone
from shapely.geometry import LineString,Point
from pyproj import Transformer
MDT=timezone(timedelta(hours=-6))
def bt(fn):  # DDHHMM in filename, Oct 2026 UTC
    m=re.search(r'_(\d{2})(\d{2})(\d{2})_',fn); d,h,mi=map(int,m.groups())
    return datetime(2026,10,d,h,mi,tzinfo=timezone.utc).astimezone(MDT)
E=[]
txt=lambda f: open(f,errors='ignore').read()
for f in sorted(glob.glob('bul/FPCN1[1-6]_CWWG_*')):
    E.append(('Public forecasts\n(all Alberta regions)',bt(f),'FPCN16' in f,os.path.basename(f)))
for f in sorted(glob.glob('bul/WWCN1[0-9]_CWWG_*')):
    E.append(('Warning bulletins\n(text, Alberta)',bt(f),'wainwright' in txt(f).lower(),os.path.basename(f)))
seen=set()
for r in csv.DictReader(open('cap_parsed.csv')):
    if r['file'] in seen: continue
    seen.add(r['file'])
    t=datetime.fromisoformat(r['sent'].replace('-00:00','+00:00')).astimezone(MDT)
    E.append(('Alert-system messages\n(CAP feed, Alberta)',t,bool(r['wainwright_areas']) and 'Wainwright' in r['wainwright_areas'],r['file']))
for f in sorted(glob.glob('bul/FTCN32_CWAO_*')):
    for m in re.finditer(r'TAF (AMD )?CYLL (\d{2})(\d{2})(\d{2})Z',txt(f)):
        d,h,mi=map(int,m.groups()[1:]); E.append(('Lloydminster aviation\nforecast (TAF)',datetime(2026,10,d,h,mi,tzinfo=timezone.utc).astimezone(MDT),True,os.path.basename(f)+(' AMD' if m.group(1) else '')))
tr=Transformer.from_crs(4326,3400,always_xy=True); site=Point(tr.transform(-110.9005,52.7695))
sig=set()
for f in sorted(glob.glob('bul/WSCN2[0-9]_CWAO_*')):
    s=txt(f); m=re.search(r'(CZEG|CZWG) SIGMET (\w+) VALID (\d{6})',s)
    if not m: continue
    key=(m.group(1),m.group(2),m.group(3))
    if key in sig: continue
    sig.add(key)
    co=re.findall(r'N(\d{2})(\d{2}) W(\d{3})(\d{2})',s); near=False
    if len(co)>=2:
        L=LineString([tr.transform(-(int(a)+int(b)/60),int(c)+int(d)/60) for c,d,a,b in co])
        near=L.distance(site)/1852<=30
    E.append(('SIGMETs\n(Edmonton & Winnipeg FIRs)',bt(f),near,os.path.basename(f)))
for f in sorted(glob.glob('bul/WACN2[0-9]_CWAO_*')):
    if re.search('CZEG|CZWG',txt(f)): E.append(('AIRMETs\n(Edmonton & Winnipeg FIRs)',bt(f),False,os.path.basename(f)))
for f in sorted(glob.glob('bul/UACN*')):
    E.append(('Pilot reports\n(Edmonton, Lloydminster)',bt(f),'CYLL' in f,os.path.basename(f)))
a,b=datetime(2026,10,16,0,0,tzinfo=MDT),datetime(2026,10,17,22,0,tzinfo=MDT)
E=[e for e in E if a<=e[1]<=b]
with open('events_oct16_17.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['lane','time_mdt','covers_or_near_cfb_wainwright','source_file']); [w.writerow([l.replace('\n',' '),t.strftime('%Y-%m-%d %H:%M'),h,f]) for l,t,h,f in E]
import collections
c=collections.Counter(e[0].replace('\n',' ') for e in E); print(len(E)); print(c); print('hit',sum(e[2] for e in E))
for e in E:
    if e[2]: print(e[0].split('\n')[0],e[1].strftime('%d %H:%M'),e[3])
