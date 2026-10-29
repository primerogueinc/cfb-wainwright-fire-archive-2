import glob,re,csv,sys
import xml.etree.ElementTree as ET
ns={'c':'urn:oasis:names:tc:emergency:cap:1.2'}
rows=[]
for f in sorted(glob.glob('cap/*.cap')):
    t=ET.parse(f).getroot()
    g=lambda e,k:(e.findtext('c:'+k,namespaces=ns) or '').strip()
    for info in t.findall('c:info',ns):
        if g(info,'language')!='en-CA': continue
        P={p.findtext('c:valueName',namespaces=ns):p.findtext('c:value',namespaces=ns) for p in info.findall('c:parameter',ns)}
        areas=[a.findtext('c:areaDesc',namespaces=ns) for a in info.findall('c:area',ns)]
        wain=[a for a in areas if re.search('Wainwright|Lloydminster|Vermilion|Provost|Edgerton|Irma|Chauvin',a or '')]
        rows.append(dict(file=f.split('/')[-1],sent=g(t,'sent'),status=g(t,'status'),msgType=g(t,'msgType'),
          event=g(info,'event'),urgency=g(info,'urgency'),severity=g(info,'severity'),certainty=g(info,'certainty'),
          effective=g(info,'effective'),expires=g(info,'expires'),headline=g(info,'headline'),
          broadcast_immediately=P.get('layer:SOREM:1.0:Broadcast_Immediately',''),
          wireless_immediate=P.get('layer:SOREM:2.0:WirelessImmediate',''),
          alert_type=P.get('layer:EC-MSC-SMC:1.0:Alert_Type',''),alert_name=P.get('layer:EC-MSC-SMC:1.1:Alert_Name','') or P.get('layer:EC-MSC-SMC:1.0:Alert_Name',''),
          n_areas=len(areas),wainwright_areas='; '.join(wain),params='|'.join(f'{k}={v}' for k,v in P.items())))
with open('cap_parsed.csv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r['sent'],r['msgType'],r['event'],r['severity'],r['urgency'],'BI=',r['broadcast_immediately'],'WI=',r['wireless_immediate'],'| W:',(r['wainwright_areas'][:90] or '-'),'|',r['file'][:40])
