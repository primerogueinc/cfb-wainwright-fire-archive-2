import csv,collections
from datetime import datetime
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.dates as md
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
plt.rcParams['font.family']='DejaVu Sans'
R=list(csv.DictReader(open('events_oct16_17.csv')))
lanes=['Public forecasts (all Alberta regions)','Warning bulletins (text, Alberta)','Alert-system messages (CAP feed, Alberta)','Lloydminster aviation forecast (TAF)','SIGMETs (Edmonton & Winnipeg FIRs)','AIRMETs (Edmonton & Winnipeg FIRs)','Pilot reports (Edmonton, Lloydminster)']
labels={'Public forecasts (all Alberta regions)':'Public forecasts\n(all Alberta regions)','Warning bulletins (text, Alberta)':'Weather warnings\n(text bulletins, Alberta)','Alert-system messages (CAP feed, Alberta)':'Alert-feed messages\n(CAP, Alberta)','Lloydminster aviation forecast (TAF)':'Aviation forecast\n(Lloydminster TAF)','SIGMETs (Edmonton & Winnipeg FIRs)':'SIGMETs: severe\nweather for aircraft','AIRMETs (Edmonton & Winnipeg FIRs)':'AIRMETs: moderate\nweather for aircraft','Pilot reports (Edmonton, Lloydminster)':'Pilot reports\n(Edmonton, Lloydminster)'}
navy,gold,red,grey='#1e3452','#d4af68','#c4552f','#a69c90'
fig,ax=plt.subplots(figsize=(13,7.2),dpi=170)
fig.patch.set_facecolor('#fbf8f2'); ax.set_facecolor('#fbf8f2')
import random; random.seed(3)
tot=collections.Counter(); hit=collections.Counter()
for r in R:
    l=r['lane']; y=len(lanes)-lanes.index(l); t=datetime.strptime(r['time_mdt'],'%Y-%m-%d %H:%M'); h=r['covers_or_near_cfb_wainwright']=='True'
    tot[l]+=1; hit[l]+=h
    j=random.uniform(-.22,.22)
    if h: ax.plot(t,y+j*.5,'o',ms=8.5,color=red,mec='white',mew=.8,zorder=4)
    else: ax.plot(t,y+j,'|',ms=13,mew=1.6,color=grey,alpha=.75,zorder=2)
for i,l in enumerate(lanes):
    y=len(lanes)-i
    note={'Lloydminster aviation forecast (TAF)':f"{tot[l]}  (Lloydminster airport,\n74 km from the fire)",'SIGMETs (Edmonton & Winnipeg FIRs)':f"{tot[l]}  ({hit[l]}: one warning,\nissued for two flight regions)",'Pilot reports (Edmonton, Lloydminster)':f"{tot[l]}  ({hit[l]} at Lloydminster)",'Public forecasts (all Alberta regions)':f"{tot[l]}  ({hit[l]} for the\nWainwright region)"}.get(l,f"{tot[l]}  ({hit[l]} covering Wainwright)")
    ax.text(1.005,y,note,transform=ax.get_yaxis_transform(),va='center',fontsize=8.5,color=navy)
ax.set_yticks(range(1,len(lanes)+1)); ax.set_yticklabels([labels[l] for l in lanes][::-1],fontsize=9)
ax.set_ylim(.4,len(lanes)+.9)
a,b=datetime(2026,10,16,0),datetime(2026,10,17,22)
ax.set_xlim(a,b)
ax.axvspan(datetime(2026,10,17,15,16),datetime(2026,10,17,15,58),color=red,alpha=.12,zorder=1)
ax.text(datetime(2026,10,17,15,37),len(lanes)+.55,'3:16 p.m. west wind hits the base (82 km/h gust)\n→ "just before 4 p.m." fire reaches soldiers',ha='center',fontsize=8,color=red,weight='bold')
ax.xaxis.set_major_locator(md.HourLocator(byhour=[0,6,12,18])); ax.xaxis.set_major_formatter(md.DateFormatter('%-I %p\n%b %-d'))
ax.tick_params(axis='x',labelsize=8.5)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.grid(axis='x',color='#e3dccf',lw=.7); ax.set_axisbelow(True)
for y in range(1,len(lanes)+1): ax.axhline(y-.5,color='#e3dccf',lw=.5)
N=len(R); H=sum(r['covers_or_near_cfb_wainwright']=='True' for r in R)
fig.text(.015,.965,f"{N} weather forecasts, warnings and alerts in 46 hours",fontsize=17,weight='bold',color=navy)
fig.text(.015,.895,f"Every Environment Canada and aviation weather message for Alberta, Oct. 16–17, 2026, by when it went out (MDT).\nRed: the {H} that covered CFB Wainwright, its forecast region or nearby Lloydminster.",fontsize=9.5,color='#3b332b')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=red,mec='white',ms=8,label='Covered CFB Wainwright, its region or Lloydminster'),Line2D([],[],marker='|',ls='',color=grey,ms=12,mew=1.6,label='Elsewhere in Alberta or the flight region')],loc='upper left',fontsize=8,frameon=False,bbox_to_anchor=(0,1.0))
fig.text(.015,.012,"None was a phone or broadcast emergency alert: every alert-feed message carried Broadcast Immediately = No and Wireless Immediate = No.\nSource: ECCC MSC Datamart (bulletins, CAP alerts), NAV CANADA/ECCC aviation bulletins. Analysis: Prime Rogue Inc., CC BY-NC 4.0.",fontsize=7.6,color='#5a5047')
try:
    im=plt.imread('/tmp/claude-0/-home-claude/11d4f8cb-573a-58d8-9723-d23e8fe40496/scratchpad/v6plan/logo.png')
    fig.figimage(im,xo=fig.bbox.xmax-im.shape[1]*0.0-0, yo=0, alpha=0) if False else None
    ab=AnnotationBbox(OffsetImage(im,zoom=60/im.shape[0]),(.965,.955),xycoords='figure fraction',frameon=False); fig.add_artist(ab)
except Exception as e: print('logo',e)
plt.subplots_adjust(left=.15,right=.86,top=.86,bottom=.13)
plt.savefig('fig_alert_flood.png',facecolor=fig.get_facecolor())
print(N,H)
