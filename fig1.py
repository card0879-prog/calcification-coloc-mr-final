import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from figstyle import V, EXPO

fig,ax=plt.subplots(figsize=(11.6,7.4)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
fig.suptitle('The pre-specified validity gate separates instrument availability from target validity',
             fontsize=14,fontweight='bold',y=0.975)

def box(x,y,w,h,text,fc,ec,tc='#222222',fs=9.5,bold=False,style='round,pad=0.012'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=style,fc=fc,ec=ec,lw=1.3,mutation_aspect=1))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color=tc,
            fontweight='bold' if bold else 'normal',zorder=5)
def arrow(x1,y1,x2,y2,c='#444444',lw=1.6,ls='-'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,lw=lw,color=c,ls=ls,zorder=1))

# ---- left: operationalizations ----
ax.text(0.135,0.90,'Operationalizations\nof a target',ha='center',va='center',fontsize=10.5,fontweight='bold',color=EXPO)
ops=['circulating enzyme activity\nor biomarker proxy','whole-blood expression',
     'disease-relevant tissue\nexpression','source-tissue expression']
oy=[0.71,0.575,0.44,0.305]
for t,y in zip(ops,oy): box(0.02,y,0.23,0.105,t,'#DCEAF6',EXPO,fs=9)
# funnel arrows into gate
for y in oy: arrow(0.25,y+0.052,0.305,0.555,c='#9FB6C9',lw=1.2)

# ---- middle: gate decision cascade ----
ax.text(0.46,0.93,'Pre-specified validity gate',ha='center',va='center',fontsize=11,fontweight='bold')
gx,gw,gh=0.315,0.29,0.082
gy=[0.78,0.625,0.47,0.315]
gtxt=['Instrumentable in this operationalization?','cis-MR association?',
      'Colocalizes with outcome?  (PP.H4)','Direction biologically coherent\n& consistent across operationalizations?']
for t,y in zip(gtxt,gy): box(gx,y,gw,gh,t,'#FFFFFF','#555555',fs=9.3,bold=True)
# yes arrows down
for i in range(3):
    arrow(gx+gw/2,gy[i],gx+gw/2,gy[i+1]+gh)
    ax.text(gx+gw/2+0.012,gy[i]-0.018,'yes',fontsize=8,color='#2E7D32',ha='left',va='center',style='italic')
ax.text(0.46,0.885,'(insufficient eQTL strength or expression \u2192 underpowered, not null)',fontsize=8,color='#888888',ha='center',va='center',style='italic')

# ---- right: verdict chips reached by 'no' ----
vx,vw,vh=0.70,0.275,0.072
verds=[('uninstrumentable',V['uninstr']),
       ('instrumented, no MR association',V['uninstr']),
       ('fail-coloc / non-colocalizing nominal',V['fail_coloc']),
       ('directionally reversed',V['noncoloc'])]
for (t,c),y in zip(verds,gy):
    cy=y+gh/2
    box(vx,cy-vh/2,vw,vh,t,c,'none',tc='white',fs=8.8,bold=True)
    arrow(gx+gw,cy,vx,cy,c='#999999',lw=1.4)
    ax.text((gx+gw+vx)/2,cy+0.013,'no',fontsize=8,color='#B0413E',ha='center',va='bottom',style='italic')

# ---- PASS (all yes) ----
arrow(gx+gw/2,gy[3],gx+gw/2,0.205+0.072)
ax.text(gx+gw/2+0.012,gy[3]-0.02,'yes',fontsize=8,color='#2E7D32',ha='left',va='center',style='italic')
box(gx,0.205,0.40,0.072,'PASS \u2014 strong (PP.H4 \u2265 0.8)   /   moderate (0.5\u20130.8)',V['PASS_strong'],'none',tc='white',fs=9.0,bold=True)

# ---- bottom calibration strip ----
ax.add_patch(FancyBboxPatch((0.02,0.035),0.95,0.085,boxstyle='round,pad=0.008',fc='#F4F7FA',ec='#CCD6DF',lw=1))
ax.text(0.04,0.078,'Calibration:',fontsize=9.5,fontweight='bold',color='#222',ha='left',va='center')
ax.text(0.04,0.055,'Lp(a) admitted;   PPi\u2013TNAP not supported under available operationalizations.',
        fontsize=9,color='#333',ha='left',va='center')

fig.savefig('Figure1.pdf',bbox_inches='tight')
fig.savefig('Figure1.png',bbox_inches='tight',dpi=200)
print("Figure 1 done")
