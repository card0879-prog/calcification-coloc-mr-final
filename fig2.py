import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
from matplotlib.lines import Line2D
from figstyle import V

# verdict -> (marker, facecolor, edgecolor, size)
ST={
 'PASS'       :('o',V['PASS_strong'],'none',190),
 'moderate'   :('o',V['PASS_moderate'],'#2E8B70',150),
 'failcoloc'  :('D',V['fail_coloc'],'none',150),
 'noncoloc'   :('^',V['noncoloc'],'none',200),
 'reversed'   :('v',V['noncoloc'],'none',200),
 'uninstr'    :('X',V['uninstr'],'none',170),
 'underpowered':('s','white','#9A9A9A',150),
 'ref'        :('o','white','#9A9A9A',150),
}
cols=['circulating\nactivity / proxy','whole-blood\nexpression','disease-tissue\nexpression\n(aorta/coronary)','source-tissue\nexpression (liver)']
rows=['LPA  (positive control)','ALPL','ENPP1','ANKH','PHOSPHO1','GGT1','UGT1A1','GOT1']
cells={
 (0,0):'PASS',
 (1,0):'failcoloc',(1,1):'noncoloc',(1,2):'uninstr',(1,3):'underpowered',
 (2,1):'noncoloc',(2,2):'failcoloc',
 (3,1):'noncoloc',(3,2):'uninstr',
 (4,2):'uninstr',
 (5,0):'ref',(5,1):'reversed',(5,3):'underpowered',
 (6,0):'ref',(6,1):'uninstr',(6,3):'underpowered',
 (7,0):'ref',(7,1):'underpowered',(7,3):'underpowered',
}
nR,nC=len(rows),len(cols)
fig,ax=plt.subplots(figsize=(10.6,7.6))
fig.suptitle('Disease-relevant operationalizations were least instrumentable in this audit',
             fontsize=13.5,fontweight='bold',y=0.975)

def yr(r): return nR-1-r  # row 0 at top
# shade the disease-tissue column (most relevant)
ax.axvspan(2-0.5,2+0.5,color='#F3F6F9',zorder=0)
# grid
for r in range(nR):
    for c in range(nC):
        ax.add_patch(plt.Rectangle((c-0.5,yr(r)-0.5),1,1,fill=False,ec='#E2E2E2',lw=0.8,zorder=1))
# markers
for (r,c),v in cells.items():
    m,fc,ec,sz=ST[v]
    ax.scatter(c,yr(r),marker=m,s=sz,facecolor=fc,edgecolor=ec,linewidths=1.4,zorder=3)
# blank cells -> faint dash
for r in range(nR):
    for c in range(nC):
        if (r,c) not in cells: ax.text(c,yr(r),'\u2013',ha='center',va='center',color='#CFCFCF',fontsize=11,zorder=2)

ax.set_xticks(range(nC)); ax.set_xticklabels(cols,fontsize=9.3)
ax.set_yticks([yr(r) for r in range(nR)]); ax.set_yticklabels(rows,fontsize=9.6)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlim(-0.6,nC-0.4); ax.set_ylim(-0.9,nR-0.4)

# group brackets on the left
def grp(r0,r1,label,color):
    x=-1.62
    ax.plot([x,x],[yr(r1)-0.42,yr(r0)+0.42],color=color,lw=2.4,clip_on=False)
    ax.text(x-0.06,(yr(r0)+yr(r1))/2,label,rotation=90,va='center',ha='right',fontsize=8.8,color=color,fontweight='bold',clip_on=False)
grp(0,0,'control','#2E5496')
grp(1,4,'PPi\u2013TNAP axis','#2E5496')
grp(5,7,'activity vs expression','#2E5496')

# convenient -> relevant arrow under header
ax.annotate('',xy=(3.35,nR-0.25),xytext=(-0.35,nR-0.25),
            arrowprops=dict(arrowstyle='->',color='#888',lw=1.3),clip_on=False)
ax.text(0.0,nR-0.05,'more convenient / powered',fontsize=8.2,color='#888',ha='left',style='italic')
ax.text(3.0,nR-0.05,'more disease-relevant',fontsize=8.2,color='#888',ha='right',style='italic')

# legend
leg=[('PASS / compatible','PASS'),('moderate / partial','moderate'),('fail-coloc','failcoloc'),
     ('non-colocalizing nominal','noncoloc'),('directionally reversed','reversed'),
     ('uninstrumentable','uninstr'),('underpowered','underpowered'),('reference exposure','ref')]
handles=[Line2D([0],[0],marker=ST[k][0],color='none',markerfacecolor=ST[k][1],
         markeredgecolor=ST[k][2] if ST[k][2]!='none' else ST[k][1],markersize=11,markeredgewidth=1.3,label=lab)
         for lab,k in leg]
ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(0.5,-0.085),ncol=4,
          frameon=False,fontsize=8.6,handletextpad=0.4,columnspacing=1.4)
ax.text(1.5,-0.78,'Each cell summarizes the verdict of using that operationalization for that target '
        '(integrating instrumentability, MR association, colocalization and direction).',
        ha='center',va='center',fontsize=7.8,color='#777',style='italic')

fig.savefig('Figure2.pdf',bbox_inches='tight')
fig.savefig('Figure2.png',bbox_inches='tight',dpi=200)
print("Figure 2 done")
