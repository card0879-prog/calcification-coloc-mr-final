import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
from matplotlib.gridspec import GridSpec
from figstyle import V

def load(f): d=np.load(f); return d['za'],d['ze'],d['r2']
panels=[
 ('data/f4_eqtlgen.npz','A   GGT: activity vs blood expression (eQTLGen)',
   'noncoloc','same lead-LD, opposite direction',
   'lead-pair $r^2$ = 0.99\ncross-variant $r$ = \u22120.89\nPP.H4 = 0.17 (broad signal)'),
 ('data/f4_gtexblood.npz','B   GGT: activity vs blood expression (GTEx)',
   'fail_coloc','strong expression signal,\nno shared activity architecture',
   'lead-pair $r^2$ = 0.01\ncross-variant $r$ = \u22120.31\nPP.H3 = 0.35, PP.H4 = 0.13'),
 ('data/f4_liver.npz','C   GGT: activity vs liver expression (GTEx)',
   'underpowered','liver eQTL under-powered\n(GTEx liver n = 208)',
   'liver lead |$Z$| \u2248 4\nconcordance not established'),
 ('data/f4_ugt.npz','D   Bilirubin vs UGT1A1 liver expression (GTEx)',
   'underpowered','illustrative only\n(structurally complex UGT1A locus)',
   'cross-variant $r$ = \u22120.86\nPP.H4 = 0.94 (complex locus;\nleads $r^2$ \u2248 0.05, illustrative)'),
]

fig=plt.figure(figsize=(11.2,8.6))
gs=GridSpec(2,2,figure=fig,hspace=0.40,wspace=0.30,left=0.075,right=0.90,top=0.88,bottom=0.085)
fig.suptitle('Activity\u2013expression concordance depends on tissue and eQTL resource',
             fontsize=13.2,fontweight='bold',y=0.965)

sc=None
for k,(f,title,vd,tag,stat) in enumerate(panels):
    ax=fig.add_subplot(gs[k//2,k%2])
    za,ze,r2=load(f)
    m=np.isfinite(za)&np.isfinite(ze)
    za,ze,r2=za[m],ze[m],np.clip(r2[m],0,1)
    order=np.argsort(r2)  # high-LD on top
    sc=ax.scatter(za[order],ze[order],c=r2[order],cmap='viridis',s=16,alpha=0.85,
                  edgecolors='none',vmin=0,vmax=1,zorder=3)
    ax.axhline(0,color='#BBBBBB',lw=0.7,zorder=1); ax.axvline(0,color='#BBBBBB',lw=0.7,zorder=1)
    # fit on strong activity variants
    s=np.abs(za)>5
    if s.sum()>5:
        b,a=np.polyfit(za[s],ze[s],1); xx=np.array([za.min(),za.max()])
        ax.plot(xx,a+b*xx,color=V[vd],lw=2,zorder=4)
    ax.set_title(title,loc='left',fontsize=10.5)
    ax.set_xlabel('GGT activity  $Z$' if k<3 else 'bilirubin  $Z$')
    ax.set_ylabel('expression  $Z$')
    ax.text(0.015,0.965,tag,transform=ax.transAxes,fontsize=8.6,va='top',ha='left',
            fontweight='bold',color=V[vd],zorder=20,
            bbox=dict(boxstyle='round,pad=0.3',fc='white',ec=V[vd],lw=1.0))
    # stat box: bottom-right for positive-slope panel C, else top-right
    if k==2:   # liver: positive slope -> bottom-right clear
        sx,sy,va,ha=0.985,0.045,'bottom','right'
    elif k==3: # UGT1A1: negative slope -> bottom-left clear
        sx,sy,va,ha=0.015,0.045,'bottom','left'
    else:      # A,B negative/shallow -> top-right clear
        sx,sy,va,ha=0.985,0.965,'top','right'
    ax.text(sx,sy,stat,transform=ax.transAxes,fontsize=7.8,va=va,ha=ha,zorder=20,
            bbox=dict(boxstyle='round,pad=0.3',fc='white',ec='#CCCCCC',lw=0.6))

cax=fig.add_axes([0.92,0.30,0.014,0.40]); cb=fig.colorbar(sc,cax=cax)
cb.set_label('$r^2$ to analyte lead variant',fontsize=8.5); cb.ax.tick_params(labelsize=7.5)

fig.savefig('Figure4.pdf',bbox_inches='tight')
fig.savefig('Figure4.png',bbox_inches='tight',dpi=200)
print("Figure 4 done")
