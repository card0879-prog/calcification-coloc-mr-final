import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
from matplotlib.gridspec import GridSpec
from figstyle import V, EXPO, OUT, GWL, badge, note, LABEL

U='data/'
def llp(p): p=pd.to_numeric(p,errors='coerce'); return -np.log10(np.clip(p,1e-320,1))

# ---------- load ----------
lpa=pd.read_csv(U+'LPA_Lpa_b37.tsv',sep='\t'); lpa['llp']=llp(lpa['P-value'])
cavs=pd.read_csv(U+'LPA_CAVS_b38.tsv',sep='\t'); cavs['llp']=llp(cavs['pval']); cavs['b37']=cavs['pos']+421032
blood=pd.read_csv(U+'axis_eQTLGenBlood_b37.tsv',sep='\t'); blood['llp']=llp(blood['Pvalue'])
ankh=blood[blood.GeneSymbol=='ANKH'].copy()
acavs=pd.read_csv(U+'axis_CAVS_b38.tsv',sep='\t'); acavs['llp']=llp(acavs['pval'])
acavs5=acavs[acavs.chrom==5].copy(); acavs5['b37']=acavs5['pos']+109
acac=pd.read_csv(U+'axis_CAC_b38.tsv',sep='\t'); acac['llp']=llp(acac['p_value']); acac6=acac[acac.chromosome==6].copy()
aorta=pd.read_csv(U+'axis_GTExAorta_b38.tsv',sep='\t'); aorta['llp']=llp(aorta['pval_nominal'])
aorta['gene']=aorta['gene_id'].str.split('.').str[0]
ae=aorta[aorta.gene=='ENSG00000197594'].copy(); ae['pos']=ae['variant_id'].str.split('_').str[1].astype(int)

def win(df,pc,lo,hi): return df[(df[pc]>=lo)&(df[pc]<=hi)]

# ---------- figure ----------
fig=plt.figure(figsize=(11.0,8.2))
gs=GridSpec(2,2,figure=fig,hspace=0.46,wspace=0.44,left=0.075,right=0.93,top=0.88,bottom=0.085)
fig.suptitle('The gate admits benchmark colocalizing signals but rejects nominal, non-colocalizing ones',
             fontsize=13.5,fontweight='bold',y=0.965)

def overlay(ax, ex_pos, ex_llp, ou_pos, ou_llp, exlab, oulab, ex_ymax, ou_ymax,
            title, verdict, txt, gwl=True):
    ax.scatter(ex_pos/1e6, ex_llp, s=9, c=EXPO, alpha=0.55, edgecolors='none', label=exlab, zorder=3)
    ax.set_ylim(0,ex_ymax); ax.set_ylabel('exposure  $-\\log_{10}P$', color=EXPO)
    ax.tick_params(axis='y', colors=EXPO)
    # lead star (exposure)
    i=np.argmax(ex_llp); ax.scatter([ex_pos.iloc[i]/1e6],[ex_llp.iloc[i]],marker='*',s=130,
            c=EXPO,edgecolors='white',linewidths=0.6,zorder=5)
    ax2=ax.twinx()
    ax2.scatter(ou_pos/1e6, ou_llp, s=9, c=OUT, alpha=0.5, edgecolors='none', label=oulab, zorder=2)
    ax2.set_ylim(0,ou_ymax); ax2.set_ylabel('outcome  $-\\log_{10}P$', color=OUT)
    ax2.tick_params(axis='y', colors=OUT)
    if gwl: ax2.axhline(7.3,ls='--',lw=0.8,c=GWL,alpha=0.6,zorder=1)
    ax.set_title(title, loc='left')
    ax.set_xlabel('chromosomal position (Mb, GRCh37)')
    badge(ax,0.015,0.97,verdict)
    note(ax,0.015,0.82,txt)
    # combined legend
    h1=plt.Line2D([],[],marker='o',ls='',c=EXPO,label=exlab,ms=5)
    h2=plt.Line2D([],[],marker='o',ls='',c=OUT,label=oulab,ms=5)
    ax.legend(handles=[h1,h2],loc='upper right',frameon=True,framealpha=0.9,borderpad=0.4)
    return ax2

# A) LPA -> CAVS  (admit)
axA=fig.add_subplot(gs[0,0])
L=win(lpa,'POS',159.8e6,161.4e6); C=win(cavs,'b37',159.8e6,161.4e6)
overlay(axA, L['POS'], L['llp'], C['b37'], C['llp'],
        'Lp(a) activity','CAVS', 345, 42,
        'A   Lp(a) \u2192 CAVS  (LPA locus)','PASS_strong',
        'shared causal variant (conditional)\nlead-pair $r^2$ = 1.00; cond. $Z$ 12.4 \u2192 0.0\n(multi-signal LPA locus)')

# B) ANKH blood -> CAVS  (reject: H1-dominant)
axB=fig.add_subplot(gs[0,1])
A=win(ankh,'SNPPos',14.6e6,15.1e6); Cc=win(acavs5,'b37',14.6e6,15.1e6)
overlay(axB, A['SNPPos'], A['llp'], Cc['b37'], Cc['llp'],
        'ANKH blood expression','CAVS', 330, 8.5,
        'B   ANKH (blood) \u2192 CAVS','noncoloc',
        'PP.H4 = 0.27   PP.H1 = 0.71\noutcome: no regional signal (max $-\\log_{10}P$ = 3.2)')

# C) ENPP1 aorta -> CAC  (reject: fail-coloc)
axC=fig.add_subplot(gs[1,0])
E=win(ae,'pos',131.6e6,132.0e6); Ca=win(acac6,'base_pair_location',131.6e6,132.0e6)
ax2=overlay(axC, E['pos'], E['llp'], Ca['base_pair_location'], Ca['llp'],
        'ENPP1 aortic expression','CAC', 8.5, 8.5,
        'C   ENPP1 (aorta) \u2192 CAC','fail_coloc',
        'PP.H4 = 0.05   instrumented, no shared signal\n(MR $\\beta$ = \u22120.09, $P$ = 0.49)')
axC.set_xlabel('chromosomal position (Mb, GRCh38)')

# D) PP.H4 decision panel
axD=fig.add_subplot(gs[1,1])
tests=[
 ('Lp(a) \u2192 CAVS',0.998,'PASS_strong',''),
 ('Lp(a) \u2192 CAC',0.05,'noncoloc','no CAC signal'),
 ('ALP activity \u2192 CAVS',0.13,'fail_coloc',''),
 ('ENPP1 aorta \u2192 CAC',0.05,'fail_coloc',''),
 ('ANKH blood \u2192 CAVS',0.27,'noncoloc','H1=0.71'),
 ('ENPP1 blood \u2192 CAC',0.04,'noncoloc','H1=0.95'),
]
ypos=np.arange(len(tests))[::-1]
for y,(name,h4,vd,ann) in zip(ypos,tests):
    axD.barh(y,h4,color=V[vd],edgecolor='white',height=0.62,zorder=3)
    axD.text(h4+0.015,y,f'{h4:.3f}'.rstrip('0').rstrip('.') if h4<1 else f'{h4}',
             va='center',ha='left',fontsize=8.3)
    if ann: axD.text(h4+0.10,y,ann,va='center',ha='left',fontsize=7.6,color=V['noncoloc'],style='italic')
axD.axvline(0.8,ls='--',lw=0.9,c='#333333',alpha=0.7); axD.text(0.79,2.6,'strong \u22650.8',rotation=90,va='center',ha='right',fontsize=7.2,color='#333333')
axD.axvline(0.5,ls=':',lw=0.9,c='#666666',alpha=0.7); axD.text(0.49,2.6,'moderate \u22650.5',rotation=90,va='center',ha='right',fontsize=7.2,color='#666666')
axD.set_yticks(ypos); axD.set_yticklabels([t[0] for t in tests],fontsize=8.0)
axD.set_xlim(0,1.18); axD.set_xlabel('colocalization posterior  PP.H4')
axD.set_title('D   The decision variable: PP.H4 across tests',loc='left',pad=10)
axD.spines['right'].set_visible(False); axD.spines['top'].set_visible(False)

fig.savefig('Figure3.pdf',bbox_inches='tight')
fig.savefig('Figure3.png',bbox_inches='tight',dpi=200)
print("Figure 3 done")
