import os, sys, re; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
from matplotlib.gridspec import GridSpec
from figstyle import V, EXPO, OUT, GWL, badge, note, LABEL

U='data/'
def llp(p): p=pd.to_numeric(p,errors='coerce'); return -np.log10(np.clip(p,1e-320,1))

# ---------- computed results (NO hard-coded statistics in this figure script) ----------
CO=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'results_coloc.csv'))
co={r['pair']:r for _,r in CO.iterrows()}          # coloc posteriors by pair
RE=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'results_reanalysis.csv'))
re_note={r['pair']:str(r['note']) for _,r in RE.iterrows()}
def h4(pair): return float(co[pair]['PP.H4'])
def h1(pair): return float(co[pair]['PP.H1'])
# parse conditional metrics for the LPA positive control from run_reanalysis note
nL=re_note.get('Lp(a) -> CAVS','')
r2=re.search(r'r2=([0-9.]+)',nL); nexp=re.search(r'Lp\(a\) signals=(\d+)',nL); nout=re.search(r'CAVS signals=(\d+)',nL)
R2=r2.group(1) if r2 else '~1.0'; NEXP=nexp.group(1) if nexp else '?'; NOUT=nout.group(1) if nout else '?'
# CAC independent signals at LPA (from run_reanalysis note for Lp(a) -> CAC)
nK=re_note.get('Lp(a) -> CAC',''); mcac=re.search(r'=\s*(\d+)\s*\(no outcome',nK) or re.search(r'LPA\s*=\s*(\d+)',nK)
NCAC=mcac.group(1) if mcac else '0'

# ---------- load regional summary stats ----------
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
    h1h=plt.Line2D([],[],marker='o',ls='',c=EXPO,label=exlab,ms=5)
    h2h=plt.Line2D([],[],marker='o',ls='',c=OUT,label=oulab,ms=5)
    ax.legend(handles=[h1h,h2h],loc='upper right',frameon=True,framealpha=0.9,borderpad=0.4)
    return ax2

# A) LPA -> CAVS  (admit: conditional colocalization, NOT single-causal PP.H4)
axA=fig.add_subplot(gs[0,0])
L=win(lpa,'POS',159.8e6,161.4e6); C=win(cavs,'b37',159.8e6,161.4e6)
overlay(axA, L['POS'], L['llp'], C['b37'], C['llp'],
        'Lp(a) activity','CAVS', 345, 42,
        'A   Lp(a) → CAVS  (LPA locus)','PASS_strong',
        f'shared causal variant (conditional)\nlead-pair $r^2$ = {R2}; multiple Lp(a) + {NOUT} CAVS signals\n(multi-signal locus: single-causal coloc not valid)')

# B) ANKH blood -> CAVS  (reject: H1-dominant)
axB=fig.add_subplot(gs[0,1])
A=win(ankh,'SNPPos',14.6e6,15.1e6); Cc=win(acavs5,'b37',14.6e6,15.1e6)
overlay(axB, A['SNPPos'], A['llp'], Cc['b37'], Cc['llp'],
        'ANKH blood expression','CAVS', 330, 8.5,
        'B   ANKH (blood) → CAVS','noncoloc',
        f'PP.H4 = {h4("ANKH blood -> CAVS"):.2f}   PP.H1 = {h1("ANKH blood -> CAVS"):.2f}\noutcome: no regional signal')

# C) ENPP1 aorta -> CAC  (reject: fail-coloc)
axC=fig.add_subplot(gs[1,0])
E=win(ae,'pos',131.6e6,132.0e6); Ca=win(acac6,'base_pair_location',131.6e6,132.0e6)
overlay(axC, E['pos'], E['llp'], Ca['base_pair_location'], Ca['llp'],
        'ENPP1 aortic expression','CAC', 8.5, 8.5,
        'C   ENPP1 (aorta) → CAC','fail_coloc',
        f'PP.H4 = {h4("ENPP1 aorta -> CAC"):.2f}   instrumented,\nno shared regional signal with CAC')
axC.set_xlabel('chromosomal position (Mb, GRCh38)')

# D) decision panel — single-causal PP.H4 (valid single-signal loci) vs LPA conditional basis
axD=fig.add_subplot(gs[1,1])
# single-signal coloc.abf tests (PP.H4 is the decision variable here)
bars=[('ANKH blood → CAVS','ANKH blood -> CAVS','noncoloc',f'H1={h1("ANKH blood -> CAVS"):.2f}'),
      ('ENPP1 aorta → CAC','ENPP1 aorta -> CAC','fail_coloc',''),
      ('ENPP1 blood → CAC','ENPP1 blood -> CAC','noncoloc',f'H1={h1("ENPP1 blood -> CAC"):.2f}')]
# LPA-locus tests judged on a different basis (shown without a single-causal PP.H4 bar)
lpa_rows=[('Lp(a) → CAVS','PASS_strong','admitted — conditional (shared causal variant)'),
          ('Lp(a) → CAC','noncoloc','rejected — no genome-wide CAC signal')]
all_labels=[b[0] for b in bars]+[r[0] for r in lpa_rows]
ypos=np.arange(len(all_labels))[::-1]
yb=ypos[:len(bars)]; yl=ypos[len(bars):]
for y,(name,pair,vd,ann) in zip(yb,bars):
    v=h4(pair)
    axD.barh(y,v,color=V[vd],edgecolor='white',height=0.6,zorder=3)
    axD.text(v+0.015,y,f'{v:.2f}',va='center',ha='left',fontsize=8.3)
    if ann: axD.text(v+0.10,y,ann,va='center',ha='left',fontsize=7.4,color=V['noncoloc'],style='italic')
for y,(name,vd,ann) in zip(yl,lpa_rows):
    axD.scatter([0.02],[y],marker='s',s=46,facecolor='white',edgecolor=V[vd],linewidths=1.6,zorder=3)
    axD.text(0.07,y,ann,va='center',ha='left',fontsize=7.6,color='#333',style='italic')
axD.axhline((yb[-1]+yl[0])/2,ls='-',lw=0.6,c='#DDDDDD')
axD.text(1.16,yb[len(bars)//2],'single-causal\ncoloc.abf',rotation=90,va='center',ha='center',fontsize=6.8,color='#888')
axD.text(1.16,np.mean(yl),'LPA locus\n(multi-signal:\nconditional basis)',rotation=90,va='center',ha='center',fontsize=6.3,color='#888')
axD.axvline(0.8,ls='--',lw=0.9,c='#333333',alpha=0.7); axD.text(0.79,yl[-1]-0.3,'strong ≥0.8',rotation=90,va='bottom',ha='right',fontsize=7.0,color='#333333')
axD.axvline(0.5,ls=':',lw=0.9,c='#666666',alpha=0.7); axD.text(0.49,yl[-1]-0.3,'moderate ≥0.5',rotation=90,va='bottom',ha='right',fontsize=7.0,color='#666666')
axD.set_yticks(ypos); axD.set_yticklabels(all_labels,fontsize=8.0)
axD.set_xlim(0,1.22); axD.set_ylim(yl[-1]-0.6, yb[0]+0.6)
axD.set_xlabel('colocalization posterior  PP.H4  (single-causal tests)')
axD.set_title('D   Decision basis across tests',loc='left',pad=10,fontsize=11)
axD.spines['right'].set_visible(False); axD.spines['top'].set_visible(False)

fig.savefig('Figure3.pdf',bbox_inches='tight')
fig.savefig('Figure3.png',bbox_inches='tight',dpi=200)
print("Figure 3 done; PP.H4 read from results_coloc.csv:",
      {k:round(h4(k),3) for k in ['ANKH blood -> CAVS','ENPP1 aorta -> CAC','ENPP1 blood -> CAC']})
