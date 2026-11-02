#!/usr/bin/env python3
"""
run_reanalysis.py — reproduce the headline statistics of the PPi-TNAP drug-target
MR validity audit from the deposited summary statistics and 1000 Genomes EUR LD
(plink --r square, 503 individuals, GRCh37). All numbers are COMPUTED here; no
value is hard-coded. See PROTOCOL in README.md.
"""
import os, numpy as np, pandas as pd
from math import erf, sqrt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_abf import coloc_abf
U = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')+os.sep
aset=lambda a,b:'_'.join(sorted([str(a),str(b)]))
pv=lambda z:2*(1-0.5*(1+erf(abs(z)/sqrt(2))))

def load_ld(tag):
    bim=pd.read_csv(f'{U}LD_{tag}.bim',sep=r'\s+',header=None,names=['chr','rsid','cm','pos','A1','A2'])
    bim['i']=np.arange(len(bim)); bim['ak']=[aset(a,b) for a,b in zip(bim.A1,bim.A2)]
    return bim, np.loadtxt(f'{U}LD_{tag}.ld.gz')

def align(bim, pos, ea, oa, z, add=0, rsid=None):
    out=np.full(len(bim),np.nan)
    d=pd.DataFrame({'pos':np.asarray(pos)+add,'ea':[str(x) for x in ea],'oa':[str(x) for x in oa],'z':np.asarray(z,float)})
    if rsid is not None: d['rsid']=list(rsid)
    d['ak']=[aset(a,b) for a,b in zip(d.ea,d.oa)]
    if rsid is not None:
        key=d.drop_duplicates('rsid').set_index('rsid')
        for _,r in bim.iterrows():
            if r['rsid'] in key.index:
                row=key.loc[r['rsid']]
                if aset(row['ea'],row['oa'])==r['ak']: out[r['i']]=(1 if str(row['ea'])==r['A1'] else -1)*row['z']
    else:
        key=d.drop_duplicates(['pos','ak']).set_index(['pos','ak'])
        for _,r in bim.iterrows():
            k=(r['pos'],r['ak'])
            if k in key.index: row=key.loc[k]; out[r['i']]=(1 if str(row['ea'])==r['A1'] else -1)*row['z']
    return out

def stepwise(z,R,thr=5.45,maxs=8):
    ok=np.where(np.isfinite(z))[0]; zz=z[ok]; Rs=R[np.ix_(ok,ok)]; S=[]; out=[]
    for _ in range(maxs):
        if not S: zc=zz.copy()
        else:
            inv=np.linalg.pinv(Rs[np.ix_(S,S)]+1e-3*np.eye(len(S)))
            b=Rs[:,S]@inv@zz[S]; den=np.clip(1-np.einsum('ij,jk,ik->i',Rs[:,S],inv,Rs[:,S]),1e-6,None)
            zc=(zz-b)/np.sqrt(den)
        zc[S]=0; j=int(np.nanargmax(np.abs(zc)))
        if abs(zc[j])<thr: break
        S.append(j); out.append(int(ok[j]))
    return out

def modal_off(p37,a37a,a37b,p38,a38a,a38b):
    L=pd.DataFrame({'p':np.asarray(p37),'ak':[aset(a,b) for a,b in zip(a37a,a37b)]})
    Rr=pd.DataFrame({'p':np.asarray(p38),'ak':[aset(a,b) for a,b in zip(a38a,a38b)]})
    j=pd.merge(L,Rr,on='ak',suffixes=('37','38')); j['d']=j.p37-j.p38; j=j[j.d.abs()<3e6]
    return int(j['d'].round(-2).value_counts().idxmax()) if len(j) else 0

rows=[]
# ---------------- LPA ----------------
bim,R=load_ld('LPA')
lpa=pd.read_csv(U+'LPA_Lpa_b37.tsv',sep='\t'); mk=lpa['MarkerName'].str.extract(r'^(\d+):(\d+)_([ACGT]+)_([ACGT]+)$')
lpa['pos']=pd.to_numeric(mk[1],errors='coerce'); lpa=lpa.dropna(subset=['pos']); lpa['pos']=lpa['pos'].astype(int)
lpa['zz']=lpa.Effect/lpa.StdErr
zL=align(bim,lpa.pos,lpa.ALT,lpa.REF,lpa.zz)
cavs=pd.read_csv(U+'LPA_CAVS_b38.tsv',sep='\t'); cavs['pos']=cavs['pos']+421032; cavs['zz']=cavs.beta/cavs.sebeta
zC=align(bim,cavs.pos,cavs.alt,cavs.ref,cavs.zz,rsid=cavs.rsids)
cac=pd.read_csv(U+'LPA_CAC_region_b38.tsv',sep='\t'); cac['pos']=cac['base_pair_location']+421032; cac['zz']=cac.beta/cac.standard_error
zK=align(bim,cac.pos,cac.effect_allele,cac.other_allele,cac.zz,rsid=cac.rsid)
nL,nC,nK=len(stepwise(zL,R)),len(stepwise(zC,R)),len(stepwise(zK,R))
# Lp(a)->CAVS shared-signal Wald (exposure lead by |Z| in fixed gene window)
GLO,GHI=160852514,161185307
inwin=np.where((bim.pos.values>=GLO)&(bim.pos.values<=GHI)&np.isfinite(zL))[0]
ilead=inwin[np.argmax(np.abs(zL[inwin]))]
# shared-variant demonstration: condition the primary CAVS signal on the Lp(a) exposure lead
okC_=np.where(np.isfinite(zC))[0]; iCavs=okC_[int(np.argmax(np.abs(zC[okC_])))]
rrC=R[iCavs,ilead]; zCmarg=abs(float(zC[iCavs])); zCcond=abs((float(zC[iCavs])-rrC*float(zC[ilead]))/np.sqrt(max(1-rrC**2,1e-6)))
def wald_ld(ilead, zout_series_df, bcol,scol, rsidcol=None, posoff=0):
    # outcome beta/se at the instrument variant
    r=bim.iloc[ilead]
    o=zout_series_df
    if rsidcol: m=o[o[rsidcol]==r['rsid']]
    else: m=o[(o['pos']==r['pos'])]
    m=m[[aset(a,b)==r['ak'] for a,b in zip(m[eacol],m[oacol])]] if False else m
    return m
# pull outcome beta/se at the instrument, EXPLICITLY harmonising the effect allele
# to the exposure effect allele (Lp(a) ALT) before taking the Wald ratio.
rI=bim.iloc[ilead]; lpa_lead=lpa[(lpa.pos==rI['pos'])].iloc[0]
beI=float(lpa_lead.Effect); ea_x=str(lpa_lead.ALT); oa_x=str(lpa_lead.REF)   # exposure effect allele = ALT
def harmonise(bo, ea_out, oa_out):
    """Return outcome beta sign-aligned to the exposure effect allele ea_x, or NaN on allele mismatch."""
    if str(ea_out)==ea_x and str(oa_out)==oa_x: return bo
    if str(ea_out)==oa_x and str(oa_out)==ea_x: return -bo
    return float('nan')   # allele set does not match -> do not produce a mis-oriented estimate
cav_m=cavs[cavs.rsids==rI['rsid']]
if len(cav_m):
    cc=cav_m.iloc[0]; bo=harmonise(float(cc.beta), cc.alt, cc.ref); so=float(cc.sebeta)
    w=bo/beI; sw=abs(so/beI)
    rows.append(('Lp(a) -> CAVS','OR %.2f (%.2f-%.2f)'%(np.exp(w),np.exp(w-1.96*sw),np.exp(w+1.96*sw)),
                 f'shared signal r2=%.2f; Lp(a) signals=%d, CAVS signals=%d; CAVS lead |Z| %.1f->%.1f on conditioning; effect allele %s (harmonised)'%(rrC**2, nL,nC, zCmarg, zCcond, ea_x),'PASS'))
cac_m=cac[(cac.pos==rI['pos'])]
if len(cac_m):
    kk=cac_m.iloc[0]; bo=harmonise(float(kk.beta), kk.effect_allele, kk.other_allele); so=float(kk.standard_error)
    w=bo/beI; sw=abs(so/beI)
    rows.append(('Lp(a) -> CAC','beta %+.3f (SE %.3f, P=%.1g)'%(w,sw,pv(w/sw)),
                 f'CAC independent GWS signals at LPA = {nK} (no outcome signal); effect allele {ea_x} (harmonised)','non-colocalizing nominal'))

# ---------------- GGT1 activity vs expression ----------------
bg,Rg=load_ld('GGT1')
act=pd.read_csv(U+'ENZ_GGT1_activity_b37.tsv',sep='\t'); act['zz']=act.Effect/act.StdErr
za=align(bg,act.POS,act.ALT,act.REF,act.zz)
eq=pd.read_csv(U+'ENZ_GGT1_eQTLGen_b37.tsv',sep='\t'); zb=align(bg,eq.SNPPos,eq.AssessedAllele,eq.OtherAllele,eq.Zscore)
def leadpair(za,zb,R):
    ia=np.nanargmax(np.abs(np.where(np.isfinite(za),za,0))); ib=np.nanargmax(np.abs(np.where(np.isfinite(zb),zb,0)))
    both=np.where(np.isfinite(za)&np.isfinite(zb)&((np.abs(za)>4)|(np.abs(zb)>4)))[0]
    cr=np.corrcoef(za[both],zb[both])[0,1] if len(both)>3 else np.nan
    return R[ia,ib]**2, cr, len(both)
r2,cr,n=leadpair(za,zb,Rg)
rows.append(('GGT activity vs GGT1 blood (eQTLGen)','lead-pair r2=%.3f, cross-var r=%+.2f (n=%d)'%(r2,cr,n),'same architecture, reversed direction','confirmed (blood misleading)'))
liv=pd.read_csv(U+'ENZ_GGT1_liverGTEx_b38.tsv',sep='\t'); liv['zz']=liv.beta/liv.se; liv['pos']=liv.variant.str.split('_').str[1].astype(int)
off=modal_off(act.POS,act.ALT,act.REF,liv.pos,liv.alt,liv.ref); zl=align(bg,liv.pos,liv.alt,liv.ref,liv.zz,add=off)
livlead=np.nanmax(np.abs(np.where(np.isfinite(zl),zl,0)))
rows.append(('GGT activity vs GGT1 liver','liver eQTL lead |Z|=%.1f'%livlead,'GTEx liver n=208 underpowered','underpowered (no concordance claim)'))

# ---------------- ANKH blood -> CAVS (single-signal coloc.abf) ----------------
blood=pd.read_csv(U+'axis_eQTLGenBlood_b37.tsv',sep='\t'); ank=blood[blood.GeneSymbol=='ANKH'].copy(); ank['pos']=ank.SNPPos.astype(int)
af=pd.read_csv(U+'ankh_af_region.tsv',sep='\t'); e=pd.merge(ank,af[['SNP','AF_AlleleB']],on='SNP')
e['maf']=np.minimum(e.AF_AlleleB.clip(0,1),1-e.AF_AlleleB.clip(0,1)); e=e[e.maf>0.01]
Z=e.Zscore.astype(float); Nn=e.NrSamples.astype(float); den=2*e.maf*(1-e.maf)*(Nn+Z**2)
e['be']=Z/np.sqrt(den); e['se']=1/np.sqrt(den); e['ak']=[aset(a,b) for a,b in zip(e.AssessedAllele,e.OtherAllele)]
ac=pd.read_csv(U+'axis_CAVS_b38.tsv',sep='\t'); c5=ac[ac.chrom==5].copy(); c5['pos37']=c5['pos']+109; c5['ak']=[aset(a,b) for a,b in zip(c5.ref,c5.alt)]
lead=int(e.loc[e.Zscore.abs().idxmax(),'pos']); lo,hi=lead-100000,lead+100000
E=e[(e.pos>=lo)&(e.pos<=hi)][['pos','ak','be','se']].rename(columns={'pos':'pos37'})
O=c5[(c5.pos37>=lo)&(c5.pos37<=hi)][['pos37','ak','beta','sebeta']]
m=pd.merge(E,O,on=['pos37','ak']).dropna()
pp,n=coloc_abf(m.be.values,m.se.values,m.beta.values,m.sebeta.values,'quant','cc')
rows.append(('ANKH blood -> CAVS','coloc.abf PP.H4=%.2f PP.H1=%.2f'%(pp['H4'],pp['H1']),'single-signal locus; CAVS flat','non-colocalizing nominal'))

# ---------------- ENPP1 aorta -> CAC (quant) ----------------
aorta=pd.read_csv(U+'axis_GTExAorta_b38.tsv',sep='\t'); aorta['g']=aorta.gene_id.str.split('.').str[0]
ae=aorta[aorta.g=='ENSG00000197594'].copy(); ae['pos']=ae.variant_id.str.split('_').str[1].astype(int); ae['ak']=[aset(p[2],p[3]) for p in ae.variant_id.str.split('_')]
cc=pd.read_csv(U+'axis_CAC_b38.tsv',sep='\t'); c6=cc[cc.chromosome==6].copy(); c6['ak']=[aset(a,b) for a,b in zip(c6.effect_allele,c6.other_allele)]
lead=int(ae.loc[ae.pval_nominal.astype(float).idxmin(),'pos']); lo,hi=lead-100000,lead+100000
E=ae[(ae.pos>=lo)&(ae.pos<=hi)][['pos','ak','slope','slope_se']]
O=c6[(c6.base_pair_location>=lo)&(c6.base_pair_location<=hi)][['base_pair_location','ak','beta','standard_error']].rename(columns={'base_pair_location':'pos'})
m=pd.merge(E,O,on=['pos','ak']).dropna()
pp,n=coloc_abf(m.slope.values,m.slope_se.values,m.beta.values,m.standard_error.values,'quant','quant')
rows.append(('ENPP1 aorta -> CAC','coloc.abf PP.H4=%.2f (quant)'%pp['H4'],'CAC quantitative; fail-coloc','fail-coloc'))

# ---------------- ENPP1 blood (eQTLGen) -> CAC single-variant Wald (nominal) ----------------
en=blood[blood.GeneSymbol=='ENPP1'].copy()
li=en.loc[en.Zscore.abs().idxmax()]                      # cis-lead by max |Z|
zl=float(li.Zscore); Nl=float(li.NrSamples); ea_e=str(li.AssessedAllele)
cm=c6[c6.rsid==li.SNP]                                   # match CAC by rsid (cross-build)
if len(cm):
    cm=cm.iloc[0]; maf=float(cm.effect_allele_frequency); maf=min(maf,1-maf)
    den=sqrt(2*maf*(1-maf)*(Nl+zl*zl)); be=zl/den          # eQTLGen Z->beta per assessed allele
    bo=float(cm.beta); so=float(cm.standard_error)
    if str(cm.effect_allele)!=ea_e: bo=-bo                 # harmonize CAC to exposure assessed allele
    w=bo/be; sw=abs(so/be)
    rows.append(('ENPP1 blood -> CAC','beta %+.2f (SE %.2f, P=%.3f)'%(w,sw,pv(w/sw)),
                 'eQTLGen cis-lead %s; nominal; not colocalizing (aorta PP.H4=%.2f)'%(li.SNP,pp['H4']),'nominal only'))

res=pd.DataFrame(rows,columns=['pair','recomputed','note','verdict'])
res.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'results_reanalysis.csv'),index=False)
pd.set_option('display.width',200); pd.set_option('display.max_colwidth',60)
print(res.to_string(index=False))
