import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
matplotlib.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8,
                            'axes.linewidth':.7,'xtick.major.width':.7,'ytick.major.width':.7})
R="/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
F="/groups/adv2105_gp/yichen/Yi/com/gastric/revision/figures"
BLUE,RED,GREY,GREEN='#2c6fbb','#c0392b','#7f8c8d','#27ae60'

# ---------------- Figure R1: prognostic benchmarking ----------------
meta=pd.read_csv(f"{R}/B4_meta.csv"); cidx=pd.read_csv(f"{R}/B4_cindex.csv")
order=['ESTIMATE(imm-str)','NetImmune','NetImmune(ssGSEA)','Immune','CD8A','Barrier']
meta=meta.set_index('score').loc[order].reset_index()
fig=plt.figure(figsize=(9.0,3.1))
gs=fig.add_gridspec(1,2,width_ratios=[1.9,1],wspace=.30,left=.13,right=.97,top=.86,bottom=.20)

ax=fig.add_subplot(gs[0]); y=np.arange(len(meta))[::-1]
for i,(_,r) in enumerate(meta.iterrows()):
    c=RED if r.HR>1 else (BLUE if r['score']=='NetImmune' else GREY)
    lw=1.6 if r['score']=='NetImmune' else 1.0
    ax.plot([r.lo,r.hi],[y[i]]*2,color=c,lw=lw,zorder=2)
    ax.plot(r.HR,y[i],'o',color=c,ms=6 if r['score']=='NetImmune' else 4.5,zorder=3)
    ax.text(1.50,y[i],f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})   p={r.p:.1e}",va='center',fontsize=6.4,clip_on=False)
ax.axvline(1,color='k',ls='--',lw=.6,zorder=1)
ax.set_yticks(y); ax.set_yticklabels(meta['score'],fontsize=7.2)
ax.set_xlim(.55,2.45); ax.set_xticks([.6,.8,1.0,1.2,1.4]); ax.set_xlabel('Pooled HR per 1 SD (95% CI)',fontsize=7.5)
ax.set_title('A   Prognostic meta-analysis, 3 cohorts (n=1,097)',fontsize=8,loc='left',pad=6)
ax.spines[['top','right']].set_visible(False)
ax.spines['bottom'].set_bounds(.55,1.45)
for lbl in ax.get_yticklabels():
    if lbl.get_text()=='NetImmune': lbl.set_fontweight('bold')

ax=fig.add_subplot(gs[1])
sub=cidx[cidx.score.isin(['NetImmune','ESTIMATE(imm-str)','CD8A'])]
coh=['TCGA','GSE84437','GSE62254']; scs=['NetImmune','ESTIMATE(imm-str)','CD8A']
w=.26; cols={'NetImmune':BLUE,'ESTIMATE(imm-str)':GREEN,'CD8A':GREY}
for j,s in enumerate(scs):
    d=sub[sub.score==s].set_index('cohort').loc[coh]
    x=np.arange(3)+(j-1)*w
    ax.bar(x,d.dC,w,color=cols[s],label=s,zorder=2)
    ax.errorbar(x,d.dC,yerr=[d.dC-d.dC_lo,d.dC_hi-d.dC],fmt='none',ecolor='k',lw=.6,capsize=1.8,zorder=3)
    for xi,(_,r) in zip(x,d.iterrows()):
        if r.p_boot<.05: ax.text(xi,r.dC_hi+.002,'*',ha='center',fontsize=9)
ax.axhline(0,color='k',lw=.6)
ax.set_xticks(range(3)); ax.set_xticklabels(coh,fontsize=7)
ax.set_ylabel('ΔC-index vs age + stage',fontsize=7.5); ax.set_ylim(-.020,.058)
ax.set_title('B   Incremental discrimination',fontsize=8,loc='left',pad=6)
ax.legend(fontsize=5.8,frameon=False,loc='lower right',bbox_to_anchor=(1.0,.02))
ax.spines[['top','right']].set_visible(False)
fig.savefig(f"{F}/Fig_R1_prognostic_benchmark.pdf"); fig.savefig(f"{F}/Fig_R1_prognostic_benchmark.png",dpi=200)
print("Fig R1 done")

# ---------------- Figure R2: ICI benchmarking ----------------
auc=pd.read_csv(f"{R}/A2_auc_ci.csv"); dl=pd.read_csv(f"{R}/A2_delong.csv"); cv=pd.read_csv(f"{R}/A3_cv_combo.csv")
q=dict(zip(dl.Comparator,dl.p_FDR))
fig=plt.figure(figsize=(7.2,3.3))
gs=fig.add_gridspec(1,2,width_ratios=[1.35,1],wspace=.45,left=.19,right=.97,top=.86,bottom=.17)
ax=fig.add_subplot(gs[0]); auc=auc.sort_values('AUC'); y=np.arange(len(auc))
for i,(_,r) in enumerate(auc.iterrows()):
    isn=r.Marker=='NetImmune'; c=BLUE if isn else GREY
    ax.plot([r.CI_low,r.CI_high],[y[i]]*2,color=c,lw=1.6 if isn else 1.0,zorder=2)
    ax.plot(r.AUC,y[i],'o',color=c,ms=6 if isn else 4.5,zorder=3)
    lab=f"{r.AUC:.3f}"
    if r.Marker in q and q[r.Marker]<.05: lab+=f"  q={q[r.Marker]:.3f}"
    ax.text(1.035,y[i],lab,va='center',fontsize=6.3,
            color=RED if (r.Marker in q and q[r.Marker]<.05) else 'k')
ax.axvline(.5,color='k',ls='--',lw=.6)
ax.set_yticks(y); ax.set_yticklabels(auc.Marker,fontsize=7)
for lbl in ax.get_yticklabels():
    if lbl.get_text()=='NetImmune': lbl.set_fontweight('bold')
ax.set_xlim(.25,1.0); ax.set_xlabel('AUC (95% bootstrap CI)',fontsize=7.5)
ax.set_title('A   Pembrolizumab response, n=44 (12 responders)',fontsize=8,loc='left',pad=6)
ax.spines[['top','right']].set_visible(False)

ax=fig.add_subplot(gs[1])
x=np.arange(len(cv)); w=.36
ax.bar(x-w/2,cv.Apparent,w,color='#d5d8dc',edgecolor='k',lw=.5,label='Apparent (in-sample)',zorder=2)
ax.bar(x+w/2,cv.CV5_mean,w,color=BLUE,edgecolor='k',lw=.5,label='5-fold CV (200 reps)',zorder=2)
ax.errorbar(x+w/2,cv.CV5_mean,yerr=[cv.CV5_mean-cv.CV5_lo,cv.CV5_hi-cv.CV5_mean],
            fmt='none',ecolor='k',lw=.6,capsize=1.8,zorder=3)
for xi,(_,r) in zip(x,cv.iterrows()):
    ax.text(xi-w/2,r.Apparent+.012,f"{r.Apparent:.3f}",ha='center',fontsize=5.8)
    ax.text(xi+w/2,r.CV5_hi+.012,f"{r.CV5_mean:.3f}",ha='center',fontsize=5.8,color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(['MSI/EBV','NetImmune','NetImmune\n+ MSI/EBV'],fontsize=6.5)
ax.set_ylim(.5,1.06); ax.set_ylabel('AUC',fontsize=7.5)
ax.set_title('B   Optimism in the combined model',fontsize=8,loc='left',pad=6)
ax.legend(fontsize=5.8,frameon=False,loc='lower right')
ax.spines[['top','right']].set_visible(False)
fig.savefig(f"{F}/Fig_R2_ICI_benchmark.pdf"); fig.savefig(f"{F}/Fig_R2_ICI_benchmark.png",dpi=200)
print("Fig R2 done")

# ---------------- Figure R3: GDSC with FDR ----------------
g=pd.read_csv(f"{R}/C1_gdsc_fdr.csv").sort_values('Barrier_r')
fig,ax=plt.subplots(figsize=(4.3,3.1)); plt.subplots_adjust(left=.29,right=.96,top=.84,bottom=.17)
y=np.arange(len(g))
cols=[RED if (r.Barrier_p<.05 and r.Barrier_r>0) else (BLUE if r.Barrier_p<.05 else GREY) for _,r in g.iterrows()]
ax.barh(y,g.Barrier_r,color=cols,height=.68,zorder=2)
for i,(_,r) in enumerate(g.iterrows()):
    if r.Barrier_p<.05:
        ax.text(0.03 if r.Barrier_r<0 else r.Barrier_r+0.02, y[i],
                f"p={r.Barrier_p:.3f}, q={r.Barrier_q:.2f}",va='center',
                ha='left',fontsize=5.8)
ax.axvline(0,color='k',lw=.6)
ax.set_yticks(y); ax.set_yticklabels(g.Drug,fontsize=6.8)
ax.set_xlim(-.62,.72); ax.set_xlabel('Spearman r (Barrier vs ln IC50), 28 GC lines',fontsize=7)
ax.set_title('GDSC2: no compound survives FDR (best q=0.24)',fontsize=7.6,loc='left',pad=6)
ax.text(-.60,len(g)-.15,'← more sensitive',fontsize=6,color=BLUE)
ax.text(.70,len(g)-.15,'more resistant →',fontsize=6,color=RED,ha='right')
ax.set_ylim(-.7,len(g)+.3)
ax.spines[['top','right']].set_visible(False)
fig.savefig(f"{F}/Fig_R3_GDSC_fdr.pdf"); fig.savefig(f"{F}/Fig_R3_GDSC_fdr.png",dpi=200)
print("Fig R3 done")
