import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7})
BL,RD,GY,GN='#2c6fbb','#c0392b','#95a5a6','#27ae60'
fig=plt.figure(figsize=(9.4,3.0))
gs=fig.add_gridspec(1,3,width_ratios=[1.15,1,1.05],wspace=.42,left=.105,right=.985,top=.84,bottom=.24)

# A: upgrades attempted (prognosis, LOCO C-index)
ax=fig.add_subplot(gs[0])
lab=['EN: abundance\n+ state (30)','EN: abundance\nonly (10)','NetImmune\n(fixed, 48 genes)','ESTIMATE\n(fixed, 282)','EN: 6 barriernet\nmodules','EN: state\nonly (20)','CD8A\n(1 gene)']
val=[0.598,0.596,0.589,0.589,0.585,0.523,0.517]
col=[BL,BL,RD,GY,BL,BL,GY]
ax.barh(range(len(val))[::-1],val,color=col,height=.66,zorder=2)
ax.axvline(0.589,color=RD,ls='--',lw=.8,zorder=3)
ax.set_yticks(range(len(val))[::-1]); ax.set_yticklabels(lab,fontsize=5.8)
ax.set_xlim(.45,.63); ax.set_xlabel('Out-of-cohort C-index (LOCO)',fontsize=7.2)
ax.set_title('A   Three upgrades vs the current score',fontsize=8,loc='left',pad=6)
for i,v in enumerate(val): ax.text(v+.003,len(val)-1-i,f"{v:.3f}",va='center',fontsize=5.8)
ax.spines[['top','right']].set_visible(False)

# B: spatial arrangement transfer
ax=fig.add_subplot(gs[1])
g=['learnable\nin spatial data','bulk prognosis\n(ΔC-index)','bulk ICI\n(AUC-0.5)','negative control\n(shuffled)']
v=[0.219,0.006,0.003,-0.003]; c=[GN,GY,GY,'#bdc3c7']
ax.bar(range(4),v,color=c,width=.62,zorder=2)
ax.axhline(0,color='k',lw=.6)
for i,x in enumerate(v): ax.text(i,x+.008 if x>0 else x-.02,f"{x:+.3f}",ha='center',fontsize=6.2)
ax.set_xticks(range(4)); ax.set_xticklabels(g,fontsize=5.8)
ax.set_ylim(-.04,.27); ax.set_ylabel('effect size',fontsize=7.2)
ax.set_title('B   Spatial arrangement transfer',fontsize=8,loc='left',pad=6)
ax.spines[['top','right']].set_visible(False)

# C: the ceiling
ax=fig.add_subplot(gs[2])
x=np.arange(2); w=.34
full=[0.600,0.754]; simp=[0.589,0.854]
ax.bar(x-w/2,full,w,color='#34495e',label='full transcriptome\n+ learned weights',zorder=2)
ax.bar(x+w/2,simp,w,color=RD,label='fixed hand-curated\nsignature',zorder=2)
for xi,(a,b_) in zip(x,zip(full,simp)):
    ax.text(xi-w/2,a+.008,f"{a:.3f}",ha='center',fontsize=6.4)
    ax.text(xi+w/2,b_+.008,f"{b_:.3f}",ha='center',fontsize=6.4,color=RD,fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(['Prognosis (LOCO C-index)\n14,386 vs 48 genes',
                                      'ICI response (CV AUC)\n31,238 vs 5 genes'],fontsize=6.2)
ax.set_ylim(.5,.92); ax.set_ylabel('performance',fontsize=7.2)
ax.set_title('C   The information ceiling',fontsize=8,loc='left',pad=6)
ax.legend(fontsize=5.6,frameon=False,loc='upper left')
ax.spines[['top','right']].set_visible(False)
fig.savefig("figures/Fig_A1_ceiling.pdf"); fig.savefig("figures/Fig_A1_ceiling.png",dpi=200)
print("done")
