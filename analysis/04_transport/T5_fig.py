import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7})
BL,RD,GY,GN,OR='#2c6fbb','#c0392b','#95a5a6','#27ae60','#e67e22'
T=pd.read_csv("results/T3_transport.csv"); S=pd.read_csv("results/T4_sensitivity.csv")
fig=plt.figure(figsize=(9.4,3.0))
gs=fig.add_gridspec(1,3,width_ratios=[1.35,1,1],wspace=.40,left=.085,right=.985,top=.83,bottom=.22)

ax=fig.add_subplot(gs[0]); x=np.arange(len(T))
ax.plot(x,T.pos_ctrl,'^',color=OR,ms=5,label='positive control (synthetic shell)',zorder=4)
ax.plot(x,T.null_shift,'_',color=GY,ms=11,mew=1.6,label='null: barrier repositioned',zorder=3)
ax.plot(x,T.null_strat,'_',color='#c8ccd0',ms=11,mew=1.6,label='null: barrier permuted in Epi strata',zorder=3)
ax.plot(x,T.detour,'o',color=BL,ms=6,label='observed',zorder=5)
for i in x: ax.plot([i,i],[T.detour[i],T.null_shift[i]],color=BL,lw=.8,alpha=.5,zorder=2)
ax.set_xticks(x); ax.set_xticklabels([f"S{i+1}" for i in x],fontsize=6.5)
ax.set_ylabel('detour ratio\n(least-cost / hop distance)',fontsize=7)
ax.set_title('A   Immune access through the stroma, per section',fontsize=8,loc='left',pad=6)
ax.legend(fontsize=5.5,frameon=False,loc='upper left',ncol=1)
ax.set_ylim(2.2,5.0); ax.spines[['top','right']].set_visible(False)

ax=fig.add_subplot(gs[1])
ax.hist(S.median_z,bins=10,color=BL,edgecolor='w',zorder=2)
ax.axvline(0,color='k',ls='--',lw=.8)
ax.set_xlabel('median z across sections',fontsize=7.2); ax.set_ylabel('parameter settings',fontsize=7.2)
ax.set_title('B   27 threshold combinations',fontsize=8,loc='left',pad=6)
ax.text(.04,.93,f"27/27 settings negative\nall Wilcoxon p=0.001",transform=ax.transAxes,
        fontsize=6.2,va='top')
ax.spines[['top','right']].set_visible(False)

ax=fig.add_subplot(gs[2])
lab=['positive\ncontrol','null\n(repositioned)','null\n(stratified)','observed']
val=[T.pos_ctrl.median(),T.null_shift.median(),T.null_strat.median(),T.detour.median()]
ax.bar(range(4),val,color=[OR,GY,'#c8ccd0',BL],width=.62,zorder=2)
for i,v in enumerate(val): ax.text(i,v+.05,f"{v:.2f}",ha='center',fontsize=6.6)
ax.set_xticks(range(4)); ax.set_xticklabels(lab,fontsize=6.2)
ax.set_ylim(0,4.8); ax.set_ylabel('detour ratio (median)',fontsize=7.2)
ax.set_title('C   More permeable than chance',fontsize=8,loc='left',pad=6)
ax.spines[['top','right']].set_visible(False)
fig.savefig("figures/Fig_A2_transport.pdf"); fig.savefig("figures/Fig_A2_transport.png",dpi=200)
print("done")
