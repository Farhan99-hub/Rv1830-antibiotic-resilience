#!/usr/bin/env python3
"""Figure 1 — FINAL, 40,276 analysis set (author-verified values)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":13,"axes.linewidth":1.1,
                     "axes.edgecolor":"#333333","svg.fonttype":"none"})
BLUE,GREEN,ORANGE,GREY,LGREY,RED="#2c6fbb","#2a9d6f","#e08214","#9a9a9a","#cfcfcf","#c0392b"

fig=plt.figure(figsize=(15.5,15.2))
gs=fig.add_gridspec(3,2,height_ratios=[1.0,0.82,1.0],hspace=0.62,wspace=0.26,
                    left=0.075,right=0.965,top=0.945,bottom=0.055)
def tag(ax,t,x=-0.14,y=1.15): ax.text(x,y,t,transform=ax.transAxes,fontsize=21,fontweight="bold",va="top",ha="left")

# (a) carriage by category
axa=fig.add_subplot(gs[0,0])
cats=["Susceptible","Other R","MDR","Pre-XDR","XDR"]; ns=[19652,5624,10049,3050,1901]
allv=[1.55,2.17,2.50,3.44,3.00]; exclv=[1.08,1.72,2.28,3.05,3.00]
x=np.arange(5); w=0.40
axa.bar(x-w/2,allv,w,color=BLUE,label="All coding variants",zorder=3)
axa.bar(x+w/2,exclv,w,color=GREEN,label="Excluding D144N",zorder=3)
for xi,v in zip(x-w/2,allv): axa.text(xi,v+0.06,f"{v:.2f}",ha="center",va="bottom",fontsize=9.6,color=BLUE)
for xi,v in zip(x+w/2,exclv): axa.text(xi,v+0.06,f"{v:.2f}",ha="center",va="bottom",fontsize=9.6,color=GREEN)
axa.set_ylim(0,4.2); axa.set_ylabel("Isolates carrying a\nprotein-changing variant (%)",fontsize=12)
axa.set_xticks(x); axa.set_xticklabels([f"{c}\n(n={n:,})" for c,n in zip(cats,ns)],fontsize=10)
axa.set_title("Carriage by drug-resistance category",fontsize=13.5,pad=26)
axa.spines[["top","right"]].set_visible(False); axa.legend(frameon=False,fontsize=10,loc="upper left")
axa.text(0.985,0.985,"Cochran–Armitage trend\n(all-variant series)\nz = 8.14,  p = 3.9 × 10$^{-16}$",
         transform=axa.transAxes,ha="right",va="top",fontsize=9.3,color="#444",linespacing=1.3)
tag(axa,"a")

# (b) forest
axb=fig.add_subplot(gs[0,1])
rows=[("Carriers excl. D144N (pooled)",2.16,1.84,2.55,GREEN,"s"),
      ("Collection 3: TB Portals",2.05,1.39,3.02,BLUE,"o"),
      ("Collection 2: CRyPTIC",1.23,0.96,1.58,BLUE,"o"),
      ("Collection 1: Genotypic WGS",2.31,1.89,2.83,BLUE,"o"),
      ("All carriers (pooled, crude)",1.69,1.47,1.95,BLUE,"s"),
      ("D144N only",0.62,0.45,0.86,ORANGE,"s")]
yy=np.arange(len(rows))
for y,(lab,est,lo,hi,col,mk) in zip(yy,rows):
    axb.plot([lo,hi],[y,y],color=col,lw=2.6,zorder=2)
    axb.scatter([est],[y],s=135,marker=mk,color=col,zorder=3,edgecolor="white",linewidth=0.8)
    axb.text(hi+0.07,y,f"{est:.2f} ({lo:.2f}–{hi:.2f})",va="center",ha="left",fontsize=9.6,color=col)
axb.axvline(1.0,ls="--",color=GREY,lw=1.4,zorder=1)
axb.set_yticks(yy); axb.set_yticklabels([r[0] for r in rows],fontsize=10.2)
axb.set_xlim(0.4,3.9); axb.set_xlabel("Odds ratio for resistance (95% CI)",fontsize=12)
axb.set_title("Association with resistance; D144N opposite",fontsize=13.5,pad=26)
axb.spines[["top","right","left"]].set_visible(False); axb.tick_params(axis="y",length=0)
tag(axb,"b",x=-0.42)

# (c) lollipop
axc=fig.add_subplot(gs[1,:]); LEN=225
doms=[(1,64,LGREY,"N-terminal domain (1–64)","#555"),
      (65,137,"#bfe3d0","DNA-binding helix-turn-helix (65–137)\n44.7% of variant observations",GREEN),
      (138,225,"#c9ddf2","C-terminal / effector domain (138–225)",BLUE)]
dom_y,dom_h=0.0,0.9
for s,e,c,lab,tc in doms:
    axc.add_patch(plt.Rectangle((s,dom_y),e-s,dom_h,facecolor=c,edgecolor="#777",linewidth=1.0,zorder=2))
    axc.text((s+e)/2,dom_y+dom_h/2,lab,ha="center",va="center",fontsize=10,color=tc,zorder=3,linespacing=1.3)
subs=[("P20L",20,11,"other",0.0),("P59L",59,13,"other",0.0),("C71R",71,13,"hth",0.0),
      ("A85V",85,37,"hth",0.0),("V90A",90,15,"hth",0.0),("A98V",98,14,"hth",0.0),
      ("S100C",100,55,"hth",0.0),("K118E",118,11,"hth",0.0),("D144N",144,151,"d144n",0.0),
      ("I148V",148,26,"other",0.0),("C160R",160,13,"other",0.0),("V167D",167,12,"other",0.55),
      ("G173S",173,14,"other",0.0),("M185V",185,17,"other",0.0)]
colmap={"hth":GREEN,"other":BLUE,"d144n":ORANGE}; stem_base=dom_y+dom_h; maxc=151; stem_max=3.1
for name,pos,cnt,grp,ldy in subs:
    hgt=stem_base+(cnt/maxc)*stem_max; col=colmap[grp]
    axc.plot([pos,pos],[stem_base,hgt],color=col,lw=1.8,zorder=2)
    axc.scatter([pos],[hgt],s=70+(cnt/maxc)*240,color=col,zorder=3,edgecolor="white",linewidth=0.8)
    lab=name if grp!="d144n" else f"{name}  (n=151)"
    if ldy>0: axc.plot([pos,pos],[hgt+0.14,hgt+0.14+ldy],color=col,lw=0.7,ls=":",zorder=2)
    axc.text(pos,hgt+0.16+ldy,lab,ha="center",va="bottom",fontsize=9.3,color=col)
axc.set_xlim(-6,LEN+32); axc.set_ylim(-0.25,stem_base+stem_max+1.0)
axc.set_xlabel("ResR (Rv1830) residue position  —  stem height and marker size scale with isolate count",fontsize=11.5)
axc.set_yticks([])
axc.set_title("Recurrent substitutions cluster in the DNA-binding region; D144N lies just outside it",fontsize=13.5,pad=10)
axc.spines[["top","right","left"]].set_visible(False)
leg=[Line2D([0],[0],marker="o",color="w",markerfacecolor=GREEN,markersize=11,label="residues 65–137 (HTH)"),
     Line2D([0],[0],marker="o",color="w",markerfacecolor=BLUE,markersize=11,label="elsewhere"),
     Line2D([0],[0],marker="o",color="w",markerfacecolor=ORANGE,markersize=11,label="D144N (most common variant)")]
axc.legend(handles=leg,frameon=False,fontsize=10,loc="upper left",bbox_to_anchor=(0.005,0.99))
tag(axc,"c",x=-0.065,y=1.13)

# (d) MIC (rifampicin) — all-quality, 40,276 set
axd=fig.add_subplot(gs[2,0])
mics=["0.03","0.06","0.12","0.25","0.5","1","2","4","8"]; xm=np.arange(9)
non=[8.5,27.8,18.0,6.5,1.4,1.0,1.1,16.5,19.3]; car=[7.6,22.9,21.1,7.6,2.9,1.1,0.7,13.8,22.2]
axd.plot(xm,non,"-o",color=GREY,lw=2.4,ms=7,label="Non-carriers (n=14,622)")
axd.plot(xm,car,"-o",color=BLUE,lw=2.4,ms=7,label="Carriers (n=273)")
axd.set_xticks(xm); axd.set_xticklabels(mics,fontsize=10.5)
axd.set_xlabel("Rifampicin MIC (mg/L)",fontsize=12); axd.set_ylabel("% within group",fontsize=12)
axd.set_ylim(0,31); axd.set_title("No significant difference in rifampicin MIC",fontsize=13.5,pad=26)
axd.spines[["top","right"]].set_visible(False); axd.legend(frameon=False,fontsize=10.5,loc="upper center",bbox_to_anchor=(0.60,1.03))
axd.text(0.405,0.66,"Medians 0.12 mg/L\nMann–Whitney p = 0.21\n(not significant)",
         transform=axd.transAxes,ha="left",va="top",fontsize=9.5,color="#444",linespacing=1.3)
tag(axd,"d")

# (e) conceptual
axe=fig.add_subplot(gs[2,1]); axe.set_xlim(0,1); axe.set_ylim(0,1); axe.axis("off")
axe.set_title("Resistance vs resilience (mechanism from refs 3, 4)",fontsize=13.5,pad=26)
def box(y0,hh,face,edge,ttl,tcol,bl):
    axe.add_patch(FancyBboxPatch((0.02,y0),0.96,hh,boxstyle="round,pad=0.012,rounding_size=0.03",
                  facecolor=face,edgecolor=edge,linewidth=1.7,transform=axe.transAxes,zorder=2))
    axe.text(0.07,y0+hh-0.05,ttl,fontsize=10.8,fontweight="bold",color=tcol,va="top",ha="left",transform=axe.transAxes,zorder=3)
    for i,b in enumerate(bl):
        axe.text(0.085,y0+hh-0.125-i*0.082,b,fontsize=9.6,color="#333",va="top",ha="left",transform=axe.transAxes,zorder=3)
box(0.60,0.38,"#fbeceb",RED,"Classical drug resistance (e.g. rpoB, katG)",RED,
    ["•  Drug target blocked or prodrug not activated","•  Raises the MIC above the clinical breakpoint","•  Growth continues under sustained drug levels"])
box(0.19,0.38,"#e8f6ef",GREEN,"Antibiotic resilience (Rv1830 / ResR)",GREEN,
    ["•  MIC not raised in engineered models (ref. 3)","•  Shortens the post-antibiotic lag; faster regrowth","•  A survival phenotype, distinct from resistance"])
axe.text(0.03,0.125,"This study:",fontsize=9.8,fontweight="bold",color="#222",va="top",ha="left",transform=axe.transAxes)
axe.text(0.03,0.075,"Rv1830 coding variants are associated with higher drug-resistance\n"
         "categories; no MIC difference was significant after correction (8 drugs).",
         fontsize=9.4,color="#222",va="top",ha="left",transform=axe.transAxes,linespacing=1.3)
tag(axe,"e",x=-0.05)

fig.savefig("/tmp/claude-0/-home-claude/02f31dba-a3e8-57ad-a9db-867a7f601012/scratchpad/figure1.png",dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig("/tmp/claude-0/-home-claude/02f31dba-a3e8-57ad-a9db-867a7f601012/scratchpad/figure1.pdf",bbox_inches="tight",facecolor="white")
print("saved png+pdf")
