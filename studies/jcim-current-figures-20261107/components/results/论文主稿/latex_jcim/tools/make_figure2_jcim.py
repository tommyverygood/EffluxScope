"""Structural-split figure from frozen inputs; no model, split or embedding refit.

The two split maps share one external-density field and one coordinate system.
The performance panel shows one shared All baseline using archived row1 intervals.
Original scaffold/cluster memberships are highlighted; 453-score lists are not used.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.font_manager import FontProperties



ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import palette
C, FS = palette.C, palette.FIGURE_STYLE
DATA = ROOT / 'data_for_figures/split_density'
OUT = ROOT / 'figures/Figure2_split_separation'

def readcsv(name):
    with (DATA/name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

emb = readcsv('fig2_emb.csv')
xy = np.array([[float(r['tx']), float(r['ty'])] for r in emb])
ref = np.load(DATA/'fig2_ref_xy.npy')
d5 = np.array([float(r['joint_model_d5']) for r in emb])
masks = [np.array([r[k]=='True' for r in emb]) for k in ('removed_scaffold','is_km_ood')]
perf = readcsv('fig2_perf.csv')
leak = json.loads((DATA/'fig2_leak.json').read_text())
hist = json.loads((DATA/'fig2_hist.json').read_text())
ci = json.loads((DATA/'fig2_ci.json').read_text())
assert len(xy)==2624 and len(ref)==443
assert [int(m.sum()) for m in masks] == [62,544]
assert np.isfinite(xy).all() and np.isfinite(ref).all()
provenance=json.loads((DATA/'provenance.json').read_text())
for name,digest in provenance['source_inputs'].items():
    assert sha(DATA/Path(name).name)==digest,name

models = ['astrazeneca_gnn_mtl_mdck_er_2025','astrazeneca_gnn_mtl_nih_mdck_er_2025',
          'internal_v0_4_primary_score','deeppk_pgp_substrate_2026_08_05',
          'admet_ai_v2_pgp_broccatelli_inhibition']
labels = ['MDCK-ER (aligned)','NIH-MDCK-ER (aligned)','Weak-label score',
          'Public substrate score','P-gp inhibition control']
colors = [C['green'], palette.tint('green', .45), FS['weak_label'],
          palette.tint('red', .45), C['teal']]
palette.assert_palette(colors)
frozen={r['comparator_id']:r for r in readcsv('abcb1_public_predictor_metrics_v1.csv') if r['population']=='all_identity_clean'}
for row in perf:
    if row['split']=='All (100%)':
        for key,column in [('bacc','balanced_accuracy'),('auc','roc_auc')]:
            assert abs(float(row[key])-float(frozen[row['cid']][column]))<1e-6

plt.rcParams.update({'font.family':'Arial','font.size':7.4,'axes.titlesize':8,
    'axes.labelsize':7.4,'xtick.labelsize':7.4,'ytick.labelsize':7.4,
    'legend.fontsize':6.7,'axes.linewidth':.7,'lines.linewidth':.85,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','text.color':'#000000',
    'axes.labelcolor':'#000000','xtick.color':'#000000','ytick.color':'#000000',
    'savefig.facecolor':'white'})
fig=plt.figure(figsize=(7.25,4.90),facecolor='white')
fig.text(.5,.974,'Structural separation and performance across data splits',
         ha='center',va='center',fontsize=12*7.25/7.2,fontweight='bold',color=FS['ink'])
BG='#FFFFFF'; RED=C['amber']; BLUE=C['blue']   # RED = excluded compounds (amber, palette warning colour)
# Frozen joint axes, common density normalisation, equal x/y units.
joint=np.vstack([xy,ref]); lo=joint.min(0);hi=joint.max(0);pad=(hi-lo)*.045;lo-=pad;hi+=pad
edges=[np.linspace(lo[j],hi[j],137) for j in range(2)]
centers=[(e[:-1]+e[1:])/2 for e in edges]
h,_,_=np.histogram2d(xy[:,0],xy[:,1],bins=edges)
sigma=1.7;q=np.arange(-6,7);kernel=np.exp(-.5*(q/sigma)**2);kernel/=kernel.sum()
h=np.apply_along_axis(lambda v:np.convolve(v,kernel,mode='same'),0,h)
h=np.apply_along_axis(lambda v:np.convolve(v,kernel,mode='same'),1,h)
z=h.T/h.max()
# A fixed grey surround and muted density bands follow the requested reference.
grey=np.interp(z,[0,.025,.065,.14,.22,1],[1,.12,.30,.80,.96,.96])
dark=1-grey; rgb=1-dark[...,None]*(1-np.array(to_rgb(FS['border'])))
cmap=ListedColormap([palette.tint('blue', a) for a in (.08,.16,.30,.45,.70)])
levels=[.22,.32,.47,.64,.82,1.001]

maps=[]
def cloud(rect,title,mask=None):
    ax=fig.add_axes(rect)
    ax.imshow(rgb,origin='lower',extent=[lo[0],hi[0],lo[1],hi[1]],interpolation='bicubic',rasterized=True)
    ax.contourf(centers[0],centers[1],z,levels=levels,cmap=cmap,antialiased=True,zorder=2)
    # Keep the reference population individually visible and distinct from density.
    ax.scatter(ref[:,0],ref[:,1],s=1.7,facecolors='white',edgecolors=C['grey'],
               linewidths=.25,zorder=3,rasterized=True)
    if mask is not None:
        ax.scatter(xy[mask,0],xy[mask,1],s=3.1,c=RED,edgecolors='white',
                   linewidths=.12,zorder=4,rasterized=True)
    ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_aspect('equal');ax.axis('off')
    maps.append(ax)
    return ax


# a/b: the shared background is retained only where it provides split context.
cloud((.035,.545,.27,.340),'Murcko scaffold exclusion\n62 compounds set aside',masks[0])
cloud((.355,.545,.27,.340),'KMeans structural split\n544 compounds set aside',masks[1])

def graph(rect,title):
    ax=fig.add_axes(rect,facecolor=BG)
    ax.grid(True,axis='y',color=FS['grid'],lw=.65);ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color(FS['border'])
    ax.tick_params(length=2,width=.5,pad=2)
    return ax

# c: unchanged D5 distribution, bins, percentile lines and split membership.
ax=graph((.715,.565,.255,.320),'Reference-relative distance')
bins=np.linspace(.25,.9,46)
ax.hist(d5,bins,color=FS['blue_fill'],edgecolor='none',label='All (2,624)')
ax.hist(d5[masks[0]],bins,histtype='step',color=RED,lw=.9,label='Scaffold-removed (62)')
ax.hist(d5[masks[1]],bins,histtype='step',color=BLUE,lw=.9,label='KMeans-OOD (544)')
for j,(name,x) in enumerate(hist['loo_q'].items()):
    ax.axvline(x,color=C['grey'],lw=.6,ls=':')
    ax.annotate(name,(x,.57-j*.15),xycoords=('data','axes fraction'),xytext=(0,0),
                textcoords='offset points',ha='center',fontsize=7.1,color=FS['ink'],
                bbox=dict(facecolor=BG,edgecolor='none',pad=.7))
ax.set_xlim(.25,.9);ax.set_xlabel(r'Structural distance $D_5$',labelpad=3)
ax.set_ylabel('No. of compounds',labelpad=3)
ax.set_xticks([.25,.5,.75]);ax.set_ylim(0,440)
ax.legend(loc='upper left',frameon=False,fontsize=7.1,handlelength=1.15,
          labelspacing=.4,borderaxespad=.3)

# d: one All baseline, then the four original split subsets.
# All uses the existing row1/scaffold interval, without averaging or resampling.
# The unused row2 All interval remains in the unmodified archived input.
groups=[('row1','All (100%)','All'),
        ('row1','Scaffold-clean','Scaffold\nclean'),
        ('row1','Scaffold-removed','Scaffold\nremoved'),
        ('row2','KMeans in-domain','KMeans\nin-domain'),
        ('row2','KMeans-OOD20','KMeans\nOOD20')]
all_ci_comparison={}
assert leak['row1|All (100%)']['n']==leak['row2|All (100%)']['n']==2624
assert leak['row1|All (100%)']['mean_sim']==leak['row2|All (100%)']['mean_sim']
for model in models:
    one=next(r for r in perf if r['row']=='row1' and r['split']=='All (100%)' and r['cid']==model)
    two=next(r for r in perf if r['row']=='row2' and r['split']=='All (100%)' and r['cid']==model)
    assert float(one['bacc'])==float(two['bacc'])
    c1,c2=ci[f'{model}|row1|All (100%)'],ci[f'{model}|row2|All (100%)']
    all_ci_comparison[model]={'bacc':float(one['bacc']),
        'displayed_row1':[c1['lo'],c1['hi']], 'retained_archive_row2':[c2['lo'],c2['hi']],
        'absolute_endpoint_differences':[abs(c1['lo']-c2['lo']),abs(c1['hi']-c2['hi'])]}

ax=graph((.090,.165,.880,.260),'Performance across structural splits')
all_intervals=[]; displayed_values=[]
for k,(model,col) in enumerate(zip(models,colors)):
    xs=np.arange(len(groups))+(k-2)*.12
    ys=[];low=[];high=[]
    for row,split,_ in groups:
        r=next(r for r in perf if r['row']==row and r['split']==split and r['cid']==model)
        value=float(r['bacc']); interval=ci[f'{model}|{row}|{split}']
        ys.append(value);low.append(value-interval['lo']);high.append(interval['hi']-value)
        all_intervals.append([interval['lo'],interval['hi']])
        displayed_values.append({'model':model,'source_row':row,'split':split,
            'n':leak[row+'|'+split]['n'],'mean_top5_similarity':leak[row+'|'+split]['mean_sim'],
            'bacc':value,'ci95_low':interval['lo'],'ci95_high':interval['hi']})
    ax.errorbar(xs,ys,yerr=[low,high],fmt='o',ms=3.1,mew=.4,mfc=col,mec='white',
                ecolor=col,color=col,elinewidth=.75,capsize=1.3,capthick=.6,zorder=3)
ax.axhline(.5,color=C['grey'],ls='--',lw=.7,zorder=1)
ax.set_ylim(.40,.87);ax.set_xlim(-.5,4.5);ax.set_yticks([.4,.5,.6,.7,.8])
ax.set_ylabel('Balanced accuracy',labelpad=3)
ax.set_xticks(range(len(groups)),[f'{label}\nn = {leak[row+"|"+split]["n"]:,}' for row,split,label in groups])
ax.tick_params(axis='x',length=0,labelsize=7.4)
for j,(row,split,_) in enumerate(groups):
    ax.text(j,-.350,f'{leak[row+"|"+split]["mean_sim"]:.3f}',transform=ax.get_xaxis_transform(),
            ha='center',va='top',fontsize=7.4,color=FS['ink'])
ax.text(-.50,-.350,'Similarity',transform=ax.get_xaxis_transform(),
        ha='right',va='top',fontsize=7.1,color=FS['ink'])
assert len(displayed_values)==25
assert len({(r['model'],r['split']) for r in displayed_values})==25
assert sum(r['split']!='All (100%)' for r in displayed_values)==20
assert min(v[0] for v in all_intervals)>.40 and max(v[1] for v in all_intervals)<.87

# Keep each panel letter and heading in one centered, baseline-aligned group.
heading_specs=[(.170,.916,'a','Murcko scaffold exclusion','62 compounds set aside'),
               (.490,.916,'b','KMeans structural split','544 compounds set aside'),
               (.8425,.916,'c','Reference-relative distance',None),
               (.530,.455,'d','Performance across structural splits',None)]
heading_artists=[]
heading_renderer=fig.canvas.get_renderer()
for cx,y,letter,title,subtitle in heading_specs:
    letter_font=FontProperties(family='Arial',size=10,weight='bold')
    title_font=FontProperties(family='Arial',size=8.4,weight='bold')
    letter_w=heading_renderer.get_text_width_height_descent(letter,letter_font,False)[0]
    title_w=heading_renderer.get_text_width_height_descent(title,title_font,False)[0]
    gap=5*fig.dpi/72
    left=cx-(letter_w+gap+title_w)/(2*fig.bbox.width)
    a=fig.text(left,y,letter,fontproperties=letter_font,ha='left',va='baseline',color=FS['ink'])
    t=fig.text(left+(letter_w+gap)/fig.bbox.width,y,title,fontproperties=title_font,
               ha='left',va='baseline',color=FS['ink'])
    heading_artists.append((letter,cx,y,a,t))
    if subtitle:
        fig.text(cx,.892,subtitle,fontsize=7.6,ha='center',va='baseline',color=FS['ink'])
map_handles=[Patch(fc=FS['blue_fill'],label='External density'),
             Line2D([],[],marker='o',ls='',ms=3,mfc='white',mec=C['grey'],mew=.6,label='Reference (443)'),
             Line2D([],[],marker='o',ls='',ms=3,mfc=RED,mec='white',mew=.2,label='Set aside')]
map_legend=fig.legend(handles=map_handles,loc='center',bbox_to_anchor=(.330,.515),ncol=3,
           frameon=False,fontsize=7.3,handlelength=1.1,handletextpad=.4,columnspacing=1.3)
handles=[Line2D([],[],ls='',marker='o',ms=3.5,mfc=c,mec='white',mew=.4,label=l) for c,l in zip(colors,labels)]
model_legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.50,.006),ncol=5,
           frameon=False,fontsize=7.1,handlelength=.7,handletextpad=.35,columnspacing=.9)

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
# No automatic tight crop: top title stays centered on the actual canvas.
text_outside=[]
for t in fig.findobj(match=Text):
    if not t.get_visible() or not t.get_text().strip(): continue
    b=t.get_window_extent(renderer)
    if b.x0<-.1 or b.y0<-.1 or b.x1>fig.bbox.width+.1 or b.y1>fig.bbox.height+.1:
        text_outside.append(t.get_text())
assert not text_outside,text_outside
legend_bounds={name:list(g.get_window_extent(renderer).transformed(fig.transFigure.inverted()).bounds)
               for name,g in [('maps',map_legend),('models',model_legend)]}
heading_geometry=[]
for letter,cx,y,a,t in heading_artists:
    ab,tb=a.get_window_extent(renderer),t.get_window_extent(renderer)
    center_error_pt=abs((ab.x0+tb.x1)/2-cx*fig.bbox.width)*72/fig.dpi
    assert center_error_pt<.1
    heading_geometry.append({'panel':letter,'heading':t.get_text(),'baseline_figure_y':y,
        'letter_and_title_baseline_equal':a.get_position()[1]==t.get_position()[1],
        'gap_pt':(tb.x0-ab.x1)*72/fig.dpi,'group_center_error_pt':center_error_pt})
for ext in ['pdf','svg','eps','png']:
    fig.savefig(OUT.with_suffix('.'+ext),dpi=450 if ext!='png' else 300,bbox_inches=None)

def array_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
receipt={'external_n':2624,'reference_n':443,'split_counts':[62,544],
         'layout':'a scaffold map; b KMeans map; c D5 histogram; d five-group shared-baseline performance',
         'figure_inches':list(fig.get_size_inches()),'density_grid':136,'gaussian_sigma_bins':sigma,
         'graph_background':'#FFFFFF','panel_heading_layout':'letter and title share one baseline; each group centered over its panel',
         'panel_heading_geometry':heading_geometry,
         'density_field':'one common external-density field reused across a/b',
         'density_normalisation':'same maximum of the full external set in both maps',
         'density_contour_levels':levels,'reference_symbols':'white circles with palette slate-blue outlines',
         'split_symbols':'amber points, original frozen memberships',
         'shared_all_ci_source':'row1/scaffold','row2_all_ci_displayed':False,
         'row2_all_ci_preserved_in_archive':True,'all_ci_recomputed_or_averaged':False,
         'all_ci_archived_comparison':all_ci_comparison,
         'all_ci_max_absolute_endpoint_difference':max(max(r['absolute_endpoint_differences']) for r in all_ci_comparison.values()),
         'old_displayed_point_count':30,'displayed_point_count':25,'nonall_displayed_point_count':20,
         'displayed_archived_intervals_unchanged':True,'displayed_values':displayed_values,
         'interval_extrema':[min(x[0] for x in all_intervals),max(x[1] for x in all_intervals)],
         'performance_ylim':[.40,.87], 'no_interval_truncation':True,
         'frozen_array_sha256':{name:array_sha(a) for name,a in [('xy',xy),('reference',ref),('D5',d5),
             ('scaffold_mask',masks[0]),('KMeans_mask',masks[1]),('density_z',z),('density_rgb',rgb),('bins',bins),
             ('hist_all',np.histogram(d5,bins)[0]),('hist_scaffold',np.histogram(d5[masks[0]],bins)[0]),
             ('hist_KMeans',np.histogram(d5[masks[1]],bins)[0])]},
         'source_input_sha256':{name:sha(DATA/Path(name).name) for name in provenance['source_inputs']},
         'text_outside_canvas':text_outside,'legend_bounds_figure_fraction':legend_bounds,
         'new_embedding_or_model_fit':False,
         'output_sha256':{OUT.with_suffix('.'+e).name:sha(OUT.with_suffix('.'+e)) for e in ['pdf','svg','eps','png']}}
(DATA/'render.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
