"""Figure 2: a snapshot of the method (a) and the results it produces (b).

Framework and typography follow the two-tier "snapshot -> derived result" layout of the reference figure
(typed nodes with a legend, labelled relations, dashed data-property boxes, a highlighted sub-graph with an arrow
down to soft gradient result containers), drawn in the manuscript palette:
dark purple = weak-label score (WL), purple tint = assay-aware score (AA), blue = retained / accepted, grey = rest.

Everything numeric is read from stored data and checked:
  * traced molecule (first external ID): features, weighted sum, sigmoid and calibration recomputed from the stored
    weights; the WL score must equal the stored score to 1e-6;
  * results panel: stored balanced accuracies (crossed_metrics.csv), coverage curves (metrics.csv), per-source real
    score-endpoint pairs and the regulatory 2 x 2 counts; n and Spearman rho must equal the paper.
Nothing is refitted.
"""
from pathlib import Path
import csv, json, hashlib, sys, io, collections
from statistics import median
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PROC = PROJECT / 'data/source_inputs'
DATA = PROJECT / 'data'
SRC = DATA / 'frozen_inputs'                    # archived frozen inputs, never fitted here
OUTPUT = PROJECT / 'outputs'
PROVENANCE = OUTPUT / 'overview_redraw_provenance.json'
PREVIOUS_PROVENANCE = None
sys.path.insert(0, str(HERE))
import palette                                           # manuscript palette, shared by all figures
from abcb1_screener_core_v0_1 import descriptor_array, morgan_array   # the frozen model's own featuriser

plt.rcParams.update({
    'font.family': 'Arial', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.cal': 'Arial:italic',
    'mathtext.sf': 'Arial', 'mathtext.tt': 'Arial', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'svg.fonttype': 'none', 'axes.linewidth': .6,
})

# ---- names (change here only) -------------------------------------------------------------
WL, AA = 'WL', 'AA'
NAME_WL, NAME_AA, POLICY = f'EffluxScope-{WL}', f'EffluxScope-{AA}', 'Full policy'

# ---- type scale: five sizes (canvas pt; the canvas is placed at 6.5 in, so x0.90 in print) ----
S_LABEL, S_TITLE, S_HEAD, S_BOLD, S_BODY = 12, 9.5, 8.6, 7.8, 7.2

# ---- colours: the reference figure's blue classes / pink entities / red relations, built from the palette ----------
T_ = palette.tint
BLUE, RED = palette.C['blue'], palette.C['red']
INK, SUB = '#1F2A44', '#000000'                           # preserve shape outlines; set supporting text in black
TEXT = '#000000'
EDGE = T_('blue', 0.45)                                   # card outlines
LANE = T_('blue', 0.55)                                   # outer panel borders
FOCUS = '#141B2D'                                         # highlighted sub-graph
HIER = T_('blue', 0.60)                                   # class hierarchy arrows (light blue)
INST = '#133166'                                          # instance arrows (dark blue, dotted)
REL = RED                                                 # object-property (relation) arrows and labels
C_WL = RED                                                # weak-label score (series colour in the results)
C_AA = T_('red', 0.45)                                    # assay-aware score
C_AA_EDGE = T_('red', 0.80)
C_ACC = BLUE                                              # full policy / selection family
C_OTHER = palette.C['grey']                               # historical baseline (not an EffluxScope score): neutral role colour (slate blue)
CLASS_FC, CLASS_EC = T_('blue', 0.30), T_('blue', 0.70)
ENT_FC, ENT_EC = T_('red', 0.26), T_('red', 0.85)         # entity fill / outline; data-property boxes use the outline
LOCK = palette.C['amber']
ENT = {k: (ENT_FC, ENT_EC) for k in ('mol', 'wl', 'aa', 'pol', 'eva', 'src')}
GRAD = [T_('red', 0.14), T_('blue', 0.18)]                # result containers: pink to blue

# ---- geometry (x in 0-100 units) --------------------------------------------------------------------
W = 100.0
W_E, STEP, X0 = 14.0, 15.9, 2.8                            # six entity slots
xk = lambda k: X0 + k * STEP
ck = lambda k: xk(k) + W_E / 2
Y_LEG = (8.0, 11.8)
Y_ROOT = (14.0, 17.6)
Y_R1 = (21.6, 25.6)
Y_R2 = (31.0, 35.0)
Y_ENT = (40.8, 45.8)
Y_PROP = (45.8, 58.4)
Y_BRK, LANE_A, LANE_B = 59.6, 61.6, 65.2
PA_BOTTOM = 68.8
PB = PA_BOTTOM + 1.8
CONT_H = 37.0
PB_H = 7.5 + CONT_H + 1.8
H = PB + PB_H + 1.2
TITLE_MARGIN = 6.0
TOTAL_H = H + TITLE_MARGIN
FIG_W_IN = 7.2

fig = plt.figure(figsize=(FIG_W_IN, FIG_W_IN * TOTAL_H / W))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(H, -TITLE_MARGIN); ax.axis('off')
RENDERER = fig.canvas.get_renderer()
FITS = []


def box(x, y, w, h, fc='white', ec=EDGE, lw=.7, r=1.1, z=2, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}',
                                fc=fc, ec=ec, lw=lw, zorder=z, ls=ls))


def T(x, y, s, size=S_BODY, weight='normal', color=TEXT, ha='center', va='center', z=6, fit=None):
    assert size in (S_LABEL, S_TITLE, S_HEAD, S_BOLD, S_BODY), size
    t = ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, zorder=z, linespacing=1.15)
    if fit:
        FITS.append((t, fit[0], fit[1]))
    return t


def text_width(t):
    bb = t.get_window_extent(RENDERER)
    (a, _), (b, _) = ax.transData.inverted().transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    return b - a


def dot(x, y, face, edge=None, hollow=False):
    if hollow:
        ax.scatter([x], [y], s=24, facecolors='white', edgecolors='#5E6C8C', linewidths=.7, zorder=7)
    else:
        ax.scatter([x], [y], s=24, color=face, edgecolors=edge or face, linewidths=.6, zorder=7)


def line(xs, ys, color=LANE, lw=.9, z=3, ls='-'):
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle='butt', zorder=z, ls=ls)


def arrow(x0, y0, x1, y1, color=LANE, lw=.9, head=(.38, .2), ls='-'):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), zorder=4,
                arrowprops=dict(arrowstyle=f'-|>,head_length={head[0]},head_width={head[1]}', color=color, lw=lw,
                                shrinkA=0, shrinkB=0, linestyle=ls))


def inset(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h + TITLE_MARGIN) / TOTAL_H, w / W, h / TOTAL_H])


def padlock(x, y, s=1.0, fc=LOCK):
    ax.add_patch(FancyBboxPatch((x + .45 * s, y), 1.2 * s, 1.5 * s, boxstyle='round,pad=0,rounding_size=.5',
                                fc='white', ec=INK, lw=.8, zorder=9))
    ax.add_patch(FancyBboxPatch((x, y + 1.0 * s), 2.1 * s, 1.6 * s, boxstyle='round,pad=0,rounding_size=.3',
                                fc=fc, ec=INK, lw=.5, zorder=10))


def ramp(t, hi='blue', amount=.50):
    """white -> palette tint ramp, t in [0, 1] (used for value-shaded cells)."""
    h = T_(hi, amount).lstrip('#'); r = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    t = float(min(max(t, 0), 1))
    return tuple(1 - t * (1 - c) for c in r)


sigmoid = lambda t: 1 / (1 + np.exp(-t))
logit = lambda p: np.log(p / (1 - p))
rd = lambda p: list(csv.DictReader(open(p, encoding='utf-8-sig', newline='')))

# ================================================================ data
rows = rd(SRC / 'biogen_frozen_scores_membership_2d.csv')
assert len(rows) == 2624
col = lambda k: np.array([float(r[k]) for r in rows])
s1, s3 = col('s1'), col('s3')
m1, m3 = col('retained_v01_margin_453').astype(bool), col('retained_v03_margin_453').astype(bool)
mf = col('retained_full_v04_453').astype(bool)
assert m1.sum() == m3.sum() == mf.sum() == 453
N_ALL, N_KEEP = len(rows), int(m1.sum())
memd = {r['external_row_id']: r for r in rows}

EX_ID = sorted(memd)[0]                                   # first external ID, chosen without scores or labels
PRED = PROC / 'traced_biogen_molecule.csv'
pred = {r['external_row_id']: r for r in rd(PRED)}
ex_smiles = pred[EX_ID]['standardized_isomeric_smiles']
ex_mol = Chem.MolFromSmiles(ex_smiles)
assert Chem.MolToInchiKey(ex_mol) == memd[EX_ID]['standardized_inchikey']
PAR = json.loads((SRC / 'features/frozen_model_parameters.json').read_text())
WEIGHTS = np.array([float(r['coefficient']) for r in rd(SRC / 'features/frozen_logistic_weights_all_2056.csv')])
assert WEIGHTS.shape == (2056,)
fp = np.asarray(morgan_array(ex_mol), float).ravel()
desc = np.asarray(descriptor_array(ex_mol), float).ravel()
zdesc = (desc - np.array(PAR['descriptor_scaler_mean'])) / np.array(PAR['descriptor_scaler_scale'])
sum_wl = PAR['logistic_intercept'] + float(np.concatenate([fp, zdesc]) @ WEIGHTS)
raw_wl = float(sigmoid(sum_wl)); u = min(max(raw_wl, 1e-6), 1 - 1e-6)
s1_ex = float(sigmoid(PAR['calibrator_alpha'] * logit(u) + PAR['calibrator_beta']))
assert abs(s1_ex - float(memd[EX_ID]['s1'])) < 1e-6
s3_ex = float(memd[EX_ID]['s3'])
flags_ex = {k: memd[EX_ID][k] == '1' for k in ('retained_v01_margin_453', 'retained_v03_margin_453', 'retained_full_v04_453')}


def mol_img(mol):
    mol = Chem.Mol(mol)
    rdDepictor.SetPreferCoordGen(True); rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    pts = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
    c = pts - pts.mean(0); ang = np.arctan2(*np.linalg.svd(c)[2][0][::-1])
    R = np.array([[np.cos(-ang), -np.sin(-ang)], [np.sin(-ang), np.cos(-ang)]])
    for i, (x, y) in enumerate(c @ R.T):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))
    d = rdMolDraw2D.MolDraw2DCairo(1800, 900)
    o = d.drawOptions(); o.useBWAtomPalette(); o.bondLineWidth = 15; o.fixedFontSize = 74; o.padding = .03
    d.DrawMolecule(mol); d.FinishDrawing()
    im = np.asarray(Image.open(io.BytesIO(d.GetDrawingText())).convert('RGB'))
    ink = np.where((im < 245).any(axis=2))
    y0, y1, x0, x1 = ink[0].min(), ink[0].max(), ink[1].min(), ink[1].max()
    return im[max(y0 - 10, 0):y1 + 11, max(x0 - 10, 0):x1 + 11]


asap = rd(PROC / 'abcb1_asap2025_sourceheldout_predictions_v1.csv')
grp = collections.defaultdict(list)
for r in asap:
    grp[r['standardized_inchikey']].append(r)
asap_y = np.array([median(float(r['log2_efflux_ratio']) for r in v) for v in grp.values()])
asap_wl = np.array([float(v[0]['v0_1_calibrated_probability']) for v in grp.values()])
asap_aa = np.array([float(v[0]['v0_3_candidate_probability']) for v in grp.values()])
assert len(grp) == 56 and round(spearmanr(asap_wl, asap_y)[0], 3) == 0.550 and round(spearmanr(asap_aa, asap_y)[0], 3) == 0.521
lang = [r for r in rd(PROC / 'abcb1_q1_direct_external_entity_predictions_v1.csv')
        if 'LANGTHALER2024' in r['source_id'] and 'WEGLER' not in r['source_id']
        and r['model_evaluation_eligible'].lower() == 'yes' and r['continuous_endpoint_median'] not in ('', 'nan')]
lang_x = np.array([float(r['frozen_probability']) for r in lang])
lang_y = np.log2(np.array([float(r['continuous_endpoint_median']) for r in lang]))
assert len(lang) == 11 and abs(spearmanr(lang_x, lang_y)[0] - 0.20956774) < 1e-6
roth = [r for r in rd(PROC / 'baseline_roth_external_predictions_v2_identity_audited.csv')
        if r['model_name'] == 'morgan_logistic' and r['transporter_gene'] == 'ABCB1'
        and r['identity_safe_external_metric_include'] == 'yes' and r['resistance_ratio'] not in ('', 'nan')]
roth_x = np.array([float(r['predicted_probability']) for r in roth])
roth_y = np.log2(np.array([float(r['resistance_ratio']) for r in roth]))
assert len(roth) == 12 and abs(spearmanr(roth_x, roth_y)[0] + 0.5709290716957522) < 1e-6
reg = json.loads((PROC / 'regulatory_counts.json').read_text(encoding='utf-8'))
REG = dict(tp=reg['tp'], fn=reg['fn'], fp=reg['fp'], tn=reg['tn'])
assert sum(REG.values()) == reg['n'] == 33
cm = {(r['mask'], r['classifier']): r for r in rd(DATA / 'crossed_metrics.csv')}
mrows = {}
for r in rd(DATA / 'metrics.csv'):
    mrows.setdefault((r['policy'], r['target']), r)
assert round(float(cm[('v03_margin', 'auxiliary_v03')]['balanced_accuracy']), 3) == .843
assert round(float(cm[('v01_margin', 'primary_v01')]['balanced_accuracy']), 3) == .637

# ================================================================ panel a: snapshot of the method
T(W / 2, -TITLE_MARGIN / 2, 'ABCB1 efflux scoring and evaluation', size=S_LABEL,
  weight='bold', color=TEXT, fit=(0, W))
box(1, 1, 98, PA_BOTTOM - 1, ec=LANE, lw=1.0, r=3.2, z=1)
T(W / 2, 5.3, 'a  A snapshot of the method', size=S_TITLE, weight='bold', color=TEXT)

# --- legend (two columns, as in the reference figure)
LEG = (2.4, 7.0, 34.6, 10.8)                                             # x, y, w, h
box(LEG[0], LEG[1], LEG[2], LEG[3], fc=T_('blue', 0.07), ec='none', r=1.6, z=2)
row_y = [LEG[1] + 2.0 + i * 2.3 for i in range(4)]
col1_x = 4.0


def legend_row(x, y, kind, label):
    if kind == 'class':
        box(x, y - 1.0, 2.2, 2.0, fc=CLASS_FC, ec=CLASS_EC, lw=.8, r=.3, z=5); xi = 3.2
    elif kind == 'entity':
        box(x, y - 1.0, 2.2, 2.0, fc=ENT_FC, ec=ENT_EC, lw=.8, r=.3, z=5); xi = 3.2
    elif kind == 'prop':
        box(x, y - 1.0, 2.2, 2.0, fc='white', ec=ENT_EC, lw=.9, r=.3, z=5, ls=(0, (2, 1.4))); xi = 3.2
    elif kind == 'lock':
        padlock(x, y - 1.5, .85); xi = 3.2
    else:
        col, sty = {'hier': (HIER, '-'), 'inst': (INST, (0, (1, 1.6))), 'rel': (REL, '-'), 'lab': (REL, (0, (3, 2)))}[kind]
        line([x, x + 3.6], [y, y], color=col, lw=1.1, ls=sty, z=5); arrow(x + 3.2, y, x + 3.9, y, color=col, head=(.3, .16))
        xi = 4.9
    t = T(x + xi, y, label, ha='left', z=6)
    return x + xi + text_width(t)


w1 = 0
for y, (kind, label) in zip(row_y, [('class', 'Class'), ('entity', 'Entity'), ('prop', 'Data property'), ('lock', 'Frozen')]):
    w1 = max(w1, legend_row(col1_x, y, kind, label))
col2_x = w1 + 2.6
w2 = 0
for y, (kind, label) in zip(row_y, [('hier', 'Class hierarchy'), ('inst', 'Instance'), ('rel', 'Object property'), ('lab', 'Labels')]):
    w2 = max(w2, legend_row(col2_x, y, kind, label))
assert w2 < LEG[0] + LEG[2] - .6, w2
lx = w2

# --- classes (grey), root to subclasses
def cls(cx, w, y, label, bold=False):
    box(cx - w / 2, y[0], w, y[1] - y[0], fc=CLASS_FC, ec=CLASS_EC, lw=.8, r=.5, z=4)
    T(cx, (y[0] + y[1]) / 2, label, size=S_BOLD if bold else S_BODY, weight='bold' if bold else 'normal', z=7,
      fit=(cx - w / 2, cx + w / 2))


def up(child_x, child_top, parent_x, parent_bottom, ymid, color=HIER):
    """child -> parent hierarchy arrow (arrowhead at the parent)."""
    line([child_x, child_x], [child_top, ymid], color=color, lw=.9, z=3)
    line([child_x, parent_x], [ymid, ymid], color=color, lw=.9, z=3)
    arrow(parent_x, ymid, parent_x, parent_bottom, color=color, lw=.9, head=(.4, .22))


ROOT_X = (ck(0) + ck(5)) / 2 + 6.0
R1 = [('Evidence', ck(0), 14.0), ('Frozen score', (ck(1) + ck(2)) / 2, 20.0),
      ('Selection rule', ck(3), 14.0), ('Evaluation', (ck(4) + ck(5)) / 2, 20.0)]
R2 = [('Biogen/Fang', 0, 0), ('Weak-label', 1, 1), ('Assay-aware', 2, 1), ('Composite policy', 3, 2),
      ('Matched-count test', 4, 3), ('Source test', 5, 3)]
cls(ROOT_X, 22.5, Y_ROOT, 'ABCB1 efflux evaluation', bold=True)
for name, cx, w in R1:
    cls(cx, w, Y_R1, name)
    up(cx, Y_R1[0], ROOT_X, Y_ROOT[1], 19.6)
for name, k, parent in R2:
    cls(ck(k), W_E, Y_R2, name)
    up(ck(k), Y_R2[0], R1[parent][1], Y_R1[1], 28.3)
    line([ck(k), ck(k)], [Y_ENT[0], Y_R2[1] + .3], color=INST, lw=1.1, ls=(0, (1, 1.6)), z=3)      # instance link
    arrow(ck(k), Y_R2[1] + 1.2, ck(k), Y_R2[1], color=INST, lw=1.1, head=(.4, .22))

# --- entities (tinted by meaning) with dashed data-property boxes
im_ex = mol_img(ex_mol); asp = im_ex.shape[1] / im_ex.shape[0]
sL = r'$s_1$ = ' + f'{s1_ex:.3f}'
sR = r'$s_3$ = ' + f'{s3_ex:.3f}'
ENTITIES = [
    ('mol', 0, 'Molecule 0001', None, False),
    ('wl', 1, NAME_WL, ['443 weak labels', 'logistic', 'calibrated', sL], True),
    ('aa', 2, NAME_AA, ['63 matched records', 'logistic', 'no calibration', sR], True),
    ('pol', 3, POLICY, ['cut-offs', 'P10–P90 band', 'chemical support', 'local rules'], True),
    ('eva', 4, 'BAcc at 453', ['list × classifier', 'ER ≥ 2 positive', 'cluster bootstrap'], False),
    ('src', 5, 'Four sources', ['WL, AA unchanged', 'ASAP · 56', 'Langthaler · 11', 'regulatory · 33', 'Roth · 12'], False),
]
for key, k, name, props, locked in ENTITIES:
    fc, ec = ENT[key]
    box(xk(k), Y_ENT[0], W_E, Y_ENT[1] - Y_ENT[0], fc=fc, ec=ec, lw=1.0, r=.5, z=4)
    T(ck(k), (Y_ENT[0] + Y_ENT[1]) / 2 + .1, name, size=S_BOLD, weight='bold', z=7, fit=(xk(k), xk(k) + W_E))
    box(xk(k), Y_PROP[0], W_E, Y_PROP[1] - Y_PROP[0], fc='white', ec=ec, lw=1.0, r=.4, z=3, ls=(0, (3, 2)))
    if locked:
        padlock(xk(k) - .6, Y_ENT[0] - 1.9, .8)
    if key == 'mol':
        wu = min(W_E - 3.0, 6.4 * asp); hu = wu / asp
        ia = inset(ck(k) - wu / 2, Y_PROP[0] + 1.0 + (6.4 - hu) / 2, wu, hu)
        ia.imshow(im_ex, aspect='equal', interpolation='lanczos'); ia.axis('off'); ia.set_zorder(9)
        for i, ln in enumerate(['2,056 features', 'ER measured']):
            T(ck(k), Y_PROP[0] + 8.7 + i * 2.2, ln, fit=(xk(k), xk(k) + W_E))
    else:
        for i, ln in enumerate(props):
            T(ck(k), Y_PROP[0] + 2.2 + i * 2.3, ln, color=TEXT,
              fit=(xk(k), xk(k) + W_E))

# --- relations (labelled), drawn in lanes under the property boxes
yb = Y_PROP[1]
line([ck(1), ck(1)], [yb, Y_BRK], color=REL, lw=1.0); line([ck(2), ck(2)], [yb, Y_BRK], color=REL, lw=1.0)
line([ck(1), ck(2)], [Y_BRK, Y_BRK], color=REL, lw=1.0)                                     # the two scores as a pair
bc = (ck(1) + ck(2)) / 2
for (x_from, x_to, lab) in [(ck(0) + 2.6, bc - 2.6, 'scoredBy'), (bc + 2.6, ck(3) - 2.6, 'usedBy'),
                            (ck(3) + 2.6, ck(4) - 2.6, 'keeps 453')]:
    first = lab == 'scoredBy'
    y_start = yb if x_from != bc + 2.6 else Y_BRK
    line([x_from, x_from], [y_start, LANE_A], color=REL, lw=1.0)
    line([x_from, x_to], [LANE_A, LANE_A], color=REL, lw=1.0)
    y_end = Y_BRK if first else yb
    arrow(x_to, LANE_A, x_to, y_end, color=REL, lw=1.0, head=(.4, .22))
    T((x_from + x_to) / 2, LANE_A + 1.7, lab, ha='center', color=REL, z=7)
line([ck(0) - 2.6, ck(0) - 2.6, ck(4) + 2.6], [yb, LANE_B, LANE_B], color=REL, lw=1.0, ls=(0, (3, 2)))   # labels path
arrow(ck(4) + 2.6, LANE_B, ck(4) + 2.6, yb, color=REL, lw=1.0, head=(.4, .22))
T(ck(2) + 6.0, LANE_B + 1.7, 'labels, evaluation only', color=REL, z=7)

# --- highlighted sub-graph (the evaluation) and the arrow down to the results
hx0, hx1 = xk(4) - 1.2, xk(5) + W_E + 1.2
ax.add_patch(FancyBboxPatch((hx0, Y_R2[0] - 1.4), hx1 - hx0, Y_PROP[1] + .9 - (Y_R2[0] - 1.4),
                            boxstyle='round,pad=0,rounding_size=.8', fc='none', ec=FOCUS, lw=1.6, zorder=8))

# ================================================================ panel b: results in gradient containers
box(1, PB, 98, PB_H, ec=LANE, lw=1.0, r=3.2, z=1)
T(W / 2, PB + 4.3, 'b  Results from the same pipeline', size=S_TITLE, weight='bold', color=TEXT)
CY0 = PB + 7.5


def container(x, w, caption):
    y, h = CY0, CONT_H
    patch = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=3', fc='none', ec='none', zorder=1.5)
    ax.add_patch(patch)
    im = ax.imshow(np.linspace(0, 1, 256)[None, :], extent=[x, x + w, y + h, y], cmap=LinearSegmentedColormap.from_list(
        'g', GRAD), aspect='auto', zorder=1.5, vmin=0, vmax=1)
    im.set_clip_path(patch)
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=3', fc='none', ec=T_('blue', 0.40), lw=.7, zorder=1.7))
    box(x + 1.3, y + 1.3, w - 2.6, h - 1.3 - 5.6, fc='white', ec='none', r=1.8, z=2)                     # inner white card
    T(x + w / 2, y + h - 2.8, caption, size=S_BODY, z=6)
    ax.set_xlim(0, W); ax.set_ylim(H, -TITLE_MARGIN)
    return x + 1.3, y + 1.3, w - 2.6, h - 6.9                                                              # card rect


CARD1 = container(2.5, 32.5, 'Equal count: 453 of 2,624')
CARD2 = container(36.8, 30.0, 'Across coverage')
CARD3 = container(68.6, 28.9, 'Other sources')

# --- b1: the evaluation grid of the method, filled with the stored balanced accuracies
cx0, cy0_, cw_card, ch_card = CARD1
LISTS = [('v01_margin', WL, C_WL, None), ('v03_margin', AA, C_AA, C_AA_EDGE), ('full_v04_pool_primary_margin', 'Full', C_ACC, None)]
CLS = [('primary_v01', WL), ('auxiliary_v03', AA)]
cw1, ch1, g1 = 8.4, 7.2, .9
cells_x = [cx0 + 12.0, cx0 + 12.0 + cw1 + .8]
for (cid, n), x in zip(CLS, cells_x):
    T(x + cw1 / 2, cy0_ + 2.2, 'by ' + n, color=SUB, fit=(x - .5, x + cw1 + .5))
vals_b1 = {}
best = max(float(r['balanced_accuracy']) for r in cm.values())
for i, (mid, n, colr, edge) in enumerate(LISTS):
    y = cy0_ + 4.3 + i * (ch1 + g1)
    dot(cx0 + 2.2, y + ch1 / 2, colr, edge); T(cx0 + 3.5, y + ch1 / 2, n + ' list', ha='left', fit=(cx0, cells_x[0]))
    for (cid, _), x in zip(CLS, cells_x):
        r = cm[(mid, cid)]; v = float(r['balanced_accuracy']); lo, hi = float(r['balanced_accuracy_ci95_low']), float(r['balanced_accuracy_ci95_high'])
        vals_b1[f'{mid}|{cid}'] = [round(v, 3), round(lo, 3), round(hi, 3)]
        tval = (min(max(v, .60), .85) - .60) / .25
        top = abs(v - best) < 1e-9
        ax.add_patch(FancyBboxPatch((x, y), cw1, ch1, boxstyle='round,pad=0,rounding_size=.5', fc=ramp(tval),
                                    ec=FOCUS if top else T_('blue', .6), lw=1.3 if top else .7, zorder=6))
        T(x + cw1 / 2, y + 2.7, f'{v:.3f}', size=S_BOLD, weight='bold', z=8, fit=(x, x + cw1))
        T(x + cw1 / 2, y + 5.3, f'{lo:.2f}–{hi:.2f}', color=TEXT, z=8, fit=(x, x + cw1))

# --- b2: balanced accuracy against coverage retained
cx0, cy0_, cw_card, ch_card = CARD2
TARGETS = ['coverage_05pct', 'coverage_10pct', 'fixed_full_coverage', 'coverage_20pct', 'coverage_25pct']
CURVES = [('v01_margin', WL, C_WL, C_WL), ('v03_margin', AA, C_AA, C_AA_EDGE), ('full_v04_pool_primary_margin', 'Full', C_ACC, C_ACC)]
AX_X, AX_Y, AX_W, AX_H = cx0 + 5.0, cy0_ + 3.4, cw_card - 12.0, ch_card - 11.0
c_ax = inset(AX_X, AX_Y, AX_W, AX_H)
ylo, yhi = .5, .95
curve_vals = {}
for pol, n, face, line_c in CURVES:
    pts = [mrows[(pol, t)] for t in TARGETS if (pol, t) in mrows]
    xs = np.array([float(p['coverage']) * 100 for p in pts]); ys = np.array([float(p['balanced_accuracy']) for p in pts])
    lo = np.array([float(p['balanced_accuracy_ci95_low']) for p in pts]); hi = np.array([float(p['balanced_accuracy_ci95_high']) for p in pts])
    curve_vals[n] = {'coverage_pct': [round(v, 2) for v in xs], 'bacc': [round(v, 3) for v in ys]}
    c_ax.fill_between(xs, lo, hi, color=line_c, alpha=.14, lw=0)
    c_ax.plot(xs, ys, color=line_c, lw=1.1, marker='o', ms=2.6, mfc=face, mec=line_c, mew=.5, zorder=4)
c_ax.axvline(100 * N_KEEP / N_ALL, color='#5E6C8C', lw=.6, ls=(0, (3, 2)), zorder=1)
c_ax.set_xlim(3, 27); c_ax.set_ylim(ylo, yhi)
c_ax.set_xticks([5, 10, 15, 20, 25]); c_ax.set_yticks([.5, .6, .7, .8, .9]); c_ax.set_yticklabels(['.5', '.6', '.7', '.8', '.9'])
c_ax.tick_params(length=2, width=.5, labelsize=S_BODY, pad=1.5, labelcolor=TEXT, color='#5E6C8C')
for sp in ['top', 'right']: c_ax.spines[sp].set_visible(False)
for sp in ['left', 'bottom']: c_ax.spines[sp].set_color('#5E6C8C')
c_ax.patch.set_alpha(0); c_ax.set_zorder(9)
xu = lambda c: AX_X + AX_W * (c - 3) / 24
yu = lambda v: AX_Y + AX_H * (1 - (v - ylo) / (yhi - ylo))
T(xu(100 * N_KEEP / N_ALL), AX_Y - 1.5, str(N_KEEP), color=SUB, fit=(cx0, cx0 + cw_card))
T(AX_X + AX_W / 2, AX_Y + AX_H + 4.8, 'coverage retained (%)', color=SUB, fit=(cx0, cx0 + cw_card))
ax.text(cx0 + 1.0, AX_Y + AX_H / 2, 'BAcc', rotation=90, fontsize=S_BODY, color=SUB, ha='center', va='center', zorder=6)
for pol, n, face, line_c in CURVES:
    yv = float(mrows[(pol, TARGETS[-1] if n != 'Full' else 'fixed_full_coverage')]['balanced_accuracy'])
    ye = yu(yv) + (1.2 if n == 'Full' else 0)
    dot(AX_X + AX_W + 1.5, ye, face, line_c if n == AA else None); T(AX_X + AX_W + 2.8, ye, n, ha='left', fit=(cx0, cx0 + cw_card))

# --- b3: real score-endpoint pairs with each source's statistic, and the real regulatory 2 x 2 counts
cx0, cy0_, cw_card, ch_card = CARD3


def real_scatter(x, y, xs, ys, colr, edge=None, size=4.6):
    ia = inset(x, y, size, size)
    pad = .08 * (ys.max() - ys.min())
    ia.scatter(xs, ys, s=4.6, color=colr, edgecolors=edge or 'none', linewidths=.4, zorder=3)
    ia.set_xlim(-.04, 1.04); ia.set_ylim(ys.min() - pad, ys.max() + pad); ia.set_xticks([]); ia.set_yticks([])
    for sp in ['top', 'right']: ia.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']: ia.spines[sp].set_color('#5E6C8C')
    ia.patch.set_alpha(0); ia.set_zorder(9)


def real_table(x, y, counts, cw_=2.3, ch_=1.7, gp=.3):
    vals = [[counts['tp'], counts['fn']], [counts['fp'], counts['tn']]]
    vmax = max(max(r) for r in vals)
    x0, y0 = x + 1.5, y + 1.6
    for i in range(2):
        for j in range(2):
            ax.add_patch(FancyBboxPatch((x0 + j * (cw_ + gp), y0 + i * (ch_ + gp)), cw_, ch_,
                                        boxstyle='round,pad=0,rounding_size=.3', fc=ramp(vals[i][j] / vmax),
                                        ec=T_('blue', .6), lw=.6, zorder=7))
            T(x0 + j * (cw_ + gp) + cw_ / 2, y0 + i * (ch_ + gp) + ch_ / 2, str(vals[i][j]), z=8)
    for j, ch in enumerate('+−'):
        T(x0 + j * (cw_ + gp) + cw_ / 2, y + .5, ch, color=SUB)
    for i, ch in enumerate('+−'):
        T(x + .5, y0 + i * (ch_ + gp) + ch_ / 2, ch, color=SUB)


rho = lambda a, b: spearmanr(a, b)[0]
fmt = lambda v: f'{v:.2f}'.replace('-', '−')
SRC_ROWS = [
    ('ASAP · WL', f'n = 56 · ρ = {fmt(rho(asap_wl, asap_y))}', 'MDCK-MDR1 ER', (asap_wl, asap_y, C_WL, None)),
    ('ASAP · AA', f'n = 56 · ρ = {fmt(rho(asap_aa, asap_y))}', 'MDCK-MDR1 ER', (asap_aa, asap_y, C_AA, C_AA_EDGE)),
    ('Langthaler · WL', f'n = 11 · ρ = {fmt(rho(lang_x, lang_y))}', 'MDCKII-MDR1 ER', (lang_x, lang_y, C_WL, None)),
    ('Roth · baseline', f'n = 12 · ρ = {fmt(rho(roth_x, roth_y))}', r'GI$_{50}$ resistance ratio', (roth_x, roth_y, C_OTHER, None)),
    ('Regulatory · WL', f'n = 33 · BAcc {reg["bacc"]:.2f}', 'strict substrate calls', None),
]
rh3 = ch_card / len(SRC_ROWS)
for i, (nm, st, ep, data) in enumerate(SRC_ROWS):
    y = cy0_ + i * rh3; yc = y + rh3 / 2
    if data:
        real_scatter(cx0 + .6, yc - 2.3, data[0], data[1], data[2], data[3])
    else:
        real_table(cx0 + .5, yc - 2.6, REG)
    T(cx0 + 7.6, y + 1.4, nm, size=S_BOLD, weight='bold', ha='left', fit=(cx0, cx0 + cw_card))
    T(cx0 + 7.6, y + 3.2, st, ha='left', fit=(cx0, cx0 + cw_card))
    T(cx0 + 7.6, y + 4.9, ep, color=SUB, ha='left', fit=(cx0, cx0 + cw_card))

# --- arrow from the highlighted evaluation to the results (cubic curve: leaves and arrives downwards)
from matplotlib.path import Path as MPath
from matplotlib.patches import FancyArrowPatch
x_s, y_s, x_e, y_e = hx1 - 4.0, Y_PROP[1] + .9, 71.0, CY0 - .4
ax.add_patch(FancyArrowPatch(path=MPath([(x_s, y_s), (x_s, y_s + 9.0), (x_e, y_e - 9.0), (x_e, y_e)],
                                        [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]),
                             arrowstyle='-|>', mutation_scale=9, color=FOCUS, lw=1.3, zorder=9, shrinkA=0, shrinkB=0))
ax.set_xlim(0, W); ax.set_ylim(H, -TITLE_MARGIN)

# ---- layout check: every registered text must lie inside its card (0.5 unit margin)
fig.canvas.draw()
rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
bad = []
for t, x0, x1 in FITS:
    bb = t.get_window_extent(rend)
    (a, _), (b, _) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    if a < x0 + .4 or b > x1 - .4:
        bad.append((t.get_text().replace('\n', ' / '), round(a, 1), round(b, 1), x0, x1))
print('overflow:', bad if bad else 'none', '| canvas height units', round(H, 1), '| legend ends at', round(lx, 1))
assert not bad, bad

OUTPUT.mkdir(exist_ok=True)
fig.savefig(OUTPUT / 'Figure_overview.pdf', dpi=400)
fig.savefig(OUTPUT / 'Figure_overview.png', dpi=240)
plt.close(fig)

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel = lambda p: str(Path(p).relative_to(PROJECT))
provenance = {
    'layout': 'snapshot of the method (typed nodes, legend, relations, highlighted sub-graph) above, results in gradient containers below',
    'names': {WL: 'primary weak-label model (formerly v0.1)', AA: 'auxiliary assay-aware model (formerly v0.3)',
              POLICY: 'composite abstention policy (formerly v0.4)'},
    'colours': {'WL': C_WL, 'AA': C_AA, 'AA_rim': C_AA_EDGE, 'retained/accepted': C_ACC, 'source': 'scripts/palette.py'},
    'traced_example': {
        'id': EX_ID, 'rule': 'first external row ID (sorted); chosen without scores or labels', 'smiles': ex_smiles,
        'wl': {'weighted_sum_plus_intercept': round(sum_wl, 4), 'sigmoid': round(raw_wl, 4), 'calibrated_recomputed': round(s1_ex, 6),
               'stored': float(memd[EX_ID]['s1']), 'equal_to_1e-6': True},
        'aa': {'stored': s3_ex}, 'retained': flags_ex,
        'smiles_source': rel(PRED), 'smiles_source_sha256': sha(PRED)},
    'results': {
        'equal_count_grid': {'source': rel(DATA / 'crossed_metrics.csv'), 'sha256': sha(DATA / 'crossed_metrics.csv'), 'values_bacc_lo_hi': vals_b1},
        'coverage_curves': {'source': rel(DATA / 'metrics.csv'), 'sha256': sha(DATA / 'metrics.csv'), 'values': curve_vals,
                            'note': 'duplicate rows removed (first kept); Full policy exists only up to its 453 pool'},
        'other_sources': {n: st for n, st, *_ in SRC_ROWS},
        'regulatory_counts': {'counts': REG, 'source': 'data/source_inputs/regulatory_counts.json'}},
    'no_refit': True,
    'presentation': {'title': 'ABCB1 efflux scoring and evaluation',
                     'title_colour': '#000000', 'title_size_pt': S_LABEL,
                     'title_and_panel_headings_centered': True,
                     'panel_letters_grouped_with_headings': True,
                     'panel_heading_bold': True,
                     'panel_heading_size_pt': S_TITLE,
                     'essential_supporting_labels_black': True,
                     'canvas_inches': list(fig.get_size_inches()),
                     'numeric_data_unchanged': True}}
if PREVIOUS_PROVENANCE is not None:
    for section in ('names', 'colours', 'traced_example', 'results', 'no_refit'):
        assert provenance[section] == PREVIOUS_PROVENANCE[section], section
PROVENANCE.parent.mkdir(exist_ok=True)
PROVENANCE.write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding='utf-8')
print('saved')
