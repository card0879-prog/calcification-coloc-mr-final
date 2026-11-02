import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- unified design system (all four figures inherit this) ----
mpl.rcParams.update({
    'font.family':'DejaVu Sans','font.size':9,
    'axes.titlesize':11,'axes.titleweight':'bold','axes.labelsize':9.5,
    'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,
    'axes.edgecolor':'#444444','axes.linewidth':0.8,
    'xtick.color':'#444444','ytick.color':'#444444',
    'axes.labelcolor':'#222222','text.color':'#222222',
    'figure.dpi':150,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42,
})

# verdict palette (terminology fixed across paper & figures)
V = {
    'PASS_strong'  : '#009E73',   # bluish green
    'PASS_moderate': '#66C2A5',   # light green
    'fail_coloc'   : '#E69F00',   # orange
    'noncoloc'     : '#D55E00',   # vermillion
    'uninstr'      : '#999999',   # grey
    'underpowered' : '#BBBBBB',   # light grey
}
# trait (series) colours — kept distinct from verdict colours
EXPO = '#0072B2'   # blue  : exposure / expression
OUT  = '#444444'   # grey  : outcome
GWL  = '#B22222'   # genome-wide line

LABEL = {
    'PASS_strong':'PASS (strong)','PASS_moderate':'PASS (moderate)',
    'fail_coloc':'fail-coloc','noncoloc':'non-colocalizing nominal',
    'uninstr':'uninstrumentable','underpowered':'underpowered',
}

def badge(ax, x, y, verdict, fs=8.5):
    ax.text(x, y, LABEL[verdict], transform=ax.transAxes, fontsize=fs,
            ha='left', va='top', color='white', fontweight='bold',
            zorder=30, bbox=dict(boxstyle='round,pad=0.3', fc=V[verdict], ec='none'))

def note(ax, x, y, txt, fs=8.2):
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=fs, ha='left', va='top',
            color='#222222', zorder=30, bbox=dict(boxstyle='round,pad=0.3', fc='white',
            ec='#CCCCCC', lw=0.6))
