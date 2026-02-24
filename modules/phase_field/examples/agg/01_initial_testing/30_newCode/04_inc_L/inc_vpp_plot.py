import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from typing import Optional, Union, Iterable
from pathlib import Path
from scipy.stats import linregress


def find_csv(
    file_part: Union[str, Path],
    subdir: Optional[Union[str, Path]] = None,
    *,
    recursive: bool = False,
    strict_single: bool = True,
) -> Path:
    """
    Find a CSV by partial name (or exact .csv), in cwd or a given subdir.

    - If `file_part` ends with ".csv": search for any file ending with that (e.g. "*foo.csv")
    - Else: search for any file containing that and ending in .csv (e.g. "*foo*.csv")

    Returns:
        Path to the match.

    Raises:
        FileNotFoundError: no matches
        FileExistsError: multiple matches and strict_single=True
    """
    base = Path(subdir) if subdir else Path.cwd()
    fp = Path(str(file_part))  # normalize inputs like 12 -> "12"
    token = fp.name  # ignore any accidental path parts in file_part

    pattern = f"*{token}" if token.lower().endswith(".csv") else f"*{token}*.csv"
    it = base.rglob(pattern) if recursive else base.glob(pattern)

    matches = sorted((p for p in it if p.is_file()), key=lambda p: p.name.lower())

    if not matches:
        # raise FileNotFoundError(f"No CSV matches pattern {pattern!r} in {str(base)!r}")
        return None

    if strict_single and len(matches) != 1:
        raise FileExistsError(
            f"Expected 1 match for {pattern!r} in {str(base)!r}, found {len(matches)}:\n"
            + "\n".join(str(m) for m in matches[:25])
            + ("\n..." if len(matches) > 25 else "")
        )

    return matches[0]


def x_at_level(df, xcol="x", ycol="contour", level=0.5, which="all"):
    """
    Return x location(s) where y(x) crosses 'level' using linear interpolation.

    Parameters
    ----------
    df : pd.DataFrame
    xcol, ycol : str
        Column names for x and the contour/field value.
    level : float
        Target contour value.
    which : {"all","first","last"}
        Which crossing(s) to return.

    Returns
    -------
    float or list[float] or None
    """
    d = df[[xcol, ycol]].dropna().sort_values(xcol)
    x = d[xcol].to_numpy()
    y = d[ycol].to_numpy()

    s = y - level

    # indices i where segment [i, i+1] crosses the level
    idx = np.where(s[:-1] * s[1:] < 0)[0]

    # also include exact hits (y == level) as "crossings"
    exact = np.where(s == 0)[0]
    xs = []

    # exact hits first
    xs.extend(x[exact].tolist())

    # interpolated crossings
    for i in idx:
        x0, x1 = x[i], x[i+1]
        y0, y1 = y[i], y[i+1]
        xi = x0 + (level - y0) * (x1 - x0) / (y1 - y0)
        xs.append(float(xi))

    if not xs:
        return None

    xs = sorted(xs)
    if which == "first":
        return xs[0]
    if which == "last":
        return xs[-1]
    return xs


def save_var(df,xy,num,var_dict,varname):
    var_dict[num] = {"r": df[xy].to_numpy(), "var": df[varname].to_numpy()}



def read_by_num(subdir,xy,num, var_dict):
    if xy == 'x':
        name = f"_horizontal_{num:04d}.csv"
    elif xy == 'y':
        name = f"_vertical_{num:04d}.csv"
    else:
        raise ValueError("xy in read_by_num needs to be 'x' or 'y'.")
    path = find_csv(name, subdir=subdir, strict_single=True)
    if path is not None:
        df = pd.read_csv(path)
        cont_05 = x_at_level(df, xcol=xy, ycol="contour", level=0.5, which="first")
        save_var(df,xy,num,var_dict,'gr1')
        if cont_05 is not None:
            return num, cont_05
        else:
            return num, 60
    else:
        return num, None

def run_subdir(subdir,xy):
    step = []
    contour = []
    gr1 = {}
    for i in range(1,101):
        id, cont = read_by_num(subdir=subdir,xy=xy,num=i,var_dict=gr1)
        if cont is None:
            continue
        cont_out = cont-60
        step.append(id)
        contour.append(cont_out)

    return np.array(step), np.array(contour), gr1


def linear_fit(x_in,y_in):
    '''
    Args:
        x_in: s
        y_in: c
    Returns:
        res (.slope, .intercept, .r, .r2)
    '''
    x = x_in*0.4
    y = y_in * y_in
    # (optional) remove NaNs/Infs so regression doesn't blow up
    mask = np.isfinite(x) & np.isfinite(y)
    x_fit = x[mask]
    y_fit = y[mask]
    res = linregress(x_fit, y_fit)
    res.r2 = res.rvalue * res.rvalue
    return res


def darken_color(color, factor=0.6):
    """
    factor < 1 makes color darker
    factor = 1 leaves it unchanged
    """
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)



s1x, c1x, g1x = run_subdir('01_bicr_iso',xy='x')
s1y, c1y, g1y = run_subdir('01_bicr_iso',xy='y')

s2x, c2x, g2x = run_subdir('02_a10_weight_constL',xy='x')
s2y, c2y, g2y = run_subdir('02_a10_weight_constL',xy='y')
s3x, c3x, g3x = run_subdir('03_a10_avg_constL',xy='x')
s3y, c3y, g3y = run_subdir('03_a10_avg_constL',xy='y')

s4x, c4x, g4x = run_subdir('04_a10_weight_nodL_kernelConstL',xy='x')
s4y, c4y, g4y = run_subdir('04_a10_weight_nodL_kernelConstL',xy='y')
s5x, c5x, g5x = run_subdir('05_a10_avg_nodL_kernelConstL',xy='x')
s5y, c5y, g5y = run_subdir('05_a10_avg_nodL_kernelConstL',xy='y')

s6x, c6x, g6x = run_subdir('06_a10_weight_anisoL_withKernel',xy='x')
s6y, c6y, g6y = run_subdir('06_a10_weight_anisoL_withKernel',xy='y')
s7x, c7x, g7x = run_subdir('07_a10_avg_anisoL_withKernel',xy='x')
s7y, c7y, g7y = run_subdir('07_a10_avg_anisoL_withKernel',xy='y')

all_lbls = ['Iso', 'A10%, Const L, WA', 'A10%, Const L, A',
            'A10%, Scaled L, WA', 'A10%, Scaled L, A',
            'A10%, Scaled L w/derivs, WA','A10%, Scaled L w/derivs, A']
all_cols_ic = ['Black','C0','C0','C1','C1','C2','C2']
all_lines = ['-','-','--','-','--','-','--']
subabc = ['ISO','a', 'b', 'c', 'd', 'e', 'f']
all_cols = all_cols_ic.copy()
for i in [2, 4, 6]:
    all_cols[i] = darken_color(all_cols[i], factor=0.8)


all_sx = [s1x, s2x, s3x, s4x, s5x, s6x, s7x]
all_cx = [c1x, c2x, c3x, c4x, c5x, c6x, c7x]
fitsx = []
fitsx1 = []
fitsx2 = []
fitsx3 = []
fitsx4 = []
for s,c in zip(all_sx, all_cx):
    fitsx.append(linear_fit(s,c))
    fitsx1.append(linear_fit(s[:45],c[:45]))
    fitsx2.append(linear_fit(s[45:90],c[45:90]))
    fitsx3.append(linear_fit(s[:15],c[:15]))
    fitsx4.append(linear_fit(s[75:90],c[75:90]))

all_sy = [s1y, s2y, s3y, s4y, s5y, s6y, s7y]
all_cy = [c1y, c2y, c3y, c4y, c5y, c6y, c7y]
fitsy = []
fitsy1 = []
fitsy2 = []
fitsy3 = []
fitsy4 = []
for s,c in zip(all_sy, all_cy):
    fitsy.append(linear_fit(s,c))
    fitsy1.append(linear_fit(s[:45],c[:45]))
    fitsy2.append(linear_fit(s[45:90],c[45:90]))
    fitsy3.append(linear_fit(s[:15],c[:15]))
    fitsy4.append(linear_fit(s[75:90],c[75:90]))

outdir = Path("pics")
outdir.mkdir(parents=True, exist_ok=True)

###########################################################################################
# r
figpath = outdir / "P01x_inc_all_radius_contour_at05.png"
plt.figure()
for i,x in enumerate(all_sx):
    tx = x * 0.4
    plt.plot(tx, all_cx[i], color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.title('Horizontal (f=1.1)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

figpath = outdir / "P01y_inc_all_radius_contour_at05.png"
plt.figure()
for i,x in enumerate(all_sy):
    tx = x * 0.4
    plt.plot(tx, all_cy[i], color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.title('Vertical (f=0.9)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

###########################################################################################
# R^2
figpath = outdir / "P02x_inc_all_r2_contour_at05.png"
plt.figure()
for i,x in enumerate(all_sx):
    tx = x * 0.4
    plt.plot(tx, all_cx[i] * all_cx[i], color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('Horizontal (f=1.1)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

figpath = outdir / "P02y_inc_all_r2_contour_at05.png"
plt.figure()
for i,x in enumerate(all_sy):
    tx = x * 0.4
    plt.plot(tx, all_cy[i] * all_cy[i], color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('Vertical (f=0.9)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

###########################################################################################
# R^2 Fits
figpath = outdir / "P03x_inc_all_r2fits.png"
plt.figure()
for i,x in enumerate(all_sx):
    tx = x * 0.4
    plt.plot(tx, tx*fitsx[i].slope + fitsx[i].intercept, color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('Horizontal (f=1.1) Fits')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

figpath = outdir / "P03y_inc_all_r2fits.png"
plt.figure()
for i,x in enumerate(all_sy):
    tx = x * 0.4
    plt.plot(tx, tx*fitsy[i].slope + fitsy[i].intercept, color=all_cols[i], linestyle=all_lines[i], label=all_lbls[i])
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('Vertical (f=0.9) Fits')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


###########################################################################################
# INDIVIDUAL FIT CHECKS
figpath = outdir / f"P04a_inc_constL_r2fits.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 1
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 2
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Constant L')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



figpath = outdir / f"P04b_inc_varL_noDeriv_r2fits.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 3
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 4
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Variable L (without derivs/kernel)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')




figpath = outdir / f"P04c_inc_varL_fullKernel_r2fits.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 5
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 6
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + fitsx[i].intercept, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + fitsy[i].intercept, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Variable L (with derivs & kernel)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')







###########################################################################################
# INDIVIDUAL FIT CHECKS- lookign at the slopes especially- int is actual not fit int
# basically project from the actual starting value
iso_r2 = all_cx[i]*all_cx[i]
intc = iso_r2[0]
figpath = outdir / f"P05a_inc_constL_r2fits_projected.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 1
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 2
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Constant L- \n' \
'Fit slope projected from IC')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



figpath = outdir / f"P05b_inc_varL_noDeriv_r2fits_projected.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 3
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 4
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Variable L (without derivs/kernel)- \n' \
'Fit slope projected from IC')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')




figpath = outdir / f"P05c_inc_varL_fullKernel_r2fits_projected.png"
plt.figure()
plt.plot(all_sx[0]*0.4, all_cx[0]*all_cx[0], color='black',linestyle='-',label='ISO X')
plt.plot(all_sy[0]*0.4, all_cy[0]*all_cy[0], color='gray',linestyle='--',label='ISO Y')
# Weighted avg
i = 5
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C0',linestyle='-',label='WA X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C0',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C1',linestyle='-',label='WA Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C1',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
# Avg
i = 6
plt.plot(all_sx[i]*0.4, all_cx[i]*all_cx[i], color='C2',linestyle='--',label='A X')
plt.plot(all_sx[i]*0.4, all_sx[i]*0.4*fitsx[i].slope + intc, color='C2',
         linestyle=':', label=f'Fit X R={fitsx[i].r2:.3f}')
plt.plot(all_sy[i]*0.4, all_cy[i]*all_cy[i], color='C3',linestyle='--',label='A Y')
plt.plot(all_sy[i]*0.4, all_sy[i]*0.4*fitsy[i].slope + intc, color='C3',
         linestyle=':', label=f'Fit Y R={fitsy[i].r2:.3f}')
#
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.title('A=10%, Variable L (with derivs & kernel)- \n' \
'Fit slope projected from IC')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')




###########################################################################################
# CHECK THE SLOPES vs GBE
gbe_iso = 4.6075e+06
gbe_v = gbe_iso * 0.9
gbe_h = gbe_iso * 1.1
slope_iso = fitsx[0].slope
slope_v = slope_iso * gbe_v / gbe_iso
slope_h = slope_iso * gbe_h / gbe_iso

all_lbls = ['Iso', 'A10%, Const L, WA', 'Const L, A',
            'Scaled L, WA', 'Scaled L, A',
            'Scaled L w/derivs, WA','Scaled L w/derivs, A']
# all_slopes_lbl = ['Iso', 'A10%, Const L', 'Scaled L Low','Iso L, High', 'Scaled L, High']
# all_slopes_fit = [fitsx[10].slope, fits[2].slope, fits[8].slope, fits[3].slope, fits[9].slope]
# all_slopes_expected = [slope_iso, slope_low, slope_low, slope_high, slope_high]

sz = [64,64,36,64,36,36,36,36]


all_cols = ['black','C0','C1','C2','C3','C4','C5']
leg_handles = [
    Line2D([0], [0], marker='x', color='black', linestyle='None',
           markersize=8, label='X- High Energy'),
    Line2D([0], [0], marker='o', color='black', linestyle='None',
           markersize=8, label='Y- Low Energy'),
    Patch(facecolor=all_cols[0], edgecolor=all_cols[0], label=all_lbls[0]),
    Patch(facecolor=all_cols[1], edgecolor=all_cols[1], label=all_lbls[1]),
    Patch(facecolor=all_cols[2], edgecolor=all_cols[2], label=all_lbls[2]),
    Patch(facecolor=all_cols[3], edgecolor=all_cols[3], label=all_lbls[3]),
    Patch(facecolor=all_cols[4], edgecolor=all_cols[4], label=all_lbls[4]),
    Patch(facecolor=all_cols[5], edgecolor=all_cols[5], label=all_lbls[5]),
    # Patch(facecolor=all_cols[6], edgecolor=all_cols[6], label=all_lbls[6]),
]



# Expected grain growth rate vs actual
figpath = outdir / "P06_inc_dr2_actualExpected.png"
ran = np.linspace(-53,-43.5,100)
xf = s2x*0.4
plt.figure()
# for i,lb in enumerate(all_slopes_lbl):
#     plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.scatter(slope_iso,fitsx[0].slope,c='black',marker='x',s=64,label=all_lbls[0])
plt.scatter(slope_iso,fitsy[0].slope,c='black',marker='o',s=64,label='_')
#
for i in [1,2,3,4,5]:
    plt.scatter(slope_h,fitsx[i].slope,color=all_cols[i],marker='x',s=sz[i],label=all_lbls[i])
    plt.scatter(slope_v,fitsy[i].slope,color=all_cols[i],marker='o',s=sz[i],label='_')

plt.plot(ran,ran,'--',c='black',label='_')
plt.plot(ran,ran[::-1],'--',c='gray',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso)')
plt.legend(handles=leg_handles, loc='upper right')
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')







# Expected grain growth rate vs actual
figpath = outdir / "P07a_inc_dr2_actualExpected_half.png"
ran = np.linspace(-55.8,-41,100)
xf = s2x*0.4
plt.figure()
# for i,lb in enumerate(all_slopes_lbl):
#     plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.scatter(slope_iso,fitsx1[0].slope,c='black',marker='x',s=64,label=all_lbls[0])
plt.scatter(slope_iso,fitsy1[0].slope,c='black',marker='o',s=64,label='_')
#
for i in [1,2,3,4,5]:
    plt.scatter(slope_h,fitsx1[i].slope,color=all_cols[i],marker='x',s=sz[i],label=all_lbls[i])
    plt.scatter(slope_v,fitsy1[i].slope,color=all_cols[i],marker='o',s=sz[i],label='_')

plt.plot(ran,ran,'--',c='black',label='_')
plt.plot(ran,ran[::-1],'--',c='gray',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso), Fit to first half of time')
plt.legend(handles=leg_handles, loc='upper right')
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


# Expected grain growth rate vs actual
figpath = outdir / "P07b_inc_dr2_actualExpected_half.png"
ran = np.linspace(-53.5,-43.5,100)
xf = s2x*0.4
plt.figure()
# for i,lb in enumerate(all_slopes_lbl):
#     plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.scatter(slope_iso,fitsx2[0].slope,c='black',marker='x',s=64,label=all_lbls[0])
plt.scatter(slope_iso,fitsy2[0].slope,c='black',marker='o',s=64,label='_')
#
for i in [1,2,3,4,5]:
    plt.scatter(slope_h,fitsx2[i].slope,color=all_cols[i],marker='x',s=sz[i],label=all_lbls[i])
    plt.scatter(slope_v,fitsy2[i].slope,color=all_cols[i],marker='o',s=sz[i],label='_')

plt.plot(ran,ran,'--',c='black',label='_')
plt.plot(ran,ran[::-1],'--',c='gray',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso), Fit to second half of time')
plt.legend(handles=leg_handles, loc='upper left')
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



# Expected grain growth rate vs actual
figpath = outdir / "P08a_inc_dr2_actualExpected_15.png"
ran = np.linspace(-63,-37,100)
xf = s2x*0.4
plt.figure()
# for i,lb in enumerate(all_slopes_lbl):
#     plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.scatter(slope_iso,fitsx3[0].slope,c='black',marker='x',s=64,label=all_lbls[0])
plt.scatter(slope_iso,fitsy3[0].slope,c='black',marker='o',s=64,label='_')
#
for i in [1,2,3,4,5]:
    plt.scatter(slope_h,fitsx3[i].slope,color=all_cols[i],marker='x',s=sz[i],label=all_lbls[i])
    plt.scatter(slope_v,fitsy3[i].slope,color=all_cols[i],marker='o',s=sz[i],label='_')

plt.plot(ran,ran,'--',c='black',label='_')
plt.plot(ran,ran[::-1],'--',c='gray',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso), Fit to first 15 steps')
plt.legend(handles=leg_handles, loc='lower left')
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


# Expected grain growth rate vs actual
figpath = outdir / "P08b_inc_dr2_actualExpected_15.png"
ran = np.linspace(-56.5,-41,100)
xf = s2x*0.4
plt.figure()
# for i,lb in enumerate(all_slopes_lbl):
#     plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.scatter(slope_iso,fitsx4[0].slope,c='black',marker='x',s=64,label=all_lbls[0])
plt.scatter(slope_iso,fitsy4[0].slope,c='black',marker='o',s=64,label='_')
#
for i in [1,2,3,4,5]:
    plt.scatter(slope_h,fitsx4[i].slope,color=all_cols[i],marker='x',s=sz[i],label=all_lbls[i])
    plt.scatter(slope_v,fitsy4[i].slope,color=all_cols[i],marker='o',s=sz[i],label='_')

plt.plot(ran,ran,'--',c='black',label='_')
plt.plot(ran,ran[::-1],'--',c='gray',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso), Fit to last 15 steps')
plt.legend(handles=leg_handles, loc='upper left')
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')
