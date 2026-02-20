import matplotlib.pyplot as plt
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


def save_var(df,num,var_dict,varname):
    var_dict[num] = {"x": df['x'].to_numpy(), "var": df[varname].to_numpy()}



def read_by_num(subdir,num, var_dict):
    name = f"_radial_{num:04d}.csv"
    path = find_csv(name, subdir=subdir, strict_single=True)
    if path is not None:
        df = pd.read_csv(path)
        cont_05 = x_at_level(df, xcol="x", ycol="contour", level=0.5, which="first")
        save_var(df,num,var_dict,'gr1')
        if cont_05 is not None:
            return num, cont_05
        else:
            return num, 60
    else:
        return num, None

def run_subdir(subdir):
    step = []
    contour = []
    gr1 = {}
    for i in range(1,101):
        id, cont = read_by_num(subdir=subdir,num=i,var_dict=gr1)
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



s1, c1, g1 = run_subdir('01_miso_low')
s2, c2, g2 = run_subdir('02_miso_high')
s3, c3, g3 = run_subdir('03_manual_gamma_low')
s4, c4, g4 = run_subdir('04_manual_gamma_high')
# L OoM
s5, c5, g5 = run_subdir('05_lowL_low')
s6, c6, g6 = run_subdir('06_lowL_high')
s7, c7, g7 = run_subdir('07_highL_low')
s8, c8, g8 = run_subdir('08_highL_high')
# Scaling L via gbe*Mgb = kappa*L with constant Mgb and kappa
s9, c9, g9 = run_subdir('09_manualL_low')
s10, c10, g10 = run_subdir('10_manualL_high')
# ISO
s11, c11, g11 = run_subdir('11_iso')

all_s = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11]
all_c = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11]
fits = []
for s,c in zip(all_s, all_c):
    fits.append(linear_fit(s,c))

outdir = Path("pics")
outdir.mkdir(parents=True, exist_ok=True)


figpath = outdir / "P01_miso_radius_contour_at05.png"
plt.figure()
plt.plot(s11*0.4,c11,c='black',label='Iso')
plt.plot(s1*0.4,c1,c='C0',label='Miso Low')
plt.plot(s2*0.4,c2,c='C1',label='Miso High')
plt.plot(s3*0.4,c3,':',c='C2',label='Man. Gamma Low')
plt.plot(s4*0.4,c4,':',c='C3',label='Man. Gamma High')
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


figpath = outdir / "P02_miso_drdt.png"
plt.figure()
plt.plot(s11[1:]*0.4,np.diff(c11)/0.4,c='black',label='Iso')
plt.plot(s1[1:]*0.4,np.diff(c1)/0.4,c='C0',label='Miso Low')
plt.plot(s2[1:]*0.4,np.diff(c2)/0.4,c='C1',label='Miso High')
plt.plot(s3[1:]*0.4,np.diff(c3)/0.4,':',c='C2',label='Man. Gamma Low')
plt.plot(s4[1:]*0.4,np.diff(c4)/0.4,':',c='C3',label='Man. Gamma High')
plt.xlabel('Time (s)')
plt.ylabel('dr/dt')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



figpath = outdir / "P03_lowhighL_radius.png"
plt.figure()
plt.plot(s11*0.4,c11,c='black',label='Iso')

plt.plot(s5*0.4,c5,'-',c='C0',label='Low L, Low')
plt.plot(s3*0.4,c3,'-',c='C1',label='L, Low')
plt.plot(s7*0.4,c7,'-',c='C2',label='High L, Low')

plt.plot(s6*0.4,c6,':',c='C3',label='Low L, High')
plt.plot(s4*0.4,c4,':',c='C4',label='L, High')
plt.plot(s8*0.4,c8,':',c='C5',label='High L, High')
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


figpath = outdir / "P04_lowhighL_drdt.png"
plt.figure()
plt.plot(s11[1:]*0.4,np.diff(c11)/0.4,c='black',label='Iso')

plt.plot(s5[1:]*0.4,np.diff(c5)/0.4,'-',c='C0',label='Low L, Low')
plt.plot(s3[1:]*0.4,np.diff(c3)/0.4,'-',c='C1',label='L, Low')
plt.plot(s7[1:]*0.4,np.diff(c7)/0.4,'-',c='C2',label='High L, Low')

plt.plot(s6[1:]*0.4,np.diff(c6)/0.4,'--',c='C3',label='Low L, High')
plt.plot(s4[1:]*0.4,np.diff(c4)/0.4,'--',c='C4',label='L, High')
plt.plot(s8[1:]*0.4,np.diff(c8)/0.4,'--',c='C5',label='High L, High')
plt.xlabel('Time (s)')
plt.ylabel('dr/dt')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


# Changing L or not
figpath = outdir / "P05_manualL_radius.png"
plt.figure()
plt.plot(s11*0.4,c11,c='black',label='Iso')

plt.plot(s3*0.4,c3,'-',c='C0',label='Iso L, Low')
plt.plot(s9*0.4,c9,'--',c='C1',label='Scaled L, Low')

plt.plot(s4*0.4,c4,'-',c='C2',label='Iso L, High')
plt.plot(s10*0.4,c10,'--',c='C3',label='Scaled L, High')
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


figpath = outdir / "P06_manualL_drdt.png"
plt.figure()
plt.plot(s11[1:]*0.4,np.diff(c11)/0.4,c='black',label='Iso')

plt.plot(s3[1:]*0.4,np.diff(c3)/0.4,'-',c='C0',label='Iso L, Low')
plt.plot(s9[1:]*0.4,np.diff(c9)/0.4,'-',c='C1',label='Scaled L, Low')

plt.plot(s4[1:]*0.4,np.diff(c4)/0.4,'--',c='C2',label='Iso L, High')
plt.plot(s10[1:]*0.4,np.diff(c10)/0.4,'--',c='C3',label='Scaled L, High')
plt.xlabel('Time (s)')
plt.ylabel('dr/dt')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



# R2
figpath = outdir / "P07_manualL_radius2.png"
plt.figure()
plt.plot(s11*0.4,c11*c11,c='black',label='Iso')

plt.plot(s3*0.4,c3*c3,'-',c='C0',label='Iso L, Low')
plt.plot(s9*0.4,c9*c9,'--',c='C1',label='Scaled L, Low')

plt.plot(s4*0.4,c4*c4,'-',c='C2',label='Iso L, High')
plt.plot(s10*0.4,c10*c10,'--',c='C3',label='Scaled L, High')
plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')



figpath = outdir / "P08_manualL_dr2dt.png"
plt.figure()
plt.plot(s11[1:]*0.4,np.diff(c11*c11)/0.4,c='black',label='Iso')

plt.plot(s3[1:]*0.4,np.diff(c3*c3)/0.4,'-',c='C0',label='Iso L, Low')
plt.plot(s9[1:]*0.4,np.diff(c9*c9)/0.4,'-',c='C1',label='Scaled L, Low')

plt.plot(s4[1:]*0.4,np.diff(c4*c4)/0.4,'--',c='C2',label='Iso L, High')
plt.plot(s10[1:]*0.4,np.diff(c10*c10)/0.4,'--',c='C3',label='Scaled L, High')
plt.xlabel('Time (s)')
plt.ylabel('d(r^2)/dt')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


# Changing L or not
figpath = outdir / "P09_manualL_dr2fit.png"
xf = s11*0.4
plt.figure()
plt.plot(xf,xf*fits[10].slope + fits[10].intercept,c='black',label=f'Iso R={fits[10].r2:.3f}')

plt.plot(xf,xf*fits[2].slope + fits[2].intercept,'-',c='C0',label=f'Iso L, Low R={fits[2].r2:.3f}')
plt.plot(xf,xf*fits[8].slope + fits[8].intercept,'--',c='C1',label=f'Scaled L, Low R={fits[8].r2:.3f}')

plt.plot(xf,xf*fits[3].slope + fits[3].intercept,'-',c='C2',label=f'Iso L, High R={fits[3].r2:.3f}')
plt.plot(xf,xf*fits[9].slope + fits[9].intercept,'--',c='C3',label=f'Scaled L, High R={fits[9].r2:.3f}')

plt.xlabel('Time (s)')
plt.ylabel('R^2')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')

# CHECK THE SLOPES vs GBE
gbe_iso = 4.6075e+06
gbe_low = 2.6293e+06
gbe_high = 4.5612e+06
slope_iso = fits[10].slope
slope_low = slope_iso * gbe_low / gbe_iso
slope_high = slope_iso * gbe_high / gbe_iso

all_slopes_lbl = ['Iso', 'Iso L, Low', 'Scaled L Low','Iso L, High', 'Scaled L, High']
all_slopes_fit = [fits[10].slope, fits[2].slope, fits[8].slope, fits[3].slope, fits[9].slope]
all_slopes_expected = [slope_iso, slope_low, slope_low, slope_high, slope_high]



# Changing L or not
ran = np.linspace(-49,-27,100)
figpath = outdir / "P10_manualL_dr2_actualExpected.png"
xf = s11*0.4
plt.figure()
for i,lb in enumerate(all_slopes_lbl):
    plt.scatter(all_slopes_expected[i],all_slopes_fit[i],label=all_slopes_lbl[i])
plt.plot(ran,ran,'--',c='black',label='_')

plt.xlabel('Expected Slope d(R^2)')
plt.ylabel('Measured Slope d(R^2)')
plt.title(r'Projected $\dfrac{d(R^2)}{dt} \propto \sigma$ (ref iso)')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


for num in list(s1):
    data1 = g1[num]
    data2 = g2[num]
    data3 = g3[num]
    data4 = g4[num]
    data9 = g9[num]
    data10 = g10[num]
    data11 = g11[num]
    figpath = outdir / f"V01_miso_gr1_frame_{num:03d}.png"
    plt.figure()
    plt.plot(data1['x'],data1['var'],label='Miso Low')
    plt.plot(data2['x'],data2['var'],label='Miso High')
    plt.plot(data3['x'],data3['var'],':',label='Man. Gamma Low')
    plt.plot(data4['x'],data4['var'],':',label='Man. Gamma High')
    plt.xlim([60,120])
    plt.ylim([-0.02,1.02])
    plt.xlabel('X')
    plt.ylabel('Gr1')
    plt.legend()
    plt.title(f'Timestep = {num:03d}')
    plt.savefig(figpath,dpi=300,transparent=True)
    plt.close('all')

    figpath2 = outdir / f"V02_misoL_gr1_frame_{num:03d}.png"
    plt.figure()
    plt.plot(data11['x'],data11['var'],c='black',label='Iso')
    plt.plot(data3['x'],data3['var'],label='Iso L, Low')
    plt.plot(data9['x'],data9['var'],label='Scaled L, Low')
    plt.plot(data4['x'],data4['var'],'--',label='Iso L, High')
    plt.plot(data10['x'],data10['var'],'--',label='Scaled L, High')
    plt.xlim([60,120])
    plt.ylim([-0.02,1.02])
    plt.xlabel('X')
    plt.ylabel('Gr1')
    plt.legend()
    plt.title(f'Timestep = {num:03d}')
    plt.savefig(figpath2,dpi=300,transparent=True)
    plt.close('all')
