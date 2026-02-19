import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Optional, Union, Iterable
from pathlib import Path


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
        return num, cont_05

def run_subdir(subdir):
    step = []
    contour = []
    gr1 = {}
    for i in range(1,101):
        id, cont = read_by_num(subdir=subdir,num=i,var_dict=gr1)
        cont_out = cont-60
        step.append(id)
        contour.append(cont_out)

    return np.array(step), np.array(contour), gr1

s1, c1, g1 = run_subdir('01_miso_low')
s2, c2, g2 = run_subdir('02_miso_high')
s3, c3, g3 = run_subdir('03_manual_gamma_low')
s4, c4, g4 = run_subdir('04_manual_gamma_high')


outdir = Path("pics")
outdir.mkdir(parents=True, exist_ok=True)


figpath = outdir / "P01_miso_radius_contour_at05.png"
plt.figure()
plt.plot(s1*0.4,c1,label='Miso Low')
plt.plot(s2*0.4,c2,label='Miso High')
plt.plot(s3*0.4,c3,':',label='Man. Gamma Low')
plt.plot(s4*0.4,c4,':',label='Man. Gamma High')
plt.xlabel('Time (s)')
plt.ylabel('Grain Radius')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


figpath = outdir / "P02_miso_drdt.png"
plt.figure()
plt.plot(s1[1:]*0.4,np.diff(c1)/0.4,label='Miso Low')
plt.plot(s2[1:]*0.4,np.diff(c2)/0.4,label='Miso High')
plt.plot(s3[1:]*0.4,np.diff(c3)/0.4,':',label='Man. Gamma Low')
plt.plot(s4[1:]*0.4,np.diff(c4)/0.4,':',label='Man. Gamma High')
plt.xlabel('Time (s)')
plt.ylabel('dr/dt')
plt.legend()
plt.savefig(figpath,dpi=500,transparent=True)
plt.close('all')


for num in list(s1):
    data1 = g1[num]
    data2 = g2[num]
    data3 = g3[num]
    data4 = g4[num]
    figpath = outdir / f"P03_miso_gr1_frame_{num:03d}.png"
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
