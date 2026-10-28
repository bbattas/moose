import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import glob
from types import SimpleNamespace


def analytic_spec(a0,t,L,kappa,gbe):
    iw = 4 * kappa / (3 * gbe)
    M_GB = L * iw * 3 / 4
    A = a0 - (2 * np.pi * M_GB * gbe * t)
    return A

def aniso_props(L0,kappa,gbe0,mu,f):
    out = SimpleNamespace()
    g = gbe0 * f / (np.sqrt(kappa * mu))
    g2 = g * g
    # pg = a1*g2^4 + a2*g2^3 + a3*g2^2 + a4*g2 + a5
    pg_coeffs = [-3.0944, -1.8169, 10.323, -8.1819, 2.0033]
    pg = np.polyval(pg_coeffs, g2)
    gamma = 1/pg
    # f0_int = 0.0788*pg^6 - 0.4955*pg^5 + ... - 0.5563*pg + 0.2907
    f0_coeffs = [0.0788, -0.4955, 1.2244, -1.5281, 1.0686, -0.5563, 0.2907]
    f0_int = np.polyval(f0_coeffs, pg)
    iw = np.sqrt(kappa/mu) * np.sqrt(1/f0_int)
    # Write outs
    out.f = f
    out.gbe = gbe0 * f
    out.g = g
    out.pg = pg
    out.gamma = gamma
    out.f0 = f0_int
    out.iw = iw
    out.L = L0 * f
    out.kappa = kappa
    out.mu = mu
    return out


def backcalculate_iso(out):
    # Not all will line up since gamma is constant, but we need L/MGB and sigma as critical ones
    # Basically inputs are (L, sigma, iw, f0)
    # Constants/Defaults:
    length_scale = 1e-9
    time_scale = 1e-9
    jtoev = 6.24150974e18
    MGB = out.L * out.iw * 3 / 4
    GBmob = MGB * jtoev * (length_scale ** 4) / time_scale
    GBenergy = out.gbe / (jtoev * length_scale * length_scale)
    kappa = 0.75 * out.gbe / out.iw
    mu = 0.75 * out.gbe / (out.f0 * out.iw)
    iso_props = SimpleNamespace(gbe=out.gbe, gamma=1.5, f0=out.f0, iw=out.iw, L=out.L, GBmob=GBmob, GBenergy=GBenergy, kappa=kappa, mu=mu,)
    # print(f'High Miso ({f_high}):')
    for name, value in vars(iso_props).items():
        value = value.item() if hasattr(value, "item") else value
        print(f"{name:>6} = {value:.8g}")


def backcalculate_iso_kappaconst(out):
    # Not all will line up since gamma is constant, but we need L/MGB and sigma as critical ones
    # Basically inputs are (L, sigma, kappa, f0), then get iw from kappa and continue on
    # Constants/Defaults:
    length_scale = 1e-9
    time_scale = 1e-9
    jtoev = 6.24150974e18
    # get iw from kappa
    iw = 0.75 * out.gbe / out.kappa
    MGB = out.L * iw * 3 / 4
    GBmob = MGB * jtoev * (length_scale ** 4) / time_scale
    GBenergy = out.gbe / (jtoev * length_scale * length_scale)
    mu = 0.75 * out.gbe / (out.f0 * iw)
    iso_props = SimpleNamespace(gbe=out.gbe, gamma=1.5, f0=out.f0, iw=iw, L=out.L, GBmob=GBmob, GBenergy=GBenergy, kappa=out.kappa, mu=mu,)
    # print(f'High Miso ({f_high}):')
    for name, value in vars(iso_props).items():
        value = value.item() if hasattr(value, "item") else value
        print(f"{name:>6} = {value:.8g}")



# GLOBAL PROPERTIES:
L0 = 2.7815e-6
kappa = 2.100e07
gbe = 1.125e7
mu = 1.305e7



f_low = 0.58133
f_high = 0.984
f_range = np.linspace(0.5,1,50)


out_full = aniso_props(L0,kappa,gbe,mu,f_range)
out_low = aniso_props(L0,kappa,gbe,mu,f_low)
out_high = aniso_props(L0,kappa,gbe,mu,f_high)

print('BOUNDS: g (0.098, 0.766), gamma [0.53, 40]')

print('Internal Properties Check for the Anisotropic Model:')
print(f'Low Miso ({f_low}):')
for name, value in vars(out_low).items():
    value = value.item() if hasattr(value, "item") else value
    print(f"{name:>6} = {value:.8g}")

print(" ")

print(f'High Miso ({f_high}):')
for name, value in vars(out_high).items():
    value = value.item() if hasattr(value, "item") else value
    print(f"{name:>6} = {value:.8g}")



fig, ax = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=False,figsize=(6,4))
ax[0].plot(out_full.f, out_full.gamma)
ax[0].plot(out_low.f, out_low.gamma,marker='*',label=f'Low Miso ({f_low})')
ax[0].plot(out_high.f, out_high.gamma,marker='o',label=f'High Miso ({f_high})')
ax[0].set_xlabel('GBE Normalized')
ax[0].set_ylabel('Gamma')
ax[1].plot(out_full.f, out_full.iw)
ax[1].plot(out_low.f, out_low.iw,marker='*',label=f'Low Miso ({f_low})')
ax[1].plot(out_high.f, out_high.iw,marker='o',label=f'High Miso ({f_high})')
ax[1].set_xlabel('GBE Normalized')
ax[1].set_ylabel('IW')
fig.suptitle('Aniso Property Space')
fig.tight_layout()
fig.savefig('P01_gbe_prop_range.png',dpi=300)
plt.close(fig)



print(' ')
print(' ')
print('GBEvolution CALC FROM INPUT L/GBE/IW/f0s')

print('')
print('Backcalculated Iso Min Miso (0.5)')
backcalculate_iso(aniso_props(L0,kappa,gbe,mu,0.5))

print('')
print(f'Backcalculated Iso Low Miso ({f_low})')
backcalculate_iso(out_low)

print('')
print('Backcalculated Iso Mid Miso (0.75)')
backcalculate_iso(aniso_props(L0,kappa,gbe,mu,0.75))

print('')
print(f'Backcalculated Iso High Miso ({f_high})')
backcalculate_iso(out_high)

print('')
print('Backcalculated Iso Max Miso (1.0)')
backcalculate_iso(aniso_props(L0,kappa,gbe,mu,1.0))


print(' ')
print(' ')
print(' ')
print('MATLAB L/g for iso matched scaling:')

for f in [0.5, f_low, 0.75, f_high, 1.0]:
    props = aniso_props(L0,kappa,gbe,mu,f)
    print(f'f = {f}, g = {props.g:.4f}, L = {props.L:.4e}')
