import matplotlib.pyplot as plt
import numpy as np



G_MIN = 0.098546
G_MAX = 0.765691
GAMMA_MIN = 0.53
GAMMA_MAX = 40.0


def poly_g(g):
    g2 = g * g
    a1 = -3.0944
    a2 = -1.8169
    a3 = 10.323
    a4 = -8.1819
    a5 = 2.0033
    poly = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5
    return poly

def f0(g):
    pol = poly_g(g)
    f = (((((0.0788 * pol - 0.4955) * pol + 1.2244) * pol - 1.5281) * pol + 1.0686) * pol - 0.5563) * pol + 0.2907
    return f

def gamma_iw(gbe,kappa,m):
    g = gbe / (np.sqrt(kappa * m))
    pg = poly_g(g)
    gamma = 1 / pg
    f0_iw = f0(g)
    iw = (np.sqrt(kappa / m)) * (np.sqrt(1 / f0_iw))
    return gamma, iw

def gbe_from_mat(g,kappa,m):
    return g * (np.sqrt(kappa * m))

def gamma_iw_range(gbe,kappa,m,a):
    amags = np.array([1-a,1+a])
    g = gbe * amags / (np.sqrt(kappa * m))
    pg = poly_g(g)
    gamma = 1 / pg
    f0_iw = f0(g)
    iw = (np.sqrt(kappa / m)) * (np.sqrt(1 / f0_iw))
    print(f'Gamma ranges: {gamma[0]:.4f} to {gamma[1]:.4f}')
    print(f'IW ranges: {iw[0]:.4f} to {iw[1]:.4f}')
    # return gamma, iw





# Solve for simple params from MATLAB input
print('MATLAB Input based properties:')
kappa = 0.3
g = np.sqrt(2)/3
L0 = 0.833
m = 0.9375
a = 0.15

gbe = gbe_from_mat(g,kappa,m)
iso_gamma, iso_iw = gamma_iw(gbe,kappa,m)

print(f'g = {g}')
print(f'gbe = {gbe}')
print(f'kappa = {kappa}')
print(f'm = {m}')
print(f'L0 = {L0}')
print(f'amag = {a}')
print(' ')
print(f'ISO Gamma = {iso_gamma}')
print(f'ISO IW = {iso_iw}')
gamma_iw_range(gbe,kappa,m,a)

print(" ")
print(" ")


# # MY PROPERTIES
# print('My MOOSE property set for the combined 0.5-1 range (cant use this one for 1+cos)')
# kappa = 2.100e07
# gbe = 1.125e7
# L0 = 2.7815e-6
# m = 1.305e7
# a = 0.15
# g = gbe / (np.sqrt(kappa * m))

# iso_gamma, iso_iw = gamma_iw(gbe,kappa,m)

# print(f'g = {g}')
# print(f'gbe = {gbe}')
# print(f'kappa = {kappa}')
# print(f'm = {m}')
# print(f'L0 = {L0}')
# print(f'amag = {a}')
# print(' ')
# print(f'ISO Gamma = {iso_gamma}')
# print(f'ISO IW = {iso_iw}')
# gamma_iw_range(gbe,kappa,m,a)


# MY PROPERTIES
print('My MOOSE property set for the second test:')
print('(An older set for 1+cos instead of the 0.5-1 range)')
kappa = 2.07337e7
gbe = 4.60748e6
L0 = 1.15382e-6
m = 5.521269e6
a = 0.15
g = gbe / (np.sqrt(kappa * m))

iso_gamma, iso_iw = gamma_iw(gbe,kappa,m)

print(f'g = {g}')
print(f'gbe = {gbe}')
print(f'kappa = {kappa}')
print(f'm = {m}')
print(f'L0 = {L0}')
print(f'amag = {a}')
print(' ')
print(f'ISO Gamma = {iso_gamma}')
print(f'ISO IW = {iso_iw}')
gamma_iw_range(gbe,kappa,m,a)
