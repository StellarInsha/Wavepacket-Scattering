
# hbar = 1, m = 1
#
# Griffiths section 2.6
# V(x) = -V0 for |x| < a,   0 elsewhere

# THREE regions (unlike delta potential which had two) so:
#
#   Region I    x < -a      psi = e^(ikx) + r*e^(-ikx)     incident + reflected
#   Region II  -a < x < a   psi = C*e^(ilx) + D*e^(-ilx)   inside well
#   Region III  x > a       psi = t*e^(ikx)                 transmitted
#
# Two wave vectors:
#   k = sqrt(2E)       outside well (free particle)
#   l = sqrt(2(E+V0))  inside well  (particle speeds up: more KE because V is negative)
#
# Four boundary conditions (psi and psi' continuous at x=-a AND x=+a)
# give Griffiths eq 2.167:
#   t = e^(2ika) / [cos(2la) - i*(k^2+l^2)/(2kl)*sin(2la)]
#
# RAMSAUER-TOWNSEND: when sin(2la)=0, denominator=1 so T=1 (perfect transmission)
# This happens at special energies 

# Each scattering state NOT normalizable alone.
# Physical solution = wavepacket integral over k, computed by VEGAS.

import numpy as np
import vegas
import matplotlib.pyplot as plt

V0    = 8.0
a     = 1.5
k0    = 3.0
sigma = 0.8
x0    = -12.0

k_min = max(0.1, k0 - 5*sigma)
k_max = k0 + 5*sigma

def phi(k):
    return np.exp(-(k - k0)**2 / (4*sigma**2)) * np.exp(-1j*k*x0)

def t_coeff(E):
    kk    = np.sqrt(2*E)
    ll    = np.sqrt(2*(E + V0))
    denom = np.cos(2*ll*a) - 1j*(kk**2+ll**2)/(2*kk*ll)*np.sin(2*ll*a)
    return np.exp(2j*kk*a) / denom

def r_coeff(E):
    kk = np.sqrt(2*E); ll = np.sqrt(2*(E+V0))
    return t_coeff(E)*np.exp(-2j*kk*a)*1j*(ll**2-kk**2)/(2*kk*ll)*np.sin(2*ll*a)

def scatt_state(x, k):
    E  = 0.5*k**2
    kk = np.sqrt(2*E); ll = np.sqrt(2*(E+V0))
    tc = t_coeff(E);   rc = r_coeff(E)
    if x < -a:
        return np.exp(1j*kk*x) + rc*np.exp(-1j*kk*x)
    elif x <= a:
        C = 0.5*tc*np.exp(1j*kk*a)*(1+kk/ll)*np.exp(-1j*ll*a)
        D = 0.5*tc*np.exp(1j*kk*a)*(1-kk/ll)*np.exp(1j*ll*a)
        return C*np.exp(1j*ll*x) + D*np.exp(-1j*ll*x)
    else:
        return tc*np.exp(1j*kk*x)

def compute_probability(x_vals, t_time):
    prob = np.zeros(len(x_vals))
    for i, x in enumerate(x_vals):
        def integrand(k_arr):
            k = k_arr[0]
            if k <= 0.01: return [0.0, 0.0]
            E   = 0.5*k**2
            val = phi(k)*scatt_state(x,k)*np.exp(-1j*E*t_time)/(2*np.pi)
            return [val.real, val.imag]
        integ = vegas.Integrator([[k_min, k_max]])
        integ(integrand, nitn=5, neval=800)
        r2 = integ(integrand, nitn=8, neval=2000)
        prob[i] = r2[0].mean**2 + r2[1].mean**2
    return prob

E0 = 0.5*k0**2
T0 = abs(t_coeff(E0))**2
R0 = abs(r_coeff(E0))**2
print(f"Finite well  V0={V0}, a={a}, k0={k0}")
print(f"At central momentum: T={T0:.3f}, R={R0:.3f}, R+T={R0+T0:.3f}")
print()

x_vals = np.linspace(-20, 20, 60)
times  = [0.0, 2.0, 4.5]
results = {}
for t in times:
    print(f"  t={t}...", flush=True)
    results[t] = compute_probability(x_vals, t)

peak = x_vals[np.argmax(results[0.0])]
print(f"\nt=0 peak at x={peak:.1f}  (expected {x0})")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f"Finite Well Scattering  V₀={V0}, a={a}\n|ψ(x,t)|²",
             fontsize=12, fontweight='bold')
titles = ['t=0  packet approaching',
          't=2  packet entering well',
          't=4.5  split: transmitted + reflected']

for ax, t, title in zip(axes, times, titles):
    ax.fill_betweenx([0,0.35], -a, a, alpha=0.12, color='orange', label='Well')
    ax.plot(x_vals, results[t], 'o-', color='C0', ms=4, lw=1.5, label='|ψ|²')
    ax.axvline(-a, color='orange', lw=1.5, ls='--', alpha=0.7)
    ax.axvline(+a, color='orange', lw=1.5, ls='--', alpha=0.7)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('x'); ax.set_ylabel('|ψ|²')
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('03_finite_well_scattering.png', dpi=150, bbox_inches='tight')
plt.show()
