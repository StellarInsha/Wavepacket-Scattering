# Two types of solutions:
# 1. Bound state (E < 0): one state, already solved in Hamiltonian repo
# 2. Scattering states (E > 0): what this file computes
#
# Scattering state for particle coming from left with momentum k:
#   x < 0:  psi = e^(ikx) + r(k)*e^(-ikx)    incoming + reflected
#   x > 0:  psi = t(k)*e^(ikx)                transmitted
#
# Griffiths eq 2.141-2.142:
#   r(k) = (i*beta)/(1 - i*beta)    beta = alpha/k
#   t(k) = 1/(1 - i*beta)
#   R + T = 1 always (probability conservation)
#
# Each scattering state is NOT normalizable alone.
# Physical solution = wavepacket integral over k, computed by VEGAS.

import numpy as np
import vegas
import matplotlib.pyplot as plt

alpha = 2.0     # delta potential strength
k0    = 5.0     # central momentum
sigma = 1.0     # momentum spread
x0    = -10.0   # initial packet position (left of delta at x=0)

k_min = k0 - 5*sigma
k_max = k0 + 5*sigma

def phi(k):
    # Gaussian in k-space, shifted to start packet at x0
    # exp(-ik*x0) shifts the packet from x=0 to x=x0 (Fourier shift theorem)
    return np.exp(-(k - k0)**2 / (4*sigma**2)) * np.exp(-1j*k*x0)

def r_coeff(k):
    beta = alpha / k
    return (1j*beta) / (1 - 1j*beta)

def t_coeff(k):
    beta = alpha / k
    return 1 / (1 - 1j*beta)

def compute_probability(x_vals, t_time):
    prob = np.zeros(len(x_vals))
    for i, x in enumerate(x_vals):
        def integrand(k_arr):
            k   = k_arr[0]
            E   = 0.5*k**2
            # Griffiths scattering state -- two regions
            if x < 0:
                psi = np.exp(1j*k*x) + r_coeff(k)*np.exp(-1j*k*x)
            else:
                psi = t_coeff(k) * np.exp(1j*k*x)
            val = phi(k) * psi * np.exp(-1j*E*t_time) / (2*np.pi)
            return [val.real, val.imag]

        integ = vegas.Integrator([[k_min, k_max]])
        integ(integrand, nitn=5, neval=800)
        r2 = integ(integrand, nitn=8, neval=2000)
        prob[i] = r2[0].mean**2 + r2[1].mean**2
    return prob

x_vals = np.linspace(-20, 20, 60)
times  = [0.0, 1.5, 3.0]

E0 = 0.5*k0**2
beta0 = alpha/k0
R = abs((1j*beta0)/(1-1j*beta0))**2
T = abs(1/(1-1j*beta0))**2
print(f"Delta potential  alpha={alpha}, k0={k0}")
print(f"At central momentum: R={R:.3f}, T={T:.3f}, R+T={R+T:.3f}")
print(f"Packet starts at x={x0}, travels right at speed {k0}")
print()

results = {}
for t in times:
    print(f"  t={t}...", flush=True)
    results[t] = compute_probability(x_vals, t)

peak = x_vals[np.argmax(results[0.0])]
print(f"\nt=0 peak at x={peak:.1f}  (expected {x0})")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f"Delta Potential Scattering  α={alpha}, k₀={k0}\n|ψ(x,t)|²",
             fontsize=12, fontweight='bold')
titles = ['t=0  packet approaching',
          't=1.5  packet at delta',
          't=3  split: transmitted (right) + reflected (left)']

for ax, t, title in zip(axes, times, titles):
    ax.plot(x_vals, results[t], 'o-', color='C0', ms=4, lw=1.5, label='|ψ|²')
    ax.axvline(0, color='red', lw=2, ls='--', alpha=0.8, label='δ(x) at x=0')
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('x'); ax.set_ylabel('|ψ|²')
    ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('02_delta_scattering.png', dpi=150, bbox_inches='tight')
plt.show()
