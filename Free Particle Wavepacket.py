# hbar=1,m=1 for simplicity
# The free particle stationary states are e^(ikx) -- plane waves.
# Each one is NOT normalizable: integral of |e^ikx|^2 dx = infinity.
# A single momentum state cannot represent a real particle.
# For this we build a wavepacket superposing many k values
# psi(x,t) = (1/2pi) * integral phi(k) * e^(ikx - ik^2*t/2) dk
# We choose a Gaussian centered at k0
# This integral IS normalizable and represents a real localized particle.
# VEGAS computes this integral numerically.

import vegas
import matplotlib.pyplot as plt
import numpy as np
k0    = 5.0    # central momentum packet moves right at speed k0
sigma = 1.0    # spread in momentum space
def phi(k):
  return np.exp(-(k - k0)**2 / (4*sigma**2))

k_min = k0 - 5*sigma    # phi is essentially zero beyond 5 sigma
k_max = k0 + 5*sigma

 compute_probability(x_vals, t):
    prob = np.zeros(len(x_vals))
    for i, x in enumerate(x_vals):
        def integrand(k_arr):
            k     = k_arr[0]
            phase = k*x - 0.5*k**2*t    # Griffiths phase: kx - hbar*k^2*t/2m
            f     = phi(k) / (2*np.pi)
            # e^(i*phase) = cos(phase) + i*sin(phase)
            # VEGAS only handles real numbers so split into two parts
            return [f*np.cos(phase), f*np.sin(phase)]
 
        integ = vegas.Integrator([[k_min, k_max]])
        integ(integrand, nitn=5, neval=800)     # warmup: VEGAS learns integrand shape
        r = integ(integrand, nitn=8, neval=2000) # production run
        prob[i] = r[0].mean**2 + r[1].mean**2   # |psi|^2 = Re^2 + Im^2
    return prob
 
def analytical(x, t):
    # Exact result from completing the square in the Fourier integral
    alpha   = 1/(4*sigma**2) + 1j*t/2
    prefac  = 1/(4*np.pi**2) * np.pi/np.abs(alpha)
    return prefac * np.exp(-2*(x - k0*t)**2 * np.real(1/(4*alpha)))
 
x_vals = np.linspace(-3, 20, 50)
times  = [0.0, 1.0, 2.5]
 
print("Free particle wavepacket (k0=5, sigma=1)")
print("Packet moves right at speed k0, spreads with time")
print()
results = {}
for t in times:
    print(f"  t={t}...", flush=True)
    results[t] = compute_probability(x_vals, t)
 
# Verify at t=0: peak should be at x=0
peak = x_vals[np.argmax(results[0.0])]
print(f"\nt=0 peak at x={peak:.2f}  (expected 0.00)")
 
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle("Free Particle Wavepacket  |ψ(x,t)|²", fontsize=13, fontweight='bold')
 
for ax, t in zip(axes, times):
    ax.plot(x_vals, results[t], 'o-', color='C0', ms=4, lw=1.5, label='VEGAS')
    x_sm = np.linspace(-3, 20, 300)
    ax.plot(x_sm, analytical(x_sm, t), 'r--', lw=2, label='Analytical')
    ax.set_title(f't = {t}', fontsize=11)
    ax.set_xlabel('x'); ax.set_ylabel('|ψ(x,t)|²')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
 
plt.tight_layout()
plt.savefig('01_free_particle_wavepacket.png', dpi=150, bbox_inches='tight')
plt.show() 
