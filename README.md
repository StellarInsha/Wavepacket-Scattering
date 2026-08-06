# Quantum Scattering

*A numerical exploration of free propagation, reflection, transmission, and finite-well scattering using VEGAS Monte Carlo integration.*

A plane wave extends forever.

It has a perfectly defined momentum, exists everywhere at once, and cannot even be normalized over all space.



A real particle is not quite like that.

So instead of treating the stationary solutions of the Schrödinger equation as the end, this project explores beyond:

> **What happens if we actually build a localized particle and watch it move?**

The answer starts with a Gaussian wavepacket and gets progressively more and more interesting.

---

## What this repository does

Following the scattering problems in **Griffiths' *Introduction to Quantum Mechanics*, Chapter 2**, this repository numerically constructs time-dependent wavepackets for three cases:

### 01 — Free Particle

We begin with the simplest possible situation: no potential at all.

Individual free-particle solutions,

[
e^{ikx},
]

aren't normalizable. Instead, many momentum states are superposed with a Gaussian weighting:

[
\psi(x,t)=\frac{1}{2\pi}\int \phi(k)
e^{i(kx-k^2t/2)},dk.
]

The result is a localized packet that **moves and spreads with time**.

This first case also gives us something valuable: an analytical solution.

That means the numerical VEGAS result can be plotted directly against the exact answer—a sanity check before moving on to actual scattering.

---

### 02 — Throw in a Delta Potential

Next, the packet encounters

[
V(x)=-\alpha\delta(x).
]

Now something distinctly quantum happens.

The incoming packet doesn't simply cross the potential or bounce away from it. It develops **reflected and transmitted components simultaneously**.

For every momentum (k), the corresponding stationary scattering state contains

[
\text{incoming}+\text{reflected}
]

on one side and

[
\text{transmitted}
]

on the other.

Those states are again individually non-normalizable, so the code superposes them into a physical wavepacket.

The numerical evolution lets us watch an initially localized packet approach (x=0), interact with the delta potential, and eventually separate into reflected and transmitted pieces.

As a check,

[
R+T=1,
]

so probability remains conserved.

---

### 03 — A Finite Square Well

Finally, things become a little less friendly.

Consider

[
V(x)=
\begin{cases}
-V_0,& |x|<a\
0,& |x|>a.
\end{cases}
]

Instead of two spatial regions, we now have three:

* before the well,
* inside the well,
* after the well.

The wavefunction and its derivative must match at **both boundaries**, and the particle carries a different wavevector inside the well because its kinetic energy changes there.

The result is interference between waves travelling through the well.

And hidden inside that interference is one of my favourite results from this simulation:

### Ramsauer–Townsend transmission

At certain energies,

[
\sin(2\ell a)=0,
]

and the well becomes **perfectly transparent**:

[
T=1.
]

The particle encounters a non-zero potential and yet passes through with no reflection.

Not because the potential is weak.

Because the phases line up.

---

## Why VEGAS?

These problems have analytical machinery behind them.

Each physical wavepacket requires evaluating an integral over momentum:

[
\psi(x,t)=\int dk;
\phi(k)\psi_k(x)e^{-iE_kt}.
]

Instead of performing that integral symbolically, this repository evaluates it numerically using **VEGAS adaptive Monte Carlo integration**.

Because VEGAS operates on real-valued integrands, the complex wavefunction is separated into its real and imaginary components, integrated independently, and reconstructed through

[
|\psi|^2=(\Re\psi)^2+(\Im\psi)^2.
]

So the project also became a small experiment in translating textbook quantum mechanics into numerical computation.

---

## The progression

The three simulations are deliberately ordered:

```text
Free particle
     ↓
Does the numerical wavepacket behave correctly?

Delta potential
     ↓
Can it reproduce reflection + transmission?

Finite square well
     ↓
What happens when boundaries and interference enter the picture?
```



---

## Tools

* **Python**
* **NumPy** — numerical operations
* **VEGAS** — adaptive Monte Carlo integration
* **Matplotlib** — visualization

Natural units are used throughout:

[
\hbar = 1,\qquad m=1.
]



Each script generates the probability density

[
|\psi(x,t)|^2
]

at several moments in time and saves the corresponding visualization.

---

## What I was trying to understand

This repository started  as an attempt to answer something that bothered me while reading scattering theory:

**If stationary scattering states extend across all space, what would an actual particle approaching the potential look like?**

Wavepackets make that distinction visible.

Instead of only calculating (R) and (T), we get to watch where those probabilities come from.

A packet approaches.

It interferes.

It spreads.

Part of it returns.

Part of it continues.

And, under just the right conditions, a potential well that should seemingly disturb the particle becomes completely transparent.

That's the part of quantum mechanics I wanted this repository to capture.

---

## Reference

D. J. Griffiths & D. F. Schroeter,
*Introduction to Quantum Mechanics*

Particularly the Chapter 2 discussions of:

* Free particles
* The delta-function potential
* Scattering states
* The finite square well

---

> **The equations tell us the transmission probability.
> The wavepacket lets us watch it happen.**
