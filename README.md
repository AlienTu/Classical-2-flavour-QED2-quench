# Classical 2-flavour QED2 quench

This repository contains a browser simulator for the classical time evolution of the bosonized version of 1+1D 2-flavour QED with closed boundary condition.

One may interpret flavour 1 as proton and flavour 2 as electron. The purpose of the quench is to test whether a 1+1D universe can nucleate hydrogen-like objects from a vacuum disturbance. In this language, positive and negative hydrogen are represented by the kink and antikink structure of phi_s.

## Physical model

The simulator evolves two real fields `phi_c(x,t)` and `phi_s(x,t)` with periodic boundary conditions. They are related to the flavour fields `phi_1` and `phi_2` by

```text
phi_c = (phi_1 + phi_2) / sqrt(2),
phi_s = (phi_1 - phi_2) / sqrt(2).
```

The classical equations of motion are

```text
∂t² phi_c - ∂x² phi_c + g² phi_c
+ (sqrt(2π)/2) kappa1 sin[sqrt(2π) (phi_c + phi_s)]
+ (sqrt(2π)/2) kappa2 sin[sqrt(2π) (phi_c - phi_s)] = 0,

∂t² phi_s - ∂x² phi_s
+ (sqrt(2π)/2) kappa1 sin[sqrt(2π) (phi_c + phi_s)]
- (sqrt(2π)/2) kappa2 sin[sqrt(2π) (phi_c - phi_s)] = 0,
```

and the potential energy density is

```text
V(phi_c, phi_s) = 1/2 g² phi_c²
                  - (kappa1/2) cos[sqrt(2π) (phi_c + phi_s)]
                  - (kappa2/2) cos[sqrt(2π) (phi_c - phi_s)].
```

## Initial condition used in the browser demo

The default quench is

```text
phi_c(x,0) = 0,
phi_s(x,0) = 0,

∂t phi_c(x,0) = A exp[-(x/sigma)²],
∂t phi_s(x,0) = A exp[-(x/sigma)²].
```

In the flavour basis this corresponds to a local Gaussian momentum kick in flavour 1, while flavour 2 starts from vacuum.

## Adjustable parameters

- `g` — charge gap
- `kappa1` — flavour-1 cosine coupling
- `kappa2` — flavour-2 cosine coupling
- `L` — box size
- `N` — grid points
- `dt` — time step
- `A` — Gaussian kick amplitude
- `sigma` — Gaussian kick width

Each parameter supports both slider input and direct number entry.

## Default parameters

- `g = 1`
- `kappa1 = 5`
- `kappa2 = 0.5`
- `L = 200`
- `N = 2000`
- `dt = 0.02`
- `A = 4`
- `sigma = 4`

## Run with GitHub Pages

If Pages is enabled for the repository, the site should be served from the root branch (`main`). The main entry point is `index.html`.
