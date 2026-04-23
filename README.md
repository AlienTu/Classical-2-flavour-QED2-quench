# Classical 2-flavour QED2 quench

Browser simulator for the classical two-field quench model.

## Physical model

The simulator evolves two real fields `phi_c(x,t)` and `phi_s(x,t)` with periodic boundary conditions. They are related to the flavour fields `phi_1` and `phi_2` by

```text
phi_c = (phi_1 + phi_2) / sqrt(2),
phi_s = (phi_1 - phi_2) / sqrt(2).
```

So `phi_c` and `phi_s` are the symmetric and antisymmetric superpositions of flavour 1 and flavour 2.

The classical equations of motion are

```text
∂t² phi_c - ∂x² phi_c + g² phi_c
+ (α/2) kappa1 sin[α (phi_c + phi_s)]
+ (α/2) kappa2 sin[α (phi_c - phi_s)] = 0,

∂t² phi_s - ∂x² phi_s
+ (α/2) kappa1 sin[α (phi_c + phi_s)]
- (α/2) kappa2 sin[α (phi_c - phi_s)] = 0,
```

with

```text
α = sqrt(2π).
```

Equivalently, the potential energy density is

```text
V(phi_c, phi_s) = 1/2 g² phi_c²
                  - (kappa1/2) cos[α (phi_c + phi_s)]
                  - (kappa2/2) cos[α (phi_c - phi_s)].
```

The dashed horizontal guide lines in the plot mark

```text
phi = n sqrt(2π),   n ∈ Z,
```

which are useful for visually tracking when `phi_s` crosses neighboring topological sectors.

## Initial condition used in the browser demo

The default quench is

```text
phi_c(x,0) = 0,
phi_s(x,0) = 0,

∂t phi_c(x,0) = A exp[-(x/sigma)²],
∂t phi_s(x,0) = A exp[-(x/sigma)²].
```

So the browser demo starts from the vacuum field profile and injects a local Gaussian velocity kick into both fields.

## Adjustable parameters

- `g` — charge gap
- `kappa1` — flavour-symmetric cosine coupling
- `kappa2` — flavour-antisymmetric cosine coupling
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
