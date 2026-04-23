# Classical 2-flavour QED2 quench

Browser simulator for the classical two-field quench model

- fields: `phi_c`, `phi_s`
- initial condition: `phi_c = phi_s = 0`
- initial velocity: `dot phi_c = dot phi_s = A exp(-(x/sigma)^2)`
- boundary condition: periodic

## Adjustable parameters

- `g` — charge gap
- `kappa1` — heavy cosine coupling
- `kappa2` — light cosine coupling
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
