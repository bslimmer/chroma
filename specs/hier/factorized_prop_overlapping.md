# Overlapping-Boundary Factorized Perambulator (draft)

Status: implemented 2026-10-16, not yet compiled or run (decisions at the end).
Target files: `lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc`,
plus one field in `Params::NamedObject_t` in the matching `.h`.
References: Ce, Giusti, Schaefer, arXiv:1601.04587, Sec. 2 and 5 (Eqs. 2.4, 5.1-5.3);
`specs/hier/factorized_prop.md`; Slimmer notes "Notes on Peardon et al. & factorized
perambulators" (Eqs. 9-17).

## Goal

Add an XML switch that makes the existing parent-lattice factorized test path build the
cross-subdomain perambulator with the overlapping-domain (Giusti) factorization instead
of the current non-overlapping Schur-complement factorization. Everything else
(sources, phasings, projector vectors, output format) stays the same.

## Current (non-overlapping) path, for reference

Field order (Lambda0 = frozen, Lambda1 = active). The code computes

```text
tau_fact(t',t) = V^dag(t') D11^-1 Pi D10 D00^-1 D01 D11^-1 V(t)
```

- `D11`, `D00` come from one solver on `u_sep` (t links zeroed on `t_start` and
  `t_end - 1` of each frozen interval, so `D00` sits one slice inside the frozen block).
- `D01`, `D10` are the full `fLinOp` (uncut `u`) followed by `restrictToTimeslices`.
- `Pi` is the per-spin Laplacian-mode projector (`num_intermediate_projectors` modes,
  optional `intermediate_projector_phase`) on the first and last slice of the sink block.
- No explicit sign: the two minus signs of the block inverse cancel.
- Three solve sets per (t_source, phasing, spin): `D11^-1`, `D00^-1`, `D11^-1`.

## Geometry (parent-lattice frame)

Same assumptions as the current path: two frozen intervals read from
`plan.child0_frozen_local_intervals`, used directly as parent-global slices. This needs
`cut0 = 0` (true for the 32^3x64 test: `cut0=0, cut1=32`).

```text
F0 = [0, fw-1]            sink side below (via wrap Lt-1 -> 0), source side above
F1 = [c1, c1+fw-1]        source side below, sink side above      (c1 = cut1)
A_src = [fw, c1-1]        source active region (t_source must lie here)
A_snk = [c1+fw, Lt-1]     sink active region
```

### Symmetric layout (revised 2026-10-21)

Each frozen interval is split into one source-only buffer slice, `fw-2` overlap slices,
and one sink-only buffer slice. The projector sits on the source-only buffer slice
(`w2`), and the hop goes from there into the first overlap slice (`w1`):

| interval | source-only (w2, projector) | overlap | sink-only | hop link zeroed for D_snk |
|---|---|---|---|---|
| F0 | `fw-1` | `[1, fw-2]` | `0` | rooted on `fw-2` |
| F1 | `c1`   | `[c1+1, c1+fw-2]` | `c1+fw-1` | rooted on `c1` |

Hop targets `w1`: `fw-2` (F0) and `c1+1` (F1). Example c1 = 32, fw = 4: projector on 3
and 32, overlap `{1, 2}` and `{33, 34}`, sink-only buffers 0 and 35. Growing fw grows
only the overlap.

Why the buffer slices: the planned fix builds `D_snk` as the exact restriction of `D`
(clover from the uncut field, only the hop dropped). Its clover term on `w1` uses the
link `w2 -> w1` and the spatial links on `w2`, i.e. one slice before `Omega_snk`, which
is the frozen source-only buffer slice. This only holds for an unsmeared Dirac operator:
stout-type smearing reaches one more slice, onto dynamical source links. The code
therefore requires `SIMPLE_FERM_STATE` in overlapping mode. The current cut-field
operators are local regardless; the check protects the exact-`D_snk` fix. `D_src` does
not need the sink-only buffer (it is the approximated piece and its cut-field version is
local); the buffer keeps the layout symmetric.

### Domains

```text
Omega_src = [1, c1+fw-2]   (A_src + all of F0, F1 except the sink-only buffer slices)
            Dirichlet: zero t links rooted on 0 and on c1+fw-2
Omega_snk = [c1+1, Lt-1] u [0, fw-2]   (A_snk + all of F0, F1 except the source-only
            buffer slices)
            Dirichlet: zero t links rooted on c1 and on fw-2   (the hop links)
```

Two separate cut fields, two solvers. Example c1=32, fw=8, Lt=64:

| field | zeroed t links (rooted on) | links cut | domain used |
|---|---|---|---|
| `u_src` | 38, 0 | F1 overlap -> sink-only buffer (38->39), F0 sink-only buffer -> overlap (0->1) | `[1, 38]` |
| `u_snk` | 32, 6 | hop in F1 (32->33), hop in F0 (6->7) | `[33, 63] u [0, 6]` |

The antiperiodic sign on 63->0 lies only in `Omega_snk`. Each cut field also has a
second, disconnected block (`[39, 63] u [0]` for `u_src`, `[7, 32]` for `u_snk`); sources
never touch it, so its solution stays zero.

The two domains overlap on `[c1+1, c1+fw-2]` and `[1, fw-2]`, i.e. `fw-2` slices per
frozen interval.

## Factorization

Exact identity (block inverse with `Omega_snk` vs. its complement, Giusti Eq. 2.4/4.9),
for `y in Omega_snk`, `x` outside it:

```text
D^-1(y,x) = - sum_{w1,w2} D_snk^-1(y,w1) D(w1,w2) D^-1(w2,x)
```

`w1 in {fw-2, c1+1}`, `w2 in {fw-1, c1}`; `D(w1,w2)` is only the temporal hop
(`P_+- U_t` term), so no clover term enters it.

Approximation (Giusti Eq. 5.1): `D^-1(w2,x) -> D_src^-1(w2,x)`. `Omega_src` extends
`d = fw-2` slices beyond each `w2` (the overlap thickness), so the error is expected
`~ exp(-M_pi (fw-2))`, the same as the previous (asymmetric) layout.

Projected form (Giusti Eq. 5.3), using the same distillation projector as now, placed
on the `w2` slices:

```text
tau_ov(t',t) = - V^dag(t') D_snk^-1 D_{w1<-w2} Pi_{w2} D_src^-1 V(t)

Pi_{w2} = sum_{t in w2} sum_{m < Nproj} phi_m(t) phi_m^dag(t) (x) 1_spin
```

`phi_m(t)` = the first `num_intermediate_projectors` Laplacian eigenvectors on slice
`t`, with `intermediate_projector_phase`, exactly as in the current code.

Left/right pieces (for the later baryon-tensor work; notes Eq. 17 analogue):

```text
L~(t')^{i,m} = - V^dag(t') D_snk^-1 D_{w1<-w2} phi_m
R~(t)^{m,j}  =   phi_m^dag D_src^-1 V(t)
```

## Algorithm per (t_source, phasing, spin_source)

1. `psi1 = D_src^-1 V(t_source)` with `PP_src` (`doInversion`, all Lt slices, as now).
2. `psi1 <- Pi_{w2} psi1`: restrict to the two `w2` slices and project. Reuse the existing
   projection code (see Implementation).
3. `psi2 = D psi1` with `fLinOp` (uncut `u`), then `restrictToTimeslices(psi2, w1)`.
4. `psi2 <- -psi2` (the single block-inverse minus sign).
5. `psi3 = D_snk^-1 psi2` with `PP_snk`.
6. `quark_solns = toSBTensor(psi3, 0, Lt)`; sink contraction and storage unchanged.
   Slices outside `Omega_snk` come out zero.

Two solve sets per (t_source, phasing, spin) instead of three; no `D00` solve.

## XML

New optional element in `NamedObject`, next to the other factorized options:

```xml
<overlapping_factorization>true</overlapping_factorization>   <!-- default false -->
```

- `false` / absent: current non-overlapping path, bit-for-bit unchanged.
- `true`: requires the factorized path to be enabled (sidecar + child id given);
  otherwise throw.
- Reuses `num_intermediate_projectors` and `intermediate_projector_phase`.
- In overlapping mode, `num_intermediate_projectors <= 0` (default -1) means no
  projector (`Pi = 1`): step 2 is only the restriction to the `w2` slices, and the
  projector colorvecs are not loaded. The non-overlapping path keeps its current
  behaviour.
- Echo to output XML: the mode, `fw`, the `w1`/`w2` slices, and the zeroed-link slices
  for each domain.

## Implementation notes (all in the .cc)

- `read`/`write` for `NamedObject_t`: optional bool, default `false`.
- In `if (factorized_geometry.enabled)`: branch on the flag only where the paths differ
  (gauge-field cuts / solver construction, and the inner chain). Source/sink colorvec
  loading, phasing loops and storage stay shared.
- Build `u_src`, `u_snk` with the existing `zeroTemporalLinksOnSlice`; two
  `SB::ChimeraSolver`s (`PP_src`, `PP_snk`). Both stay alive for the whole loop (two MG
  setups in memory).
- Move the existing per-slice projection block into a file-local helper, e.g.
  `projectOnSlices(fields, slices, colorvec_proj)`, and call it from both paths. For the
  non-overlapping path it is called with the same two sink-edge slices as now, so the
  result is identical.
- Geometry checks (throw with a clear message): 2 frozen intervals, equal widths,
  `fw >= 3` (at least one overlap slice), `F0` starts at 0, every `t_source` in `A_src`.
- Fermion state check: `FermionAction/FermState/Name` must be `SIMPLE_FERM_STATE`
  (absent counts as simple, as in Chroma's reader); otherwise throw.
- Do not add new full-field `.print()` debug dumps.
- Header: add `bool overlapping_factorization;` to `Params::NamedObject_t`.
- No changes to `factorized_distillation.{h,cc}`.

## Known caveats

- `D_snk` exactness (not fixed by the buffer slice). The exact identity needs `D_snk`
  to be `D` restricted to `Omega_snk`, with the clover term on `w1` built from the full
  field. `PP_snk` is built from `u_snk`, where the hop link is zeroed before the clover
  term is built, so the clover term on `w1` differs from the true restriction. The
  `D_snk` solve is sourced exactly on `w1`, so this is an O(1) local error, not
  suppressed by the overlap (toy: ~0.97 per quark for clover only). The buffer slice
  keeps the exact fix local; it does not remove this error. The fix (clover from the
  uncut field, hop-only cut) needs an operator the stock `ChimeraSolver` does not
  provide. Validation 2 measures its size.
- Expected effect on the signal deficit: the project review (2026-10-14) found that the
  Laplacian projector, not the `S^-1 -> D00^-1` step, dominates the loss (~0.4 per
  quark at Nproj = 128). This mode replaces the Schur step and moves the projector to
  the source-only buffer slice; it is not expected on its own to remove the projector loss.

## Validation

1. Flag off: output identical to the current code on one config (diff the SDB).
2. Exact-identity check (dev only): replace `D_src^-1` by the full `D^-1` (solver on
   `u`) and skip `Pi`. `tau_ov` must equal the exact perambulator on `Omega_snk` slices
   to solver precision, up to the caveat above. This tests sign, hop restriction and
   `D_snk`.
3. Approximation check: `num_intermediate_projectors <= 0` (`Pi = 1`), `D_src^-1` as specified. Ratio to the exact
   correlator vs fw (4, 8); expect it to approach 1 like `exp(-M_pi (fw-2))`.
4. Production comparison: `<C_ov>/<C_true>` vs the non-overlapping path at the same
   Nproj (64, 128) and fw (4, 8), 100 cfgs, same analysis scripts.

## Decisions

1. (Revised 2026-10-21, replaces the 2026-10-16 layouts) Symmetric layout: one
   source-only buffer slice (projector, `c1` and `fw-1`), `fw-2` overlap slices, one
   sink-only buffer slice. Hop from the projector slice into the first overlap slice.
   `d` = `fw-2`. Requires `SIMPLE_FERM_STATE` (enforced).
2. Projector on the source side of the hop (`w2`), as in Giusti Eq. 5.3. Any `fw >= 3`.
3. The flag is a new `NamedObject_t` member (one line in the `.h`).
4. Overlapping mode supports `num_intermediate_projectors <= 0` as "no projector" for
   validation.

## Open

- Whether to free `PP_src` before building `PP_snk` if two MG setups do not fit in
  memory (would need the step-2 results for all spins/vectors buffered first).
- Exact-identity check (Validation 2) is a temporary dev edit, not an XML option.
