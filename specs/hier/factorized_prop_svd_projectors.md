# Edge-SVD Boundary Projectors for the Factorized Perambulator

Status: implemented 2026-10-22 (working tree, not committed). The .cc compiles to an
object file (g++ / OpenMPI, QDP++ `eloy/localbinarydb` scalar, superbblas master,
`BUILD_SB`). The helper routines were run standalone on a 4^3x8 scalar lattice (see
Verification). It has not yet been run inside Chroma with real solvers.
Target files: `lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc`,
plus new fields in `Params::NamedObject_t` in the matching `.h`.
References: `specs/hier/factorized_prop.md`, `specs/hier/factorized_prop_overlapping.md`;
Ce, Giusti, Schaefer, arXiv:1601.04587 (D-adapted boundary vectors);
Halko, Martinsson, Tropp, arXiv:0909.4061 (randomized range finder);
toy comparison `fact_prop_testing/toy_overlap_vs_nonoverlap/basis_cmp.py`.

## Goal

Add an XML option that replaces the vectors used to build the boundary projector `Pi`.
Today `Pi` is a per-spin projector onto Laplacian eigenvectors (`distillation_vecs`).
The new option (`singular_vectors`) uses the dominant singular vectors of the
edge-to-edge block that carries the quark across the frozen region. These are
spin-colour vectors, and they are computed once per run from frozen-region operators only.

Only the generation of the projection vectors changes, plus the projection helper
that applies them. The chain, cut fields, production solvers, signs, slices, phasings,
and output format are unchanged.

Toy motivation (4^3x16, clover, rough field, per-quark ratio vs. the same chain with
`Pi = 1`, rank R on a 768-dim slice): at R = 32, Laplacian 0.37 / 0.41, edge-SVD
0.71 / 0.77 (non-overlap / overlap). The frozen-only overlap block `G_F` matched the
true-`D_snk` block `G` to ~0.01.

## XML

New optional elements in `NamedObject`, next to the other factorized options:

```xml
<projector_type>singular_vectors</projector_type>  <!-- distillation_vecs (default) | singular_vectors -->
<num_intermediate_projectors>128</num_intermediate_projectors>  <!-- existing element -->
<svd_oversampling>16</svd_oversampling>            <!-- optional, default 16 -->
<svd_power_iterations>2</svd_power_iterations>     <!-- optional, default 2 -->
<svd_seed>11</svd_seed>                             <!-- optional, default 11 -->
```

- `projector_type` absent or `distillation_vecs`: current behaviour, bit-for-bit
  unchanged, in both paths. Any other value: throw. The four new elements are only
  echoed in the output XML for the non-default type.
- `singular_vectors`:
  - `num_intermediate_projectors` = R = the number of spin-colour singular vectors per
    projector slice (total rank of `Pi` on that slice). For comparisons at equal
    rank with `distillation_vecs` (rank 4N per slice), use R = 4N.
  - Non-overlapping path: R must be > 0, otherwise throw.
  - Overlapping path: R <= 0 keeps its current meaning (`Pi = 1`, no vectors are computed).
  - `intermediate_projector_phase` must be all zero, otherwise throw (phasing is a
    property of the Laplacian basis only).
  - R must be <= 6 V3 (the spin half that the temporal hop passes); otherwise throw.
- Output XML group `EdgeSVDProjectors`: rank, oversampling, power iterations, seed,
  time, and per projector slice (`elem`): `projector_slice`, `far_slice`, all k Ritz
  `singular_values`, and `orthonormality_error` (max |phi^dag phi - 1|).

## Geometry and notation

Same assumptions as the current code: two frozen intervals from
`plan.child0_frozen_local_intervals`, `F0 = [0, fw-1]`, `F1 = [c1, c1+fw-1]`.
`R_t` = restriction to slice(s) t (`restrictToTimeslices`). `D` = `fLinOp` (uncut `u`),
`D^dag` = `fLinOp` with `MINUS`. `g5` = `Gamma(Ns*Ns-1)`. `S_X = X^-1` via
`SB::doInversion` with solver X. Inverse adjoints use g5-hermiticity:
`(X^-1)^dag = g5 X^-1 g5`.

For every case below, the basis is the top-R **right** singular vectors of an operator
`M : slice p -> slice q`, where `p` is the projector slice.

### Non-overlapping path

Projector slices (as now): `p = c1+fw-1` (F1) and `p = 0` (F0). `D00` slabs:
`[c1+1, c1+fw-2]` and `[1, fw-2]` (from `PP_sep`). The quark-line block is
`T = R_p D R_{e_p} S_00 R_{e_s} D R_q`, with `Pi` acting on its output, so the basis
is the left singular vectors of `T` = right singular vectors of `M = T^dag`:

```text
M   y = R_q D^dag R_{e_s} g5 S_00 g5 R_{e_p} D^dag y      (y on p)
M^dag x = R_p D R_{e_p} S_00 R_{e_s} D x                  (x on q)
```

| interval | p (Pi) | e_p | e_s | q (source buffer) | solver |
|---|---|---|---|---|---|
| F1 | c1+fw-1 | c1+fw-2 | c1+1 | c1 | `PP_sep` (existing) |
| F0 | 0 | 1 | fw-2 | fw-1 | `PP_sep` (existing) |

Example c1 = 32, fw = 4: F1 `p=35, e_p=34, e_s=33, q=32`; F0 `p=0, e_p=1, e_s=2, q=3`.
Requires `F0` to start at 0 and both widths >= 3 (checked).

### Overlapping path

Projector slices (as now): `w2 = c1` (F1) and `w2 = fw-1` (F0); hop targets
`w1 = c1+1`, `fw-2`. The block uses a new frozen-only operator `D_F`:

```text
u_F = u with t links zeroed on slices {c1, c1+fw-1, fw-2, Lt-1}
      -> Dirichlet slabs [c1+1, c1+fw-1] and [0, fw-2]
PP_F = ChimeraSolver(fermact, invParam, u_F)

M   x = R_q S_F R_{w1} D x                    (x on p = w2)
M^dag y = R_p D^dag R_{w1} g5 S_F g5 R_q y      (y on q)
```

| interval | p = w2 (Pi) | w1 | q (sink-only buffer) | solver |
|---|---|---|---|---|
| F1 | c1 | c1+1 | c1+fw-1 | `PP_F` (new, temporary) |
| F0 | fw-1 | fw-2 | 0 | `PP_F` (new, temporary) |

Example c1 = 32, fw = 4: F1 `p=32, w1=33, q=35`; F0 `p=3, w1=2, q=0`.

`D_F` is not `D_snk`. Using `PP_snk` would make the basis depend on sink-active links,
which is not valid for two-level. `PP_F` is built only when needed, and it is destroyed
after the basis is computed, so no third MG setup stays alive during the production loop.
While the basis is being computed, three MG setups (`PP_src`, `PP_snk`, `PP_F`) are in memory.

### Both intervals in one batch

`M` is block-diagonal over the two intervals: their p, q and slab slices are disjoint and
the slabs are decoupled. One application of `M` (or `M^dag`) to a vector supported on
both p slices (or both q slices) applies both blocks at once. All solves are therefore
shared between the two intervals. Only the small dense linear algebra is done
separately per slice. If the two intervals end up with different column counts (only
after rank truncation), the shorter set repeats earlier columns and those outputs are
dropped, so no solve ever gets a zero source.

## Algorithm (randomized SVD, fixed cost)

k = R + `svd_oversampling`, qp = `svd_power_iterations`.
"Per slice" means separately for each of the two p (or q) slices.

1. Save the QDP RNG state, seed with `svd_seed`. Draw k Gaussian `LatticeFermion`s and
   keep their values on the two q slices: `Omega`. Restore the RNG state.
2. `Y = M^dag Omega` (lands on the p slices).
3. Repeat qp times: `Y <- orth(Y)`; `Y <- M^dag (M Y)`.
4. `Q = orth(Y)` per p slice (k orthonormal columns, or fewer if rank-deficient;
   throw if fewer than R).
5. `Z = M Q` (lands on the q slices).
6. Per slice: `C = Z^dag Z` (k x k), `SB::svd(C) -> W, lambda` (descending).
   Singular values `sigma = sqrt(lambda)`. Basis on p: `phi = Q W[:, 0:R]`.

`orth(Y)` per slice (Cholesky-QR via the Gram matrix, done twice for stability):
`C = Y^dag Y`, `SB::svd(C) -> W, lambda`, keep `lambda_i > 100 eps lambda_0`
(eps of the build precision), `Y <- Y W_keep diag(lambda_keep^-1/2)`; repeat once.

Cost: k (2 + 2 qp) applications of `M` / `M^dag`, each one slab solve (shared by both
intervals). Example: R = 128, k = 144, qp = 2 gives 864 solves. The slab solves run on the
full lattice (the cut operator is block-diagonal), but they are well conditioned.
`M` is applied in batches of `max_rhs` columns using the existing multi-RHS
`SB::doInversion(PP, out, in, max_rhs)`.

Memory: `Omega`, `Y`, `Q`, `Z` are stored as slice tensors (`toSBTensor(psi, t, 1)`, order
`cxyztXns`, t extent 1, columns on `n`), not as `LatticeFermion` vectors. Full
`LatticeFermion`s exist only for the batch being solved (at 32^3x64 one
`LatticeFermion` is ~0.4 GB; one slice column is ~6 MB in double).

Default qp = 2 (changed from 1 in the draft). Toy check of the randomized basis against
the exact SVD basis (per-quark ratio relative to `Pi = 1`, rough 4^3x16 clover field,
mean of 2 seeds, oversampling 16):

| path, R | exact SVD | qp = 0 | qp = 1 | qp = 2 | qp = 3 |
|---|---|---|---|---|---|
| non-overlap, 8 | 0.407 | 0.190 | 0.348 | 0.383 | 0.401 |
| non-overlap, 32 | 0.713 | 0.461 | 0.699 | 0.713 | 0.714 |
| non-overlap, 64 | 0.844 | 0.592 | 0.815 | 0.845 | 0.849 |
| overlap, 8 | 0.461 | 0.302 | 0.458 | 0.461 | 0.461 |
| overlap, 32 | 0.767 | 0.577 | 0.764 | 0.770 | 0.768 |
| overlap, 64 | 0.889 | 0.746 | 0.879 | 0.888 | 0.889 |

## Existing Chroma / superbblas pieces used (nothing reimplemented)

- `SB::svd<Nu,Ns,Nv>(x, rows, cols, label)` (superb_contractions.h; used in mgproton for
  null-vector blocking): all dense SVDs, which here are only the small k x k Gram
  matrices. superbblas factors each matrix once and the result is then replicated
  (`make_sure(none, OnHost, OnEveryoneReplicated)`), so every rank uses identical vectors.
  The small coefficient matrices are copied to a tensor with the distribution of the Gram
  matrix before contracting with the slice tensors (the mgproton
  `getDeflationProj` pattern).
- `SB::contract`: Gram matrices `Y^dag Y` over `cxyztXs` (does the MPI reduction) and
  products `Y W`.
- `SB::doInversion` (multi-RHS `MultipleLatticeFermions` form), `fLinOp` with
  `PLUS` / `MINUS`, `restrictToTimeslices`, `toSBTensor`, `returnNLatticeFermions`,
  `zeroTemporalLinksOnSlice`: all applications of `M`.
- Not used: PRIMME via `SB::getInexactEigensolverGD` on `M^dag M`. It would work, but
  its matvec count is not fixed, it is harder to batch, and its operator domain is the
  full lattice. It stays a fallback if the randomized basis turns out to be inaccurate.

## Implementation (as written)

- `.h`: `std::string projector_type; int svd_oversampling; int svd_power_iterations;
  int svd_seed;` in `NamedObject_t`; `read` / `write` with the defaults above.
- File-local helpers in the `.cc`, next to `projectOnSlices`:
  - `projectOnSlicesSV(psi, slices, basis)`: `out|_t = phi_t (phi_t^dag psi|_t)`,
    contracting `c,x,y,z,X,s`; zero on other slices. `basis` is
    `std::map<int, SB::Tensor<Nd+4, SB::Complex>>` (slice -> `cxyztX s n`, t extent 1).
  - `computeEdgeSVDBasis(M, Mdag, p_slices, q_slices, R, oversampling, power_its, seed,
    batch, sing_vals, ortho_err)`: the algorithm above.
  - `edgeSVDApply`, `edgeSVDOrthonormalize`, `edgeSVDGram`, `edgeSVDGramEig`,
    `edgeSVDTimesSmall`, `newEdgeSVDSliceTensor`: batching, Cholesky-QR, Gram, SVD and
    small-matrix products.
- In `func`: `sv_projector` flag and the R / phase / 6V3 checks next to the solver
  setup; the basis is computed right after `fLinOp` is built and before the `t_source`
  loop, with `M` / `M^dag` lambdas built per path from the formulas above. `PP_F` is
  freed right after. `colorvec_proj` is not loaded in `singular_vectors` mode. The two
  existing projector call sites call `projectOnSlicesSV` instead of `projectOnSlices`
  when `sv_projector` is set.
- No changes to `factorized_distillation.{h,cc}`. No new full-field `.print()` dumps.

## Known caveats

- `D_F` is built like the other cut operators (links zeroed before the clover term), so
  on its edge slices it is not the exact restriction of `D`. This is fine for a basis
  (it only has to span the right directions). The t links rooted on c1+fw-1 and Lt-1
  are zeroed, so the clover leaves on the slab edges that would reach the active
  slices c1+fw and Lt-1 vanish. On the slab, `D_F` therefore depends only on frozen links.
- Non-overlapping path with stout smearing (not `SIMPLE_FERM_STATE`): `fLinOp` hops use
  smeared links, so the basis is not strictly frozen-only. This doesn't matter for
  full-config tests, but it does for two-level.
- This is a deterministic truncation, so a bias remains. Check it against the
  `Pi = 1` overlap run (Validation 3).
- The basis is per configuration. In this measurement it is recomputed every run; it is not
  stored (see Open).

## Verification done (2026-10-22)

1. Compile: the modified `.cc` compiles to an object file with `-Wall` (only the
   sign-compare warnings the file already had), against QDP++ `eloy/localbinarydb`
   (scalar), superbblas master and OpenMPI, with `BUILD_SB` on.
2. Standalone runtime test of the helpers (`t_edge_svd.cc`, scalar QDP, 4^3x8): a known
   block-diagonal operator on two slice pairs with random site spin and colour matrices.
   - Adjoint consistency of the test operator ~1e-15 relative.
   - Exact reference: `SB::svd` of the dense block. With a full-rank sketch (k = 768),
     singular values agree to 1e-14, eigen-residuals are 1e-14, and
     |phi^dag phi - 1| is 1e-14.
   - With k = R + 16 on this flat spectrum, the top-R singular values are within ~5% of exact.
   - `projectOnSlicesSV`: Pi Pi = Pi to 1e-27, zero off the projector slices,
     Pi phi_0 = phi_0.
3. Toy (Python, the toy Wilson-clover operators): the C++ formulas for `M` and `M^dag` in
   both paths are adjoint to 1e-15. The randomized basis reproduces `Pi = 1` exactly at full
   rank (R = 6 V3), and matches the exact-SVD basis at the default settings (table above).

Not verified: anything multi-rank (the standalone test is single-process); the real
`ChimeraSolver` slab solves; GPU builds.

## Validation still to do

1. `projector_type` absent: SDB identical to the current code on one config, both paths.
2. Full-rank check (small lattice, e.g. 8^3x16): with R = 6 V3 (the whole spin half that
   the hop passes), `singular_vectors` must reproduce the `Pi = 1` result to solver
   precision. For the overlap, compare with `num_intermediate_projectors <= 0`; for the
   non-overlap, compare with a dev-only run without the projector. This tests the slices,
   the adjoints, the g5 handling and the orthonormalization with real solvers and MPI.
3. On one 32^3 config, fw = 4 and 8: singular-value decay per slice, and the plateau
   `<C_fact>/<C_true>` vs R (64, 128, 256) against `distillation_vecs` at equal total rank
   (N = R/4) and against overlap `Pi = 1`.
4. Production: 100 cfgs at the chosen R, same analysis scripts as `128_bw4` / `bw4_overlap`.

## Decisions

1. Overlap block: frozen-only `G_F` (new temporary `PP_F`), not `G` through `PP_snk`.
2. Both paths get the option. Non-overlap uses `M = T^dag` with `PP_sep`.
3. `num_intermediate_projectors` = total singular vectors per projector slice.
4. Basis computed in-line once per run; singular values logged to XML. Not written to disk.
5. Randomized SVD with Gram-matrix orthonormalization. Dense steps use `SB::svd`.
6. (2026-10-22) Default `svd_power_iterations` = 2, from the toy table above.

## Open

- XML tag name: the request was `Projector_Type`. Existing `NamedObject` tags are
  lower snake case, so this uses `projector_type`. Change it if the capitalized form is preferred.
- Read slice for the overlap block: q = c1+fw-1 (sink-only buffer, every path to the sink
  crosses it) vs c1+fw (first sink-active slice, more filtering, but `D_F` would then
  reach an active slice). The implementation uses c1+fw-1.
- Slab-solve tolerance: the basis probably does not need the production tolerance.
  It uses `invParam` as-is for now.
- Two-level use: storing the basis (with the frozen-field checksum) is a follow-up.
