# Child-Local Factorized Baryon Tensors Spec (Draft)

Status: draft for discussion. Items marked **[DECIDED]** were settled in review
(see "Resolved Decisions"). Items marked **[PROVISIONAL]** are first-version
choices that can still change. The remaining questions are under "Open
Questions".

## Resolved Decisions

1. **Projector placement.** `phi_m` lives on frozen slices adjacent to the
   source child's active region. `D01` moves from `L~` into `R~`, and `D00^-1`
   stays in `L~`. (See "Locality Audit".)
2. **Spin content.** In the first version all six spins stay open, and a size
   guard (`max_output_gb`) protects against oversized outputs. Reducing spins
   for production is follow-on work.
3. **Solve strategy.** Both sides are solved from `phi_m` sources only:
   `M * Ns = 8 * Nproj` right-hand sides per solve stage. LEFT has two stages
   (`D00^-1`, then `D11^-1`) and RIGHT has one (`D11^-1`). One solve chain
   covers every sink time (LEFT) or every source time (RIGHT) in the child's
   active region. There is
   no per-time-source inversion and no `t_sources` input.
4. **Projector eigenvectors.** Both children read one shared frozen-slice
   colorvec file. It is generated once and keyed by parent-global time, so
   both sides use identical `phi_m` (same phases, same ordering).
5. **Time labels.** Every time given in the XML (baryon `t_slices`) is
   parent-global. The measurement converts to child-local time through the
   sidecar's `local_to_global_t`. The output is parent-global too.
6. **Fermion action scope.** The first version targets clover with a stout
   fermion state and **no smearing in the time direction**
   (`orthog_dir = decay_dir`). With no temporal smearing, the smeared links used
   by `D01` and `D10` never reach across a child edge, so the factorization at
   the child boundaries is preserved. The cut-link `D00_sep`/`D11_sep`
   approximation carries over from the prototype.
7. **Temporal BC mapping.** The parent temporal phase goes onto the
   child-internal link that crosses parent `Lt-1 -> 0` (see "Uncut field").
8. **Prototype refactor.** The prototype's helpers move into a shared unit that
   both measurements use (see "Code Reuse Plan").

## Goal

Add one inline measurement that runs on a single child sublattice and computes
one of the two bracketed tensors in Eq. (18) of the factorized-perambulator
notes (B. Slimmer, "Notes on Peardon et al. & factorized perambulators",
July 2026):

```text
C~_B(t',t) = [ Phi^{ijk} L~_d^{i,l} L~_u^{j,m} L~_u^{k,n} ](t') [ R~_d^{l,i'} R~_u^{m,j'} R~_u^{n,k'} Phi^{i'j'k'*} ](t)
           - [ Phi^{ijk} L~_d^{i,l'} L~_u^{j,m'} L~_u^{k,n'} ](t') [ R~_d^{l',i'} R~_u^{m',k'} R~_u^{n',j'} Phi^{i'j'k'*} ](t)
```

- the **LEFT** tensor is the sink bracket at `t'`, built on the sink child
- the **RIGHT** tensor is the source bracket at `t`, built on the source child

The measurement is given one child configuration plus the split sidecar, the
usual baryon-elemental input, the usual distillation propagator input, the
number of intermediate projectors, and a `side` selector (`LEFT` or `RIGHT`).
It writes the selected tensor for that child. Contracting a LEFT tensor from
one child with a RIGHT tensor from the other child gives the factorized
baryon two-point function; that contraction step is not part of this
measurement (see "Downstream Contraction").

The implementation should lean heavily on the parent-lattice prototype in
`lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc`
(the `factorized_geometry.enabled` branch), which already carries out the
full `D11^-1 D10 D00^-1 D01 D11^-1` chain with intermediate projectors on one
parent lattice. This measurement splits that chain across the two children.

## Relationship To Existing Specs And Code

- Split geometry, sidecar format, child-local time maps:
  `specs/hier/gauge_subdomain_split.md`,
  `lib/util/gauge/gauge_subdomain_split.{h,cc}`
- Exact (projector-free) factorized propagator outline:
  `specs/hier/factorized_prop.md`
- Parent-lattice prototype with intermediate projectors:
  `lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc`
- Baryon elementals (input schema and the colour contraction to reuse):
  `lib/meas/inline/hadron/inline_baryon_matelem_colorvec_superb_w.cc`
  (`SB::doMomDisp_colorContractions`)

## Notation

Follow the notes: `Lambda_0` is the frozen boundary region, `Lambda_1` is the
active (unfrozen) region. With the field ordered `(Lambda_0, Lambda_1)`:

```text
D = [ D00  D01 ]      D01 : active -> frozen   (row frozen, column active)
    [ D10  D11 ]      D10 : frozen -> active   (row active, column frozen)
```

The factorized perambulator (notes Eq. 12, with `S00^-1 ~ D00^-1`) is:

```text
tau_fact(t',t) = V^dag(t') D11^-1 D10 D00^-1 D01 D11^-1 V(t)
```

The notes write `V(t')` on the sink; this spec uses `V^dag(t')`.

Split geometry from `gauge_subdomain_split.md`:

```text
child0 local ordering : F0 + A + F1      (frozen widths fw, active length lenA)
child1 local ordering : F1 + B + F0      (frozen widths fw, active length lenB)
child-local frozen intervals: [0, fw-1] and [Lc-fw, Lc-1]
```

`F0` starts at parent time `cut_left`, `F1` at `cut_right`. Both frozen blocks
appear in both children with identical link data.

## Locality Audit (Why Eq. 17 Needs One Adjustment)

A child run only has its own active-region links plus the shared frozen
links. Checking each piece of the chain for Wilson-type hopping
`-1/2 [(1-g_mu) U_mu(x) psi(x+mu) + (1+g_mu) U_mu^dag(x-mu) psi(x-mu)]`:

| Piece | Links it needs | Available on |
|---|---|---|
| `D11^-1` (source side) | source active links | source child |
| `D01` (source active -> frozen) | `U_t(A[last])` (rooted in source active), `U_t(F0[fw-1])` (frozen) | source child only |
| `D00^-1` (cut-link version, see below) | frozen links only | either child |
| `D10` (frozen -> sink active) | `U_t(F1[fw-1])` (frozen), `U_t(B[last])` (rooted in sink active) | sink child only |
| `D11^-1` (sink side) | sink active links | sink child |
| Laplacian eigenvectors on a slice | spatial links on that slice | the child owning the slice; both if frozen |

Consequences:

1. `D01` depends on a link rooted in the source child's active region, so it
   cannot sit in the LEFT (sink) tensor. Notes Eq. (13)/(17) place `D01` in
   `L`; on children it has to move to `R`.
2. The projector vectors `phi_m` have to be the same objects in both children,
   or `sum_m L~^{i,m} R~^{m,j}` does not insert a projector. The only slices
   whose eigenvectors both children can build identically are frozen slices.
   Eq. (17) as written puts `phi_m` on source-active slices (they appear as
   `D01 phi_m`); the prototype puts them on sink-active slices (after `D10`).
   Neither works child-locally.

**[DECIDED]** child-local factorization:

```text
tau~(t',t) = sum_m L~^{i,m}(t') R~^{m,j}(t)

L~^{i,m}(t') = V^dag(t') D11^-1 D10 D00^-1 phi_m        (sink child)   LEFT
R~^{m,j}(t)  = phi_m^dag D01 D11^-1 V(t)                 (source child) RIGHT
```

with `phi_m` supported on the frozen slices adjacent to the source child's
active region. This keeps the notes' grouping of the frozen inverse `D00^-1`
with `L` and moves only `D01` to `R`.

A side benefit is that `D01 D11^-1 V` on a frozen slice is a pure temporal hop,
so its spin content is rank 2 (`(1 +/- gamma_t)` projected). That allows a later
factor-of-2 reduction of the projector spin index (see "Later Optimizations").

## Projector Index

For each frozen boundary `b in {0, 1}` (`b = 0` means `F0`, `b = 1` means `F1`;
the labels are parent-global, not child-local block order), use the
first `Nproj = num_intermediate_projectors` Laplacian eigenvectors
`phi_{b,v}`, `v < Nproj`, on one projector slice `g_b`:

| Source child | `g_0` (in `F0`) | `g_1` (in `F1`) |
|---|---|---|
| child0 (active `A`) | `cut_left + fw - 1` | `cut_right` |
| child1 (active `B`) | `cut_left` | `cut_right + fw - 1` |

(parent-global times, mod `Lt`). For `fw = 1` both rows coincide.

In child-local coordinates, the projector slices are
- on the source child (RIGHT run): local `fw-1` and `Lc-fw`
- on the sink child (LEFT run): local `0` and `Lc-1`

The combined projector index is `m = (b, v)` with size `M = 2 * Nproj`, stored
boundary-major (`m = b * Nproj + v`). Spin is carried separately, as in the
prototype (projection is colour x space only, `phi phi^dag` acting per spin).

**[PROVISIONAL]** A LEFT run on child `c` assumes the source child is `1 - c`.

### Projector Eigenvector File

**[DECIDED]** `phi_{b,v}` come from one colorvec file computed **once on the
parent lattice** with the ordinary `CREATE_COLORVECS_SUPERB`. Both the LEFT and
RIGHT runs pass that same file as `NamedObject/projector_colorvec_files`, and no
conversion step is needed.

The only subtlety is how the file is read. The parent file is indexed by parent
time `0 .. Lt-1`, while a child run's layout has child time `0 .. Lc-1`. The
stock `SB::getColorvecs` indexes the storage tensor with the current layout's
time (`s3t.kvslice_from_size({{'t', t_slice}, ...})`, where `t_slice` is a
child-local slice). Pointed at a parent file, it would silently read the wrong
parent slice. The measurement therefore reads projector slices through a small
helper:

```text
for each boundary b:
  g_b   = projector slice (parent-global)
  t_loc = the child-local slice with local_to_global_t[t_loc] == g_b
  read  parent-file slice t = g_b, vectors n < Nproj
  place into the child-local lattice tensor at t = t_loc
```

This is the same per-slice read, natural-to-red-black reorder and copy that
`getColorvecs` performs; only the file index changes. Only the spatial extents
have to match between parent and child, and they always do.

Requirements on the helper and the file:

- Check that the file's `lattSize` equals the sidecar's `parent_nrow`.
- Read raw stored vectors only. Skip the `fingerprint` recomputation branch
  of `getColorvecs`, which would recompute from the child gauge field with
  child-local indexing.
- A null (missing) vector on `g_b` is an error.
- Record the file name and a checksum of the `g_b` vectors in the output
  metadata, so the downstream contraction can refuse LEFT/RIGHT pairs built
  from different projector files.

Because frozen spatial links never change during child evolution, the parent
colorvecs on the frozen slices stay valid for every child update from the same
parent sample. Any parent-sized colorvec file from that sample works, including
the per-stitched-config `eigs_N.sdb` files already in use, since they agree on
the frozen slices up to eigensolver phase conventions. For strict
reproducibility, the recommended practice is still one file per parent sample.
The colorvec link smearing must be spatial only (`no_smear_dir = 3`), which is
already the convention.

`colorvec_files` remain child-local: they hold the child's own eigenvectors,
used for `V(t)` and `Phi`.

## Child-Local Operators

These follow the prototype's construction exactly, applied to the child layout
instead of the parent.

### Cut-link field `u_sep`

Starting from the child gauge field, zero the `t_dir` links rooted on:

- local `fw - 1` (`F_first | active` interface)
- local `Lc - fw - 1` (`active | F_second` interface)
- local `Lc - 1` (the child's periodic wrap, which is not a physical link of the
  parent)

This is the prototype's loop over `frozen_local_intervals` calling
`zeroTemporalLinksOnSlice(t_start - 1)` and `zeroTemporalLinksOnSlice(t_end)`
with `Lt` replaced by the child extent `Lc`. With `u_sep`, the Dirac operator is
block diagonal:

```text
D_sep = diag(D00_sep[F_first], D11_sep[active], D00_sep[F_second])
```

and one `SB::ChimeraSolver` on `u_sep` provides both `D00_sep^-1` and
`D11_sep^-1` (as in the prototype). With fermion-state smearing and clover
built from `u_sep`, `D00_sep` depends only on frozen links and is identical on
both children.

Note: for clover and stout fermion states, `D00_sep` and `D11_sep` are not
exactly the blocks `D_{Lambda00}` and `D_{Lambda11}` of the full operator;
the clover term and smeared links near the cut differ. This is the same
approximation the prototype already makes, and the spec keeps it.

### Uncut field `u_bc` for the coupling blocks

`D01` and `D10` are applied as in the prototype: apply the full child `LinOp`
built from the uncut field, then restrict with `restrictToTimeslices`. The
child's own periodic-wrap hop only connects frozen slices to frozen slices,
and the restriction discards it.

**[DECIDED]** Temporal fermion BC mapping. The parent's temporal BC phase
(for example `-1` from `SIMPLE_FERMBC boundary 1 1 1 -1`) sits on the parent
link `Lt-1 -> 0`. On a child, that link may be internal (for example
`B[last] -> F0[0]` in child1 when `cut_left = 0`), where the child's own FermBC
does not reach it. The measurement multiplies the `t_dir` links rooted on the
child-local slice `t*` with `local_to_global_t[t*] = Lt - 1` and
`local_to_global_t[t* + 1] = 0` by the parent temporal phase, in both `u_bc`
and `u_sep`. The FermionAction XML is taken as the parent action; its phase on
the child wrap link is harmless because that link is cut in `u_sep` and
discarded in `u_bc` after restriction. The prototype does not need this
because it runs on the parent.

**[DECIDED]** No temporal smearing. The fermion state must not smear in the time
direction. For `STOUT_FERM_STATE`, that means `orthog_dir == decay_dir`: the
`t_dir` links are left unsmeared, and spatial links get only spatial staples
(`smear_in_this_dirP[t_dir] = false` in `stout_utils.cc`). An unsmeared
fermion state is also allowed. Every link entering `D01`, `D10`, `D00_sep` and
`D11_sep` then depends only on links rooted on its own time slice, or is a
thin temporal link. So no piece reaches across a child edge or the child wrap,
for any `fw >= 1` and any `n_smear`. The measurement parses the FermState and
rejects any temporal smearing (for example `orthog_dir = -1`) with a clear
error. Note: the current real-lattice test XMLs use `orthog_dir = -1` and must
be changed to `3`.

## Algorithms

**[DECIDED]** Both sides are driven only by `phi_m` sources. Each run does
`M * Ns = 2 * Nproj * Ns` solves per solve stage, and one chain covers every
time slice of the child's active region. The prototype's step of looping
over `t_sources` and inverting `V(t)` is not used.

All solves go through `SB::doInversion(PP, ...)` with the
`MultipleLatticeFermions` overloads and `max_rhs` batching, exactly as in the
prototype. Throughout, `R_act` restricts to the child's active slices and
`R_proj` restricts to the projector slices.

### LEFT (sink child)

For every `(b, v, mu)` (`M * Ns` sources, batched):

1. `eta = phi_{b,v} (x) e_mu` on slice `g_b` (sink-child local `0` or `Lc-1`)
2. `chi = D00_sep^-1 eta` (solve on `u_sep`; prototype `contract2`)
3. `y = R_act( D_bc chi )` (prototype `y_boundary2`, `y_boundary2_rs`)
4. `zeta = D11_sep^-1 y` (solve on `u_sep`; prototype `contract3`)
5. `L~^{(I,P),(b,v,mu)}(t') = v_I(t')^dag zeta_P(t')` for every active `t'`
   (prototype sink-colorvec `elems.contract`)

A single chain (two solve stages of `M * Ns` right-hand sides each) covers
every sink time `t'` in the child's active region. The number of solves does
not depend on the number of sink times.

### RIGHT (source child)

**[DECIDED]** Adjoint, `phi_m`-sourced route. Because `D^dag = g5 D g5`
holds blockwise, `(D01)^dag = g5 D10 g5` and
`(D11_sep^-1)^dag = g5 D11_sep^-1 g5`. Therefore

```text
R~^{m,j}(t)^* = v_j^dag(t) g5 D11_sep^-1 D10 g5 phi_m
```

For every `(b, v, mu)`:

1. `eta = g5 (phi_{b,v} (x) e_mu)` on slice `g_b` (source-child local `fw-1`
   or `Lc-fw`)
2. `y = R_act( D_bc eta )`
3. `zeta = g5 D11_sep^-1 y`
4. `R~^{(b,v,mu),(J,beta)}(t) = conj( v_J(t)^dag zeta_beta(t) )` for every active `t`

This is the LEFT chain with `D00^-1` removed and `g5` sandwiches added. It
costs `M * Ns` solves and covers all source times at once. The
prototype-style direct route
(`phi^dag R_proj D_bc D11_sep^-1 V(t)`, prototype `contract1` and
`y_boundary1`) is not part of the measurement. It survives only as an
independent cross-check in the focused test (validation item 4).

### Baryon Elementals And Tensor Assembly

Compute `Phi^{IJK}(t)` on the child for the requested times, displacements,
momenta and phasings with the same `SB::doMomDisp_colorContractions` call and
`ColorContractionFn` callback pattern as `BARYON_MATELEM_COLORVEC_SUPERB`.
The child's colorvecs and smeared links are used. Elementals are timeslice
local, so the result matches the parent elemental on that child's active
slices.

Inside the callback, contract per `(t, d, mom, h)` chunk:

```text
LEFT_{ijk,pqr,PQR}(t')  = sum_{IJK} Phi^{IJK}(t') L~^{(I,P),(i,p)} L~^{(J,Q),(j,q)} L~^{(K,R),(k,r)}
RIGHT_{ijk,pqr,PQR}(t)  = sum_{IJK} R~^{(i,p),(I,P)} R~^{(j,q),(J,Q)} R~^{(k,r),(K,R)} Phi^{IJK}(t)^*
```

Contract one slot at a time (`K` first, then `J`, then `I`) with `SB::contract`,
chunked by `max_tslices_in_contraction` and `max_moms_in_contraction`.

**[PROVISIONAL]** The quarks are mass-degenerate (single `mass_label`), so
`L~_u = L~_d` and `R~_u = R~_d`. Then one tensor per side serves both Wick
terms in Eq. (18): the second term is a permutation of the first (see
"Downstream Contraction").

Default time support:
- LEFT: all active slices of the sink child
- RIGHT: all active slices of the source child
- An optional baryon `t_slices` list (parent-global) narrows which slices are
  contracted and written; it does not change the solves. Entries that do not
  map onto the child's active region are rejected.

## Spin Content And Output Size

**[DECIDED]** All spins stay open in the first version. Written in full, each quark slot carries a
projector-side spin (`p, q, r`) and a baryon-operator-side spin (`P, Q, R`).
Per `(t, d, mom, h)`, a tensor then holds `M^3 * Ns^6 = 4096 M^3` complex
doubles:

| `Nproj` | `M = 2 Nproj` | All spins open | Projector spin reduced to 2 |
|---|---|---|---|
| 4 | 8 | 34 MB | 4 MB |
| 8 | 16 | 268 MB | 34 MB |
| 16 | 32 | 2.1 GB | 268 MB |
| 64 | 128 | 137 GB | 17 GB |

(These sizes are for one time slice, one displacement, one momentum.) With
fully open spins, validation on small `Nproj` is feasible; production values
like the prototype's `Nproj = 64` are not. The first version therefore keeps
fully open spins with a size guard: the run aborts before any solve if the
estimate exceeds `max_output_gb`. Contracting the baryon-side spins with a
user-supplied spin structure, and reducing the projector spin, are the paths
to production.

## Output Format

One SUPERB `SB::StorageTensor` per run (`use_superb_format = true` required in
the first version). **[PROVISIONAL]** labels and order:

```text
order = "ijkpqrPQRtdmh"
  i, j, k : projector index m = b*Nproj + v for quark slots 1, 2, 3   (size M)
  p, q, r : projector-side spin for slots 1, 2, 3                      (Ns)
  P, Q, R : baryon-operator-side spin for slots 1, 2, 3                (Ns)
            (sink spin for LEFT, source spin for RIGHT)
  t       : parent-global time slice                                   (Lt_parent)
  d, m, h : displacement, momentum, phasing (as in the baryon elemental)
```

Slots 1, 2 and 3 follow the elemental's `left, middle, right` displacement
slots. Storage is `SB::Sparse`, indexed by parent-global time so that LEFT and
RIGHT files from different children share one time axis.

Required metadata (`DBMetaData`):

- `id = factorizedBaryonTensorSuperb`, `side`, `tensorOrder`
- `child_id`, `source_child_id`, `sink_child_id`
- the full `GaugeSubdomainSplitPlan` (sidecar contents) and the sidecar path
- `local_to_global_t`, active local and global slices
- projector slices `g_0, g_1` (global), `num_intermediate_projectors`,
  and the projector colorvec file identity
- `num_vecs`, `mass_label`, `Params` (propagator and fermion action)
- `displacements_left_middle_right`, `moms`, `phasings` (as in the elemental)
- `Config_info` and `lattSize` (child)

**[PROVISIONAL]** An optional `factorized_intermediate_file` also writes
`L~` or `R~` themselves, labelled `"nspmt"`-style (distillation index,
distillation-side spin, projector index, projector spin, time). This supports
the validation checks below and is cheap compared with the tensors.

## XML Input

**[PROVISIONAL]** Measurement name `FACTORIZED_BARYON_TENSOR_DISTILLATION_SUPERB`,
files `lib/meas/inline/hadron/inline_factorized_baryon_tensor_distillation_superb_w.{h,cc}`,
namespace `InlineFactorizedBaryonTensorDistillationSuperbEnv`.

```xml
<elem>
  <Name>FACTORIZED_BARYON_TENSOR_DISTILLATION_SUPERB</Name>
  <Frequency>1</Frequency>
  <Param>
    <!-- Prototype Contract_t schema; t_sources, Nt_forward and Nt_backward are not used -->
    <Contractions>
      <mass_label>U-0.2050</mass_label>
      <num_vecs>64</num_vecs>
      <decay_dir>3</decay_dir>
      <max_rhs>4</max_rhs>
      <use_superb_format>true</use_superb_format>
      <output_file_is_local>false</output_file_is_local>
    </Contractions>

    <!-- Same schema as the prototype's Propagator (ChromaProp_t); parent action -->
    <Propagator> ... </Propagator>

    <!-- Same schema as BARYON_MATELEM_COLORVEC_SUPERB's Param; num_vecs and
         decay_dir are taken from Contractions and must match if given -->
    <BaryonElemental>
      <use_derivP>false</use_derivP>
      <max_tslices_in_contraction>4</max_tslices_in_contraction>
      <max_moms_in_contraction>1</max_moms_in_contraction>
      <max_vecs>8</max_vecs>
      <combos>
        <elem><phase>0 0 0</phase><mom_list><elem>0 0 0</elem></mom_list></elem>
      </combos>
      <displacement_list>
        <elem><left>0</left><middle>0</middle><right>0</right></elem>
      </displacement_list>
      <LinkSmearing>
        <LinkSmearingType>STOUT_SMEAR</LinkSmearingType>
        <link_smear_fact>0.08</link_smear_fact>
        <link_smear_num>10</link_smear_num>
        <no_smear_dir>3</no_smear_dir>
      </LinkSmearing>
    </BaryonElemental>

    <FactorizedTensor>
      <side>LEFT</side>                                  <!-- LEFT | RIGHT -->
      <num_intermediate_projectors>8</num_intermediate_projectors>
      <max_output_gb>64</max_output_gb>                  <!-- size guard -->
    </FactorizedTensor>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>             <!-- child config -->
    <colorvec_files><elem>child1_eigs.sdb</elem></colorvec_files>
    <projector_colorvec_files><elem>frozen_eigs.sdb</elem></projector_colorvec_files>
    <frozen_boundary_sidecar_file>parent.sidecar.xml</frozen_boundary_sidecar_file>
    <frozen_boundary_child_id>1</frozen_boundary_child_id>
    <factorized_tensor_file>left_child1.sdb</factorized_tensor_file>
    <factorized_intermediate_file>Ltilde_child1.sdb</factorized_intermediate_file>  <!-- optional -->
  </NamedObject>
</elem>
```

The `<nrow>` and `<Cfg>` of the enclosing chroma input are the child's.

Input validation, before any solve:

- sidecar present; `child_id` is 0 or 1; `Layout::lattSize()` equals
  `childN_nrow`; sidecar `t_dir == decay_dir == Nd - 1`
- `side` is `LEFT` or `RIGHT`; `num_intermediate_projectors > 0`
- baryon `t_slices`, if given, are parent-global and map through
  `local_to_global_t` onto this child's active slices
- `projector_colorvec_files` holds every projector slice `g_b` with at least
  `Nproj` vectors
- `use_superb_format == true`; `Nc == 3`, `Ns == 4`
- estimated output size does not exceed `max_output_gb`
- the fermion state has no temporal smearing (`orthog_dir == decay_dir`, or
  no smearing)

The prototype keeps `num_intermediate_projectors` under `NamedObject`. This
draft moves it to `Param/FactorizedTensor`, because it is a physics parameter
rather than a named object.

## Code Reuse Plan

**[DECIDED]** Move the prototype's reusable helpers into a shared unit,
for example `lib/util/ferm/factorized_distillation.{h,cc}`, and have both the
prototype and the new measurement include it:

- `FactorizedPropGeometry`, `loadFactorizedPropGeometry`,
  `buildActiveIntervals`, `appendIntervalSlices`,
  `writeFactorizedPropGeometry`. Generalize to use child-local geometry for
  either child; the prototype hardwires `child0` intervals on the parent.
- `zeroTemporalLinksOnSlice`, `restrictToTimeslices`, `toLatticeFermions`,
  `toSBTensor`, `returnNLatticeFermions`
- new: `buildCutLinkField(u, geometry)`, `applyParentTemporalPhase(u, geometry, phase)`,
  `makeProjectorSources(...)`, `applyCouplingBlock(linop, psi, keep_slices)`

Copy the XML readers and writers for `Contract_t`, `Phasing_t` and
`ChromaProp_t` from the prototype. Copy `Displacement_t`, the
`combos`/`mom_list`/`phases` handling and `normalizeDisplacements` from the
baryon elemental, or factor them out the same way.

Wiring, per `AGENTS.md`: register in
`lib/meas/inline/hadron/inline_hadron_aggregate.cc`; add sources to both
`lib/Makefile.am` and `lib/CMakeLists.txt`; guard with `BUILD_SB`.

## Downstream Contraction (Reference Only)

Not implemented by this measurement, but the tensor layout is chosen for it.
With `Gsnk_{PQR}` and `Gsrc_{P'Q'R'}` the baryon operators' spin tensors:

```text
C~(t',t) = sum Gsnk_{PQR} Gsrc*_{P'Q'R'} sum_{ijk,pqr}
             LEFT_{ijk,pqr,PQR}(t') [ RIGHT_{ijk,pqr,P'Q'R'}(t) - RIGHT_{ikj,prq,P'Q'R'}(t) ]
```

The exchange term swaps the projector index and projector spin of slots 2 and
3 in RIGHT, while the operator-side spins stay attached to the elemental slots.
This is Eq. (18) with `L~_u = L~_d`. Signs and slot conventions have to be
checked against the prototype and elemental path in the validation below.

## Validation

A focused checker `mainprogs/tests/t_factorized_baryon_tensor.cc` (wired into
both test build systems) reads one LEFT file and one RIGHT file, performs the
contraction above for a fixed operator (for example the nucleon
`(C g5)` structure with a positive-parity projector), and compares the result
against a reference. Smoke XML goes under `tests/factorized_baryon_tensor/`,
and run artifacts go under `cfgs/`.

1. Geometry: projector slices `g_b` agree between a LEFT run on child `c`
   and a RIGHT run on child `1-c` (compare metadata); the parent BC link
   `t*` is found correctly for `cut_left = 0` and `cut_left != 0` splits.
2. Projector completeness, checked at the perambulator level with the
   `factorized_intermediate_file` outputs, because full-rank baryon tensors are
   far too large. On a small lattice (for example `4^3 x 16`), set
   `Nproj = Nc * Ls^3` so that `Pi` is the identity on `g_b`. `D01 D11^-1 V`
   is supported only on the slices `g_b`, so `sum_m L~ R~` then has to equal
   the prototype's projector-free chain
   `V^dag D11^-1 D10 D00^-1 D01 D11^-1 V` (the commented `contract3_db` path)
   on the unevolved parent, to solver precision, for any `fw`.
3. Unevolved split: split one parent config, then run RIGHT on child0 and LEFT
   on child1 with no child HMC. The contracted `C~` has to match a baryon
   correlator built from the prototype-style parent perambulator (same
   projector placement) and parent elementals.
4. Adjoint cross-check. For one source time `t`, the focused test builds
   `R~(t)` directly as `phi^dag R_proj D_bc D11_sep^-1 V(t)`, the prototype's
   `contract1`/`y_boundary1` route. It then compares that result with the
   measurement's `phi`-sourced `R~(t)`. This tests the `g5` hermiticity of the
   blocks and the conjugation convention. The matching LEFT check compares one
   sink slice against a direct `V(t')`-sourced solve.
5. Wick exchange: the permutation form of the second term matches an explicit
   evaluation from `L~`/`R~` intermediates.
6. Real lattice (`32^3 x 64` clover): convergence of `C~` in `Nproj` toward
   the exact correlator, in the style of `bw_testing`.

## Later Optimizations (Not First Version)

- Projector-spin reduction: `D01 D11^-1 V` on `g_b` lies in the image of the
  temporal hop projector, so `p, q, r` can use a 2-dimensional spin basis
  (8x smaller tensors).
- Solve `D00_sep^-1` on the frozen block only rather than the whole child
  lattice.
- Contract baryon-side spins with supplied operators at write time (see Open
  Question 2).

## Out Of Scope

- the LEFT x RIGHT production contraction code and its redstar integration
  (beyond the checker)
- two-level outer and inner averaging estimators and bias corrections
- more than two children
- non-degenerate quark masses

## Open Questions

No questions currently block implementation. The remaining **[PROVISIONAL]**
items (the measurement name, the output dimension labels, the optional
intermediate file, and the assumption that a LEFT run on child `c` has source
child `1 - c`) can be settled during implementation review.
