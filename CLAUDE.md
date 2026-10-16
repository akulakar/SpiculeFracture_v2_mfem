# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MPI-parallel C++17 finite element code (MFEM + hypre) for regularized variational (AT2 phase-field) fracture of
spicules: a notched cylinder driven by a prescribed boundary displacement plus a volumetric swelling eigenstrain.
This directory is a modular rewrite of the older single-file driver in
`/Users/Aashutosh/Academic/Software/MFEM/projects/SpiculeFracture`; the FE spaces and integrators are meant to match it.
There is no git repo, no test suite (`tests/` is empty) and `README.md` is empty.

## Build and run

```bash
cmake -S . -B build && cmake --build build          # Release by default
cd build && mpirun -np 4 ./SpiculeFracture ../input/inputParameters/SpiculeFracture.json [override.json]
```

- Run from `build/`: the mesh path, `output.path` default (`../results`) and the JSON paths are relative to the CWD.
- `argv[2]` is an optional override JSON merged onto the base with RFC 7386 merge-patch (`json::merge_patch`),
  handy for short smoke runs, e.g. `{"Test_Name":"smoke","Boundary_Conditions":{"nonzero_displacement":{"final_displacement":2e-5}}}`.
- `CMakeLists.txt` globs `source/*.cpp` into static lib `spicule` and `drivers/*.cpp` into one executable per file.
  Dependency paths are hardcoded Homebrew/`~/Academic/Software` locations (`MFEMPLUS_DIR`, METIS, nlohmann-json,
  ADIOS2); override with `-D` cache variables.
- mfemplus (`~/Academic/Software/MFEM_v2/mfemplus`) is pulled in via `add_subdirectory` and linked as
  `mfemplus::fracture`; headers are included as `<mfemplus/SolidMechanics/mfemplus_solids.hpp>`.

## Architecture (include/spicule + source, namespace `spicule`)

`drivers/SpiculeFracture.cpp` wires the pipeline: Config -> serial mesh -> ParMesh (refine on the ParMesh) ->
`FESpaces` -> `Fields` -> `Forms` (+ `forms.AddIntegrators(config)`) -> `BoundaryConditions` -> `PhaseFieldSolver` ->
`OutputWriter`.

- **Config**: JSON -> typed structs (`Get` throws on missing key, `GetOr` has a fallback). Keys are
  case-sensitive and nested (`Section.sub_section.key`); a name mismatch in a `GetOr` silently falls back to the
  default, so check key spelling against `input/inputParameters/*.json` when a setting seems ignored.
  `fe.ordering` is matched exactly (`byNODES` / `byVDIM`); anything else silently becomes `byNODES`.
- **Setup**: `FESpaces` (one H1 `FECollection` shared by a vector displacement space and a scalar damage space, both
  `byNODES`; makes the mesh nodal in the displacement space), `Fields`, `BoundaryConditions`.
- **Fields**: `u, d` are the accepted state; `utemp, dtemp` are the working copies the integrators actually read.
  `Update()` copies temp -> accepted, `Revert()` the reverse (used when a load step is rejected).
- **Forms**: mfemplus integrators, added in `Forms::AddIntegrators(config)` (the constructor only stores coefficients).
  `Kdisp(d)` damaged elasticity (reads `dtemp`), `Kdmg1` = mass(Gc/eps) + diffusion(Gc*eps) assembled once, `Kdmg2(u)` + `fdmg(u)` strain-energy driven (read `utemp`), `fdisp(d)` eigenstrain
  body force. mfemplus damage forms can be exactly zero at first assembly: assemble with `Assemble(0)`/`Finalize(0)`.
- **PhaseFieldSolver**: staggered (alternate minimisation) loop inside adaptive displacement-load stepping.
  Per iteration: assemble/solve displacement, assemble/solve damage, then irreversibility `d = max(d_new, d_old)`;
  convergence is the global max-norm of the damage update. Failed steps `Revert()` and shrink the increment.
  Preconditioner and Krylov objects are persistent (`std::unique_ptr<mfem::Solver>`); only the operator is swapped
  each iteration (the AMG hierarchy is rebuilt every iteration). Both solves warm-start (`iterative_mode = true`).
- **Linear solvers**: each of `Solvers.elasticity_solver` / `Solvers.damage_solver` independently picks
  `library` (`mfem` | `hypre`) and `type` (`pcg` | `gmres`) -> `CGSolver` / `GMRESSolver` / `HyprePCG` / `HypreGMRES`,
  built by `MakeKrylov` in `PhaseFieldSolver.cpp`. `preconditioner` is `amg_systems` (`SetSystemsOptions`, byNODES),
  `amg_elasticity` (`SetElasticityOptions`, needs `byVDIM`), `amg` (HMIS, aggressive coarsening 2, relax 6, V-cycle,
  strong threshold R 0.25) or `none`.
- **Output**: ADIOS2 BP5 via `mfem::adios2stream` (or ParaView), plus the resolved input JSON and a details file in
  `<path>/<testName>/`. `SaveStep` callback in the driver decides which steps are written.

## Comparing against the old driver

The older single-file driver builds as `Debug` (`-g -O0`) while this project defaults to Release, yet the old one was
observed faster in some runs. Compare per-stagger-iteration time and total stagger iterations/load steps, and make sure
solver library, tolerance, `max_iterations` (old: 300) and load parameters match before concluding anything.

## Conventions and gotchas

- BCs come from `Boundary_Conditions.dirichlet_attributes` (attributes, components, value | `"nonzero_displacement"`);
  the string value means "follows the load controller", whose range is in `Boundary_Conditions.nonzero_displacement`.
  `BCEntry::attrs` are 0-based indices into the mesh's boundary attributes (input `0` = mesh attribute 1), and
  `BoundaryConditions` marks `marker[a]` directly; MFEM attributes themselves are 1-based.
- **`tolerance` means different things per library**: mfem solvers (`CGSolver`, `GMRESSolver`) reduce the residual
  relative to the *initial* preconditioned residual; hypre solvers are relative to `||b||`. With warm starts an mfem
  solver at the same tolerance is far tighter (and slower) than a hypre one. Hypre solvers expose no convergence flag,
  so `Converged()` treats reaching `max_iterations` as failure.
- `SetSystemsOptions(dim, true)` in the displacement AMG assumes `byNODES` ordering; keep them consistent if the
  ordering becomes configurable.
- The eigenstrain is currently a constant `beta * deltaC`; the old driver used `beta * DeltaC * 0.01` (and a commented-out
  `ln(r/a)` radial variant), neither of which is ported. `SuppressDamage` and stress/strain output are also not ported.
- Style: 4-space indent, contents of `namespace spicule` indented, `mfem::` qualified in headers (no
  `using namespace` in headers); C++ is auto-formatted on save by the VS Code C/C++ extension.
