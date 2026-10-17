# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MPI-parallel C++17 finite element code (MFEM + hypre) for regularized variational (AT2 phase-field) fracture of
spicules: a notched cylinder driven by a prescribed boundary displacement plus a volumetric swelling eigenstrain.
This directory is a flattened rewrite of the older single-file driver in
`/Users/Aashutosh/Academic/Software/MFEM/projects/SpiculeFracture`; the FE spaces and integrators match it.
There is no git repo, no test suite (`tests/` is empty) and `README.md` is empty.
`old_modules/` holds the previous modular version (Config/Setup/Forms/Output/Solver as header+cpp pairs and its driver).
It is reference only and is not built.

## Build and run

```bash
cmake -S . -B build && cmake --build build          # Release by default
cd build && mpirun -np 4 ./SpiculeFracture ../input/inputParameters/SpiculeFracture.json [override.json]
```

- Run from `build/`: the mesh path, `Output.path` default (`../results`) and the JSON paths are relative to the CWD.
- `argv[2]` is an optional override JSON merged onto the base with RFC 7386 merge-patch (`json::merge_patch`),
  handy for short smoke runs, e.g. `{"Test_Name":"smoke","Boundary_Conditions":{"nonzero_displacement":{"final_displacement":1e-2}}}`.
- `CMakeLists.txt` globs `drivers/*.cpp` into one executable per file; there is no project library. mfemplus
  (`~/Academic/Software/MFEM_v2/mfemplus`) is added with `add_subdirectory` and linked as `mfemplus::fracture`
  (headers: `<mfemplus/SolidMechanics/mfemplus_solids.hpp>`). Dependency paths are hardcoded Homebrew/`~/Academic/Software`
  locations (`MFEMPLUS_DIR`, METIS, nlohmann-json, ADIOS2); override with `-D` cache variables.
- The editor's clang diagnostics (e.g. "unknown type name constexpr") come from a language-server config without C++17
  flags, not from the real build; trust `cmake --build`.

## Architecture

`drivers/SpiculeFracture.cpp` reads the input as a plain `nlohmann::json` (no typed config), builds the ParMesh (serial
mesh in a scope, `SetCurvature(order,false)`, partition, refine on the ParMesh), then drives the solver through small
init steps, each reading only its own JSON section:

`InitializeFESpaces -> InitializeFields -> InitializeForms -> InitializeBCs -> InitializeSolvers`, then
`OutputWriter`, then `RunStaggeredSolve(callback)`.

All classes live in `include/spicule/PhaseFieldFracture.hpp` (header-only, definitions in-class, namespace `spicule`):

- **PhaseFieldSolver** owns everything: the H1 collection and two spaces (vector displacement, scalar damage, ordering from
  `Finite_Element.ordering`), grid functions `u, d` (accepted) and `utemp, dtemp` (working copies the integrators read),
  the mfemplus forms, the BCs and the solvers. Each `Initialize*` throws if an earlier step was skipped.
  - Forms: `Kdisp(d)` damaged elasticity, `Kdmg1` = mass(Gc/eps) + diffusion(Gc*eps) assembled once, `Kdmg2(u)` + `fdmg(u)`
    strain-energy driven, `fdisp(d)` eigenstrain body force (`volEig = beta * Delta_concentration`). mfemplus damage forms can
    be exactly zero at first assembly, so assemble with `Assemble(0)`/`Finalize(0)`.
  - Staggered loop: per iteration assemble/solve displacement, assemble/solve damage, irreversibility `d = max(d_new, d_old)`,
    convergence = global max-norm of the damage update. Failed steps `Revert()` and shrink the increment; accepted steps `Accept()`.
    Both solves warm-start (`iterative_mode = true`). AMG is rebuilt every iteration. Four timers (disp/damage assembly and solve)
    are printed at the end.
  - Convergence checks of the linear solves are currently disabled (`SolveDisplacement/SolveDamage` return `true`); the
    `Converged()` helper is kept for re-enabling.
- **BoundaryConditions**: from `Boundary_Conditions.dirichlet_attributes` (attributes, components, value | `"nonzero_displacement"`).
  Attributes are 0-based indices into the mesh's boundary attributes (input `0` = mesh attribute 1). The string value means
  "follows the load controller"; the load range is `Boundary_Conditions.nonzero_displacement`.
- **OutputWriter**: ADIOS2 BP5 via `mfem::adios2stream` (or ParaView), plus the resolved `input.json` and `SimulationDetails.txt`
  in `<Output.path>/<Test_Name>/`. The driver's callback applies `Output.save_frequency`.

## Linear solvers

Each of `Solvers.elasticity_solver` / `Solvers.damage_solver` independently picks `library` (`mfem` | `hypre`) and `type`
(`pcg` | `gmres`) -> `CGSolver` / `GMRESSolver` / `HyprePCG` / `HypreGMRES`. `preconditioner` is `amg_systems`
(`SetSystemsOptions`, flag follows the space ordering), `amg_elasticity` (`SetElasticityOptions`, requires `byVDIM`),
`amg` (HMIS, aggressive coarsening 2, relax 6, V-cycle, strong threshold R 0.25) or `none`.

- **`tolerance` means different things per library**: mfem solvers reduce the residual relative to the *initial* preconditioned
  residual; hypre solvers are relative to `||b||`. With warm starts an mfem solver at the same tolerance is far tighter (and
  slower) than a hypre one.
- Hypre solvers expose no convergence flag, so `Converged()` treats reaching `max_iterations` as failure.

## Conventions and gotchas

- Keys are case-sensitive and nested; `.value("key", default)` silently falls back on a misspelled optional key, `.at("key")`
  throws on a missing required one. Check spelling against `input/inputParameters/*.json` when a setting seems ignored.
- Unknown `Finite_Element.basis` / `ordering`, `Output.format`, solver `library`/`type`/`preconditioner` throw at startup.
- The eigenstrain is a constant `beta * Delta_concentration`; the old driver used `beta * DeltaC * 0.01` (and a commented-out
  `ln(r/a)` radial variant), neither is ported. `SuppressDamage` and stress/strain output are also not ported.
- The older single-file driver builds as `Debug` (`-g -O0`) while this project defaults to Release. When comparing timings,
  match solver library, tolerance, `max_iterations` (old: 300), load parameters, and compare stagger iterations per step too.
- Style: 4-space indent, contents of `namespace spicule` indented, `mfem::` qualified in headers (no `using namespace` in headers).
