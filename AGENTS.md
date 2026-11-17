# PostTopOpt — Codex Instructions

## Project overview

PostTopOpt is a reusable C++20 static library for extracting surfaces from regular-grid scalar fields and element volume fractions (VF). Implemented backends: own DC, emilk-inspired DC, classic MC, optional CGAL DC and ParaView AMR. Smoothing, remeshing and a topology optimizer are not implemented.

Build: CMake, Eigen; optional VTK/ParaView/CGAL. Keep GCC/MinGW and MSVC compatibility. Local release/debug presets use MSYS2 UCRT64.

## Structure and documentation

- include/src: public API and library implementations.
- apps/main.cpp: PostTopOptRun, shared input preparation, benchmark and export.
- apps/AnalyticModels.cpp: the primary analytic model bank in one file.
- apps/VolumeIO.cpp: VF CSV and legacy PTOV1 I/O.
- apps/legacy: retained direct analytic examples and specialized programs.
- data/models and data/volumes: PLY and saved VF arrays.
- results: one directory per run; generated outputs are ignored by Git.
- build/release and build/debug: active presets; build/legacy: moved caches, do not use.
- tests, docs, third_party: checks, documentation, attribution and licenses.

Read docs/architecture.md for data contracts, DC/QEF and backend differences; docs/running.md for current CLI and metric meanings; docs/benchmark.md for measured results and proposed work; docs/implementations.md for the dated external implementation review. Keep measured facts separate from proposals. Do not restore results deleted by the user.

Both own DC variants share a QEF solver; emilk changes planes. Defaults: regularization=0.001, regularizationEnabled=true, constrainToCell=true, bisection=true; all remain configurable. CGAL has a separate solver; do not map our regularization to it. MC uses no QEF. ParaView receives native element VF through its adapter. Benchmark is on by default; false means one extraction per method and no summary/runs CSV. Report times are seconds with three decimal places.

## Development rules

1. Prefer clean, readable, modern C++20.
2. Keep the implementation modular and reusable.
3. Separate numerical algorithms from mesh data structures and file I/O.
4. Use Eigen for linear algebra when appropriate.
5. Avoid introducing unnecessary dependencies.
6. Preserve compatibility with GCC and MSVC.
7. Do not modify unrelated code.
8. Prefer small, focused changes over large refactorings.
9. Add or update tests when changing algorithmic behavior.
10. Explain important mathematical and architectural decisions.

## Dual Contouring guidelines

The current extractors use a regular Cartesian grid.

The classic DC pipeline is:

1. Scalar field representation and isovalue selection.
2. Detection of sign-changing grid edges.
3. Hermite data collection: intersection points and normals.
4. QEF minimization for vertex placement.
5. Mesh connectivity generation.
6. Mesh validation and optional export.

Pay particular attention to numerical robustness, degenerate configurations, mesh topology, sharp features, and preservation of thin structural members.

Do not assume that a basic Dual Contouring implementation always produces manifold, watertight, or self-intersection-free meshes.

## Collaboration rules

- Before significant modifications, briefly explain the proposed approach.
- Do not make major architectural decisions without discussing alternatives.
- Do not commit or push changes unless explicitly requested.
- Do not silently introduce new third-party libraries.
- When a mathematical assumption is uncertain, state it explicitly.
- Keep project documentation consistent with implemented behavior.


## Workflow constraints

Do not run builds, applications or tests until the user authorizes them again. Read-only inspection and documentation editing are allowed. Preserve source model files and third-party attribution. Do not commit or push without an explicit request.
