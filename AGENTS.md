# PostTopOpt — Codex Instructions

## Project overview

PostTopOpt is a C++20 project for studying, implementing, and testing geometry extraction and post-processing algorithms in the context of topology optimization and CAE software.

Implemented: regular-grid DC, emilk-inspired corner-plane DC, and classic Marching Cubes. Smoothing and remeshing remain future work.

The project is a reusable C++ static library with separate example executables.

## Technology stack

- Language: C++20 (required by the current build and examples)
- Build system: CMake
- Linear algebra: Eigen
- Testing: CTest
- Platforms: Windows with GCC/MinGW and MSVC
- IDEs: Visual Studio Code and Visual Studio 2022/Visual Studio 2019
- Version control: Git

## Project structure

- `include/` — public C++ headers
- `src/` — library implementations
- `tests/` — numerical and integration tests
- `docs/` — algorithm explanations, mathematical notes, and architectural decisions

- examples/ — analytic fields and method-selection CLI
- results/ — exported OBJ geometry, tracked separately from build artifacts
- third_party/ — MC table attribution and license

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

Read `docs/decisions.md` and `docs/algorithm.md` when they are relevant and available.
## Current implementation and documentation

Read docs/mathematics.md for QEF modes, docs/method_comparison.md for measurements
and defects, and docs/experiment_plan.md for proposals that are not yet implemented.
Both DC variants share SVD and four combinations of regularization/cell constraints.
Defaults: regularization=0.001, regularizationEnabled=true, constrainToCell=true, bisection=true. The emilk-inspired variant
changes planes, not the solver. MC uses no QEF. Constrained QEF and MC already exist.
Separate visual observations from measurements and proposals from accepted dependencies.
Preserve attribution and license notices. Examples default to project-root results/;
build caches remain excluded from Git.


## Дополнение: CGAL, 10 ноября 2026

Добавлен необязательный метод cgal-dc. Собственные DC, emilk и MC сохранены.
CGAL имеет отдельный решатель QEF и соглашение о нуле; описание общего
решателя собственных DC к нему не относится. [Сборка, различия и проверки](docs/cgal.md).

## Структура после объединения CLI, 13 ноября 2026

Основная программа — apps/main.cpp (PostTopOptRun); аналитический банк — apps/AnalyticModels.cpp; CSV/PTOV ввод — apps/VolumeIO.cpp. Прежние программы сохранены в apps/legacy. PLY находятся в data/models, массивы — data/volumes, результаты — в отдельном подкаталоге results. Presets создают build/release и build/debug. build/legacy — архив, его кэши после переноса не использовать. Актуальные команды: docs/running.md. Генерируемые результаты не добавлять в Git автоматически. Не восстанавливать удалённые пользователем результаты. Сборки и тесты после реорганизации выполняет пользователь.
