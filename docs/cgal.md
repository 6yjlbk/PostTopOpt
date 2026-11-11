# CGAL Dual Contouring

Добавлено 10 ноября 2026 года: необязательный четвёртый метод `cgal-dc`.
Адаптер: include/CgalExtraction.h и src/CgalExtraction.cpp. Остальные методы сохранены.

## Сборка и запуск

На этой машине установлена CGAL 6.2 из MSYS2 UCRT64 (пакет mingw-w64-ucrt-x86_64-cgal).
CMake требует CGAL >= 6.1. Для текущего MinGW:

```powershell
cmake -S . -B build -DPOSTTOPOPT_WITH_CGAL=ON -DCMAKE_PREFIX_PATH=C:/msys64/ucrt64 -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
./build/PostTopOptSphereUnion.exe --method cgal-dc
./build/PostTopOptDrilledBox.exe --method cgal-dc
./build/PostTopOptGyroid.exe --method cgal-dc
```

По умолчанию опция выключена. `-DPOSTTOPOPT_WITH_CGAL=OFF` позволяет собирать
проект без CGAL; запрос cgal-dc тогда выдаёт понятную ошибку.
Бинарные зависимости UCRT64 нельзя смешивать с MSVC.
Все пять примеров поддерживают ключ; OBJ по умолчанию записываются в results/.

## Отличия от собственного DC

Передаются те же поля, AABB и число ячеек; CGAL получает число узлов = ячейки + 1.
Пересечения линейные, нормали — аналитический градиент либо центральные разности.
Запуск последовательный, выходные грани — quads без дополнительного исправления.
CGAL 6.2 решает QEF через SVD нормальной матрицы 3x3 с порогом 1e-3 и смещением
относительно среднего пересечений. Ограничение ячейкой обрезает координаты;
наши DC используют SVD исходной матрицы и перебор 27 состояний ограниченной QEF.
Наша regularization в CGAL не передаётся. Отдельный CgalContourOptions содержит
isovalue, constrainToCell и gradientStep.

CGAL считает внутри f <= isovalue, собственные методы — f < isovalue.
Это влияет на активные ячейки при попадании поверхности в узлы, особенно у блока.
Округление координат также влияет на значения около нуля. Поэтому этот эксперимент
сравнивает реализации целиком, а не изолированно только решатели QEF.

## Проверка

Сборки ON и OFF проверены. ON: 6 выбранных тестов прошли, включая три фигуры;
OFF: тест адаптера и существующий тест методов прошли (2/2).
Тест адаптера проверяет анизотропную сетку, ненулевой изоуровень, ориентацию
плоскости, аналитическую нормаль, пустой результат и ошибочные входы.

| Файл в results/ (CTest) | Вершины | Quads |
|---|---:|---:|
| SphereUnion_cgal-dc.obj | 7462 | 7462 |
| DrilledBox_cgal-dc.obj | 17878 | 17886 |
| Gyroid_cgal-dc.obj | 63952 | 63990 |

CLI использует имена sphere_union_cgal-dc.obj, drilled_box_cgal-dc.obj,
gyroid_cgal-dc.obj. Сфера и тор в этом прогоне не рассчитывались.
Успех тестов не доказывает отсутствие самопересечений или превосходство по точности.

## Источники

[Руководство CGAL](https://doc.cgal.org/latest/Isosurfacing_3/index.html).
Детали решателя сверены с установленным internal/dual_contouring_functors.h CGAL 6.2.
Заголовки пакета указывают GPL-3.0-or-later OR LicenseRef-Commercial;
условия распространения необходимо учитывать при включении в продукт.

## DC с octree (11 ноября 2026)

Добавлен экспериментальный метод `--method cgal-dc-octree` в общий запуск примеров.
Обычный `cgal-dc` сохраняет регулярную сетку. Сборка и расчёты нового варианта
не выполнялись по просьбе пользователя.

Настройки в `include/CgalExtraction.h`: `octreeMinDepth=3`, `octreeMaxDepth=7`.
До минимальной глубины область делится полностью, затем до максимальной
уточняются ячейки с разными знаками в углах (равенство считается внутри).
На глубине d размер по каждой оси равен соответствующей стороне AABB / 2^d.
Grid задаёт AABB через origin + spacing*cells; разрешением octree управляют
глубины, а не cells напрямую. Это уточнение около поверхности, не по ошибке QEF
и не специально около вершины конуса. При minDepth=maxDepth дерево равномерное.
Параметры биссекции, градиента и ограничения вершины используются общие с CGAL DC.

Пример запуска после сборки:
```powershell
.\build\PostTopOptConeNumeric.exe --method cgal-dc-octree
.\build\PostTopOptConeAnalytic.exe --method cgal-dc-octree
```
Результаты: `results/cone_numeric_cgal-dc-octree.obj` и
`results/cone_analytic_cgal-dc-octree.obj`.

Алгоритм уточнения следует идее [официального примера CGAL](https://doc.cgal.org/latest/Isosurfacing_3/Isosurfacing_3_2contouring_octree_8cpp-example.html).
Проверка только углов может пропустить поверхность внутри ячейки; разная глубина
около поверхности может давать трещины. Это не реализация гарантированно
согласованного адаптивного DC и не гарантия устранения дефекта конуса.
