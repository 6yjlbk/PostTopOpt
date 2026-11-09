# Геометрические примеры

Актуализировано 9 ноября 2026 года. Поля считают отрицательные значения
внутренней областью материала. Размеры фигуры и сетка извлечения независимы.

| Цель | Файл | Геометрия | AABB | Ячейки | Шаг |
| --- | --- | --- | --- | --- | --- |
| PostTopOptDemo | sphere.cpp | Сфера R=2 | [-3,5] по каждой оси | 80³ | 0.1 |
| PostTopOptTorus | torus.cpp | Тор R=1.05, r=0.38, ось z | [-1.8,1.8]³ | 60³ | 0.06 |
| PostTopOptSphereUnion | sphere_union.cpp | Объединение трёх сфер | [-1.8,1.8]³ | 60³ | 0.06 |
| PostTopOptDrilledBox | drilled_box.cpp | Куб с половиной стороны 1.15, каналы r=0.46 | [-1.5,1.5]³ | 60³ | 0.05 |
| PostTopOptGyroid | gyroid.cpp | Полосовой гироид внутри куба | [-1.6,1.6]³ | 80³ | 0.04 |

Сфера использует значения Grid по умолчанию из заголовка и аналитический градиент.
У остальных фигур Grid задан явно через examples::grid(), нормали DC вычисляются
центральными разностями. MC градиент не запрашивает.

Центры объединённых сфер: (-0.55,0,0), (0.55,0,0), (0,0.65,0.25);
радиусы 0.9, 0.9, 0.75. Объединение через min, без сглаживания.

Гироид:
$g=\sin(\pi x)\cos(\pi y)+\sin(\pi y)\cos(\pi z)+\sin(\pi z)\cos(\pi x)$.
Материал: $|g|<0.28$ внутри куба с половиной стороны 1.35.
Порог 0.28 не означает постоянную физическую толщину.
Куб входит в поле; это не автоматическое закрытие границ AABB алгоритмом.

## Запуск и выходные файлы

Из корня проекта, PowerShell:

~~~powershell
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
.\build\PostTopOptGyroid.exe --method dc
.\build\PostTopOptGyroid.exe --method emilk
.\build\PostTopOptGyroid.exe --method mc
~~~

Заменой имени EXE выбирается фигура. Без --method используется dc.
--help показывает синтаксис; неизвестный метод считается ошибкой.
Без пути создаются results/gyroid.obj, results/gyroid_emilk.obj,
results/gyroid_mc.obj в корне проекта независимо от текущей рабочей папки.
Каталог создаётся автоматически; повторный запуск перезаписывает тот же файл.

~~~powershell
.\build\PostTopOptSphereUnion.exe results/experiment_01/spheres_emilk.obj --method emilk
~~~

Явный путь имеет приоритет. Относительный путь отсчитывается от текущей папки.
Примеры создают родительские каталоги; библиотечный writeObj() этого не делает.
Путь по умолчанию задаёт CMake при сборке. На другом компьютере или после
перемещения проекта нужно выполнить конфигурацию и сборку заново.

Для Visual Studio EXE находятся в каталоге конфигурации, например
build-msvc/Release; CTest получает -C Release.

## Отдельная настройка AABB и шага

Вызов examples::grid(1.8,60) задаёт куб [-1.8,1.8]³ с 60 ячейками на ось:
шаг 3.6/60=0.06. Второй аргумент меняет разрешение при фиксированной области.

Для прямоугольной области можно явно создать Grid в файле фигуры:

~~~cpp
geometry::Grid grid;
const geometry::Vec3 minimum(-2.0, -1.5, -1.0);
const geometry::Vec3 maximum( 2.0,  1.5,  1.0);
grid.origin = minimum;
grid.cells = {80, 60, 40};
grid.spacing = (maximum - minimum).cwiseQuotient(
    geometry::Vec3(grid.cells[0], grid.cells[1], grid.cells[2]));
return examples::run(argc, argv, "torus.obj", grid, field);
~~~

Связь: maximum = origin + spacing ⊙ cells. Границы, шаг и число ячеек нельзя
назначать произвольно одновременно. Изменение AABB не масштабирует фигуру.
Оставлять запас до границы области. При уменьшении всех шагов вдвое число
объёмных ячеек растёт примерно в восемь раз.
После правки: cmake --build build --target PostTopOptTorus.
CLI для AABB/шага пока отсутствует.

## Настройки QEF и проверки

Четыре сочетания regularization и constrainToCell доступны через API обоих DC.
ExampleRunner использует настройки по умолчанию; --method не меняет λ
или ограничение ячейкой. [Математика и пример API](../docs/mathematics.md).

~~~powershell
ctest --test-dir build -R "^(ExtractionMethodsTest|InvalidMethodCLI|Compare.*)$" --output-on-failure
~~~

Это серия для объединения сфер, блока с каналами и гироида тремя методами.
Сравнительные CTest-запуски используют явные имена results/SphereUnion_dc.obj,
results/SphereUnion_emilk.obj и т. д.; они отличаются от автоматических имён.
Полный ctest также запускает старые тесты/примеры, включая тор.

DC экспортирует quads, MC — triangles. На стыках min/max поле недифференцируемо.
Корректность топологии и отсутствие самопересечений не гарантированы.
[Измеренные результаты](../docs/method_comparison.md).
[План экспериментов и требования к импорту](../docs/experiment_plan.md).
