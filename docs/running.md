# Запуск и результаты

Изменения подготовлены 13 ноября 2026. Конфигурация CMake, сборки, приложения и тесты после реорганизации не запускались по указанию пользователя. Старые команды в исторических заметках относятся к прежним путям. Актуальная точка входа — PostTopOptRun.

## Сборка

PowerShell из корня проекта:

```powershell
$env:PATH="C:\msys64\ucrt64\bin;$env:PATH"
cmake --preset release
cmake --build --preset release
# Полная пересборка:
cmake --build --preset release --clean-first
# Debug:
cmake --preset debug
cmake --build --preset debug
```

Presets требуют MSYS2 UCRT64 по пути C:/msys64, GCC, mingw32-make и установленные SDK ParaView/VTK/CGAL. build/legacy содержит перенесённые старые сборки с недействительными абсолютными путями; их не использовать. Старые файлы сохранены до проверки новой конфигурации.

## Общая команда

```powershell
.\build\release\PostTopOptRun.exe --input SOURCE --cells 64 --method dc --benchmark true --export true --output-dir results/RUN
```

SOURCE — PLY, имя аналитического задания либо CSV/PTOV1. Для файлов плотностей убрать --cells и --complexity: геометрия уже содержится в файле. Все пути считаются относительно текущего рабочего каталога. Новый запуск требует новой или пустой папки результатов.

Аналитический банк находится в apps/AnalyticModels.cpp: sphere, torus, sphere-union, drilled-box, cone, gyroid, channels, lattice, shells. Область [-1,1]^3; cells — число ячеек на каждой оси. Параметры фигур фиксированы внутри модуля; complexity=1..12 меняет частоту четырёх периодических заданий. Формулы прежних прямых аналитических примеров в apps/legacy сохранены для совместимости; они не являются тождественными по размерам новому банку.

Для PLY cells задаёт разрешение вдоль самой длинной непаддированной стороны; итоговые размеры включают запас. width-cells=(0,4] задаёт ширину перехода плотности в шагах сетки, по умолчанию 2. Открытые PLY текущий адаптер отклоняет.

## Примеры

```powershell
# DC, PLY, автоматический OBJ:
.\build\release\PostTopOptRun.exe --input data/models/Armadillo.ply --cells 32 --method dc --output-dir results/armadillo-dc-32
.\build\release\PostTopOptRun.exe --input data/models/Armadillo.ply --cells 64 --method dc --output-dir results/armadillo-dc-64

# ParaView:
.\build\release\PostTopOptRun.exe --input data/models/Armadillo.ply --cells 64 --method paraview-amr --output-dir results/armadillo-pv-64

# Сравнение только DC и ParaView, пять повторов:
.\build\release\PostTopOptRun.exe --input data/models/Armadillo.ply --cells 64 --method compare --benchmark true --repeats 5 --output-dir results/armadillo-compare

# Аналитическое задание; сохранить плотности для повторного использования:
.\build\release\PostTopOptRun.exe --input gyroid --cells 64 --complexity 4 --method compare --benchmark true --save-volume data/volumes/gyroid-64.csv --output-dir results/gyroid-64

# Повторить сравнение на сохранённых данных:
.\build\release\PostTopOptRun.exe --input data/volumes/gyroid-64.csv --method compare --benchmark true --output-dir results/gyroid-loaded

# Экспорт без повторных benchmark-измерений:
.\build\release\PostTopOptRun.exe --input data/volumes/gyroid-64.csv --method compare --benchmark true --export true --output-dir results/gyroid-meshes

# Все пять методов:
.\build\release\PostTopOptRun.exe --input sphere --cells 64 --method all --benchmark true --output-dir results/sphere-all

# Настраиваемые параметры DC:
.\build\release\PostTopOptRun.exe --input sphere --cells 64 --method dc --lambda 0.001 --regularization true --constrained true --bisection true --cache true --qef true --output-dir results/sphere-dc
```

Методы: dc, emilk, mc, cgal-dc, paraview-amr. compare — DC и ParaView; all требует все пять доступных backend. Порог материала задаётся --threshold (по умолчанию 0.5). Выключение экспорта: --export false. Булевы параметры принимают true/false и 1/0. Сохранение входных плотностей выполняется только при --save-volume; без него CSV массива автоматически не дублируется.

## Одинаковые метрики

settings.txt содержит источник, сетку и параметры. При benchmark=true (по умолчанию) runs.csv содержит каждый повтор, summary.csv — статистику и полный доступный профиль. При benchmark=false эти два CSV не создаются; каждый выбранный метод выполняется ровно один раз, независимо от --repeats. OBJ именуются dc.obj, paraview-amr.obj и т. д.

- preparation_vf_s: время с начала main до готовой сетки и валидированного массива volume fraction (VF), включая разбор аргументов, чтение/генерацию и проверки. Измеряется один раз для всех методов.
- isosurf_median_s и isosurf_min_s - isosurf_max_s: публичный extract, включая реконструкцию поля и внутренние адаптеры, до Mesh. Проверки результата, уничтожение Mesh и экспорт исключены.
- build_isosurf_s: preparation_vf_s + isosurf_median_s. Это расчётная сумма этапов, не длительность всего процесса с прогревами и повторениями.
- export_s: отдельно запись OBJ. Опциональная запись массива плотностей исключена из подготовки и извлечения.
- vertices/quads/triangles: размеры сетки. Проверяются конечность координат и допустимость индексов. Это не метрики геометрической точности или manifold.

При benchmark=true выполняются прогрев каждого метода, repeats=2..50 замеров с циклической сменой порядка методов, затем отдельный диагностический вызов для профиля и экспорта. Микротаймеры не влияют на обычные замеры. При benchmark=false выполняется один вызов каждого выбранного метода; summary.csv и runs.csv отсутствуют. Неприменимые профильные столбцы оставляются пустыми. Для MC и CGAL внутреннее профилирование пока не реализовано: доступны общее время и размеры Mesh, внутренние столбцы пусты. Для DC/emilk записываются все доступные времена и счётчики, для ParaView — подготовка VTK, фильтр и конвертация. profile_isosurf_s — полное время отдельного диагностического вызова; его этапы не обязаны суммироваться до обычной медианы.

## Формат плотностей

CSV версии 1 хранит схему, геометрию сетки и затем строки x,y,z,rho в порядке z-fastest:

```text
PostTopOptVolume,1
cells,Nx,Ny,Nz
origin,Ox,Oy,Oz
spacing,Hx,Hy,Hz
x,y,z,rho
0,0,0,0.5
...
```

Значения относятся к центрам элементов, диапазон [0,1], координаты узлов вычисляются по origin и spacing. Файл с именем фигуры — данные конкретного разрешения, не описание её формулы. Сохраняются 17 значащих цифр. Прежний PTOV1 читается для совместимости; новые массивы записываются в CSV.

## Проверки пользователем

```powershell
ctest --preset release
# Только численные и API-проверки:
ctest --test-dir build/release --output-on-failure -R '^(DualContouringTest|ExtractionMethodsTest|CgalBackendTest|ArchitectureApiTest|OptimizationRegressionTest|ApplicationDataTest)$'
```

ApplicationDataTest добавлен для проверки генерации, roundtrip CSV и отклонения повреждённого заголовка; агентом не запускался. Старые интеграционные примеры CTest сохраняют OBJ в build/release/test-results, а не в пользовательские results.

## Лимит разрешения

В едином CLI --cells принимает 16..500. Аналитический банк и загрузка CSV/PTOV1 допускают до 500³ элементов. Для PLY --cells означает разрешение вдоль длинной непаддированной стороны; внешний запас добавляется отдельно, поэтому итоговый размер по оси может быть больше 500. Ограничение разрешения PLY-адаптера — 500.

## Представление отчётов

summary.csv транспонирован: metric и description задают идентификатор и русское описание показателя; последующие столбцы — выбранные методы. Во всех временных показателях используется слово «время». Все времена в summary.csv, runs.csv и консоли выражены в секундах с тремя знаками после десятичной точки. Для CSV разделитель полей — запятая; десятичный разделитель — точка. summary.csv имеет UTF-8 BOM для чтения русских описаний в Excel.

Минимум и максимум объединены в строку isosurf_min_s - isosurf_max_s со значением вида 1.661 - 1.706. Медиана называется «Медианное время построения изоповерхности». runs.csv остаётся журналом повторов; время каждого построения — isosurf_s. Счётчики сохраняются целыми числами. Округление выполняется только при выводе: статистика вычисляется с полной точностью. Пустые ячейки означают неприменимые показатели или отключённый экспорт.

Внутренние профили API библиотеки пока хранят миллисекунды для совместимости; CLI переводит их в секунды. Исторические отчёты и legacy-программы автоматически не переписываются.

## Общее время алгоритма

В summary.csv добавлены две явные итоговые строки: algorithm_without_vf_s — общее время алгоритма без подготовки VF (медиана одного построения изоповерхности); algorithm_with_vf_s — общее время алгоритма с подготовкой VF (preparation_vf_s плюс эта медиана). Экспорт OBJ, сохранение входного CSV, прогревы, диагностический вызов и сумма повторов в эти значения не входят. По определению они совпадают с isosurf_median_s и build_isosurf_s соответственно; отдельные строки предназначены для явного итогового сравнения методов. Единицы — секунды с тремя знаками после точки.
