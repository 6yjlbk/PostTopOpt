# Архитектура PostTopOpt: единая библиотека и backend

Статус на 12 ноября 2026: изменения подготовлены в architecture/unified-backends. Сборка, конфигурация CMake и тесты агентом не запускались. Наличие кода не подтверждает работоспособность ParaView-интеграции.

## Границы ответственности

Одна библиотечная цель `PostTopOpt`, общий API `PostTopOpt.h`. Необязательные реализации включаются опциями CGAL, VTK и ParaView. Внешние типы не входят в публичные сигнатуры. `PostTopOptVtk` — совместимое имя (ALIAS) той же цели, а не отдельная библиотека. Исходные имена старых API сохранены.

- `Field.h` — регулярная геометрия, Field и Gradient, проверка и индексация ячеек.
- `CellScalarVolume.h` — значения на ячейках, выборка поля и линейная реконструкция.
- `DualContouring.h`, `CgalExtraction.h`, `ParaViewExtraction.h` — собственные параметры backend.
- `PostTopOpt.h` — выбор backend и два входа: аналитический Field или CellScalarVolume.
- `ParaViewExtraction.cpp` — единственное место использования типов ParaView; при выключенной опции доступен явный отказ с диагностикой.
- PLY-адаптеры и OBJ сохраняют отдельную ответственность за вход/выход.

Eigen остаётся публичной зависимостью: Vec3 входит в API. Внешние реализации неизбежно добавляют зависимости линковки при включении, но их типы скрыты от пользователя библиотеки.

## Данные на элементах

CellScalarVolume владеет Grid и массивом double. Origin — физический нижний узел расчётной сетки, значения расположены в центрах ячеек. Индекс `(x*Ny+y)*Nz+z`, ось z меняется быстрее всего. validate проверяет размеры, геометрию, длину массива и конечность значений. Функция sampleCells позволяет явно получить такой массив из аналитической функции; это точечная выборка центров, не интегрирование объёмной доли по элементу.

linearCellField выполняет трилинейную интерполяцию между центрами. За пределами выпуклой оболочки центров используется постоянное продолжение по соответствующей оси. Данные не становятся узловыми значениями посредством переименования; положения отсчётов учитываются со сдвигом половины шага. Этот метод не сохраняет интегральную массу автоматически и не закрывает поверхность на физической границе расчётной области.

Встраиваемый продукт в будущем может предоставить Field/Gradient из функций формы напрямую. Ядро не требует реализации функций формы. Реконструкция центров — начальный независимый адаптер, а не утверждённая модель интерполяции конечного продукта. Лямбда удерживает shared_ptr<const CellScalarVolume>; данные нельзя изменять через другие алиасы в период использования callback.

Для высокой плотности внутри densityField задаёт $f=\tau-\rho$, изоуровень извлечения равен нулю. При входе CellScalarVolume аргумент threshold является порогом плотности и заменяет изоуровни backend. При аналитическом входе используются собственные isovalue каждого backend. Плотность и физическая объёмная доля не отождествляются автоматически.

## Параметры

ExtractionSettings содержит dc, cgal, paraview и mcIsovalue. Все значения изменяемы для каждого вызова, глобальных постоянных настроек нет. DC сохраняет исходное состояние: коэффициент 0.001, regularizationEnabled=true, constrainToCell=true, bisection=true. При regularizationEnabled=false коэффициент хранится, но не используется в QEF. Старый способ отключения коэффициентом 0 также остаётся допустимым.

CGAL использует собственную QEF, поэтому коэффициент нашего DC не передаётся. Старый диспетчер исправлен: bisection теперь передаётся явно и в обычный, и в octree-вариант CGAL. Эксперимент gradient_comparison приведён в соответствие с его назначением: бисекция действительно включена. MC не использует QEF.

## ParaView AMR

POSTTOPOPT_WITH_PARAVIEW включает реальный вызов vtkAMRDualContour из ParaView::VTKExtensionsAMR. Вход — один блок vtkUniformGrid в одном уровне vtkNonOverlappingAMR. Значения [0,1] копируются в cell data с преобразованием z-fastest → x-fastest. Используется SetIsoValue в единицах массива, без неявного умножения на 255 и без оболочки vtkPVAMRDualContour. Опции закрытия границы и слияния точек задаются отдельно; MPI-обмен выключен.

Multiblock-результат обходится по листьям, треугольники копируются в Mesh со смещением индексов. Метаданные частиц не переносятся, глобальное слияние разных блоков не выполняется. Вариант пока рассчитан на один блок; полноценный AMR и ghost-слои требуют отдельной проверки на выбранной версии ParaView. Для первого прогона использовать поле с пустым запасом около границ, затем отдельно проверить capping и положение поверхности. Топологические гарантии и поведение на границе пока не подтверждены.

Аналитический вход нельзя передать в ParaView напрямую: сначала явно вызывается sampleCells для функции плотности в [0,1]. Произвольный SDF не является объёмной долей. Собственные DC/CGAL используют восстановленное поле, ParaView — исходные значения ячеек; это различие должно фиксироваться в сравнениях.

SDK ParaView и standalone VTK не должны смешиваться. Если включены обе опции, PLY-модуль использует VTK из ParaView; наличие IOPLY проверяется. Компилятор, архитектура и конфигурация SDK должны совпадать с проектом. Установка ParaView не выполнялась.

## Пример API

```cpp
#include "PostTopOpt.h"
geometry::ExtractionSettings settings;
settings.dc.regularization = 0.001;
settings.dc.regularizationEnabled = true;
settings.dc.constrainToCell = true;
settings.dc.bisection = true;
auto mesh = geometry::extract(geometry::Backend::Dc, grid, analyticField, {}, settings);
// Для данных на элементах:
auto result = geometry::extract(geometry::Backend::Dc, volume, 0.5, settings);
// При включённом SDK:
auto comparison = geometry::extract(geometry::Backend::ParaViewAmr, volume, 0.5, settings);
```

## Ручная проверка пользователем

Исходники проверок добавлены в tests/architecture.cpp: воспроизведение линейного поля, порядок индексов, размер массива, знак плотности, сохранение параметров, аналитический API, отключённый ParaView. Они не запускались.

В отдельном каталоге с MinGW сначала проверить конфигурацию без внешних backend:

```powershell
cmake -S . -B build-architecture -G "MinGW Makefiles" -DBUILD_TESTING=ON -DPOSTTOPOPT_WITH_VTK=OFF -DPOSTTOPOPT_WITH_CGAL=OFF -DPOSTTOPOPT_WITH_PARAVIEW=OFF
cmake --build build-architecture -j 4
ctest --test-dir build-architecture --output-on-failure
.\build-architecture\PostTopOptCellDensity.exe --method dc --output results/cell-density-dc.obj
```

CellDensity принимает --lambda, --regularization 0|1, --constrained 0|1 и --bisection 0|1. Результат не перезаписывается. Затем отдельно проверить текущие PLY и CGAL-сборки. Для ParaView нужна новая конфигурация с POSTTOPOPT_WITH_PARAVIEW=ON и ParaView_DIR от SDK; использовать совместимый компилятор. В примере метод paraview-amr вызывает тот же API библиотеки.

Некоторые старые CTest-сценарии по-прежнему пишут в results/. Перенос их артефактов и обобщение CLI остаются отдельным небольшим этапом, не скрытым изменением текущего поведения.

Источники API: [vtkAMRDualContour](https://www.paraview.org/paraview-docs/latest/cxx/classvtkAMRDualContour.html), [исходники ParaView AMR](https://github.com/Kitware/ParaView/tree/master/VTKExtensions/AMR).

## Установленные компоненты — 12 ноября 2026

Через MSYS2 UCRT64 установлены ParaView 6.1.1-4 и зависимости, включая CGAL 6.2-1; VTK 9.6.2-8 уже установлена. Найдены vtkAMRDualContour.h, библиотека vtkPVVTKExtensionsAMR и CMake-конфигурация ParaView в `C:/msys64/ucrt64/lib/cmake/paraview`. Пакет ParaView использует системную VTK UCRT64. Инициализация AMR адаптирована к объявленному в заголовках VTK 9.6 интерфейсу. Это статическая сверка интерфейса, а не подтверждение сборки или исполнения.

Для ручной проверки пользователем после добавления `C:/msys64/ucrt64/bin` в PATH:

```powershell
cmake -S . -B build-paraview -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/c++.exe -DCMAKE_MAKE_PROGRAM=C:/msys64/ucrt64/bin/mingw32-make.exe -DCMAKE_PREFIX_PATH=C:/msys64/ucrt64 -DParaView_DIR=C:/msys64/ucrt64/lib/cmake/paraview -DPOSTTOPOPT_WITH_PARAVIEW=ON -DPOSTTOPOPT_WITH_VTK=OFF -DBUILD_TESTING=ON
cmake --build build-paraview -j 4
ctest --test-dir build-paraview --output-on-failure
.\build-paraview\PostTopOptCellDensity.exe --method paraview-amr --output results/cell-density-paraview.obj
```

Агент эти команды не выполнял. Для первой проверки ParaView PLY-модуль выключен, чтобы проверять интеграцию поэтапно.

## Проверка после разрешения сборок — 12 ноября 2026

Сборка build-paraview с включёнными CGAL и ParaView успешна. Повторно пройдены 4 выбранных теста: DualContouringTest, ExtractionMethodsTest, CgalBackendTest, ArchitectureApiTest. Последний теперь запускает реальное ParaView-извлечение при включённом backend и проверяет непустой результат, конечность вершин и индексы граней. Полный набор CTest в этой проверке не запускался.

Исправлены три проблемы:

- В установленной VTK вызов шаблонного SetValue не линковался; использован публичный SetTuple1.
- Регрессионный тест сравнивал линейную интерполяцию, не задавая режим; теперь явно отключает бисекцию только для этой проверки.
- AMR helper ParaView обращается к глобальному контроллеру. Адаптер временно устанавливает vtkDummyController и восстанавливает предыдущий через RAII. Вызовы адаптера сериализуются mutex. Это не делает безопасной одновременную работу постороннего VTK-кода с глобальным контроллером; такой режим встраивания не поддерживается и требует согласования с хостом.

Ручные запуски одного синтетического поля: ParaView — 1896 вершин и 3788 треугольников; DC и CGAL — по 1832 вершины и 1830 quads. OBJ сохранены только в build-paraview/verified-*.obj. Это подтверждение исполнения, а не измерение геометрической точности или топологических гарантий.
