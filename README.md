# PostTopOpt

C++20-библиотека извлечения поверхности для постобработки топологической оптимизации. Eigen используется для линейной алгебры; CGAL, VTK и ParaView подключаются опционально. Методы: dc, emilk, mc, cgal-dc, paraview-amr.

## Единый запуск

Программа `PostTopOptRun` принимает PLY, имя встроенной аналитической фигуры либо сохранённые плотности CSV/PTOV1. Все источники приводятся к одной регулярной сетке с плотностями элементов. DC реконструирует поле трилинейно между центрами; ParaView получает cell data напрямую. Аналитический Field по-прежнему доступен через API библиотеки и старые примеры.

```powershell
$env:PATH="C:\msys64\ucrt64\bin;$env:PATH"
cmake --preset release
cmake --build --preset release
.\build\release\PostTopOptRun.exe --input data/models/Armadillo.ply --cells 64 --method compare --benchmark true --output-dir results/armadillo-64
.\build\release\PostTopOptRun.exe --input gyroid --cells 64 --complexity 4 --method compare --benchmark true --output-dir results/gyroid-64 --save-volume data/volumes/gyroid-64.csv
.\build\release\PostTopOptRun.exe --input data/volumes/gyroid-64.csv --method dc --output-dir results/gyroid-loaded
```

Preset release настроен на установленный MSYS2 UCRT64 в C:/msys64. Для другой машины пути и доступность SDK следует настроить отдельно. Eigen 3.4.1 загружается при конфигурации. Сборки и программы после этой реорганизации не запускались по указанию пользователя.

## Структура

| Папка | Назначение |
| --- | --- |
| include, src | Публичный API и реализация единой библиотеки |
| apps | Единый CLI, аналитические модели, ввод/вывод плотностей |
| apps/legacy | Сохранённые прежние примеры и специальные эксперименты |
| data/models | Исходные полигональные модели |
| data/volumes | Сохранённые массивы плотностей |
| results | Отдельная папка на запуск: настройки, CSV и OBJ |
| build/release, build/debug | Новые сборки через presets |
| build/legacy | Архив старых кэшей и промежуточных файлов; не использовать для сборки |
| tests, docs, third_party | Тесты, описание и лицензии |

Обычные результаты и большие новые данные исключены из Git. Ранее версионированные PLY сохраняются при переносе. Удалённые пользователем результаты не восстанавливаются.

## Параметры и проверка

Рабочие параметры DC: lambda=0.001, regularization=true, constrained=true, bisection=true; все изменяемы. `--benchmark true` включён по умолчанию и включает прогрев, повторные замеры и отдельный профиль. `--export true` включён по умолчанию. При `--benchmark false` каждый метод выполняется один раз; summary.csv и runs.csv не создаются. Выходной каталог должен быть новым или пустым.

Команды, формат данных и смысл метрик: [docs/running.md](docs/running.md). Математика: [docs/mathematics.md](docs/mathematics.md). Архитектура: [docs/architecture.md](docs/architecture.md). Предыдущие измерения: [docs/benchmark.md](docs/benchmark.md).

## Ограничения и происхождение

Плотности, полученные из PLY и аналитических фигур, синтетические: это выборка поля в центрах элементов, а не точное интегрирование занятого объёма. Гарантии manifold, отсутствия самопересечений и автоматического закрытия AABB не предоставляются. Качество разных backend требует отдельной оценки; одинаковое число элементов не доказывает эквивалентность поверхностей.

Базовый DC написан самостоятельно по постановке, сопоставленной с BorisTheBrave/mc-dc; это не дословный перенос Python. Emilk-подобный вариант сохраняет общий собственный решатель. Таблица MC взята из scikit-image v0.25.2 с сохранением [атрибуции и лицензии](third_party/README.md). Runtime-зависимости от Python нет.
