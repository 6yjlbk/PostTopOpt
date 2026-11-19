// Reference implementation for the proposed vtkop contract.
// Requires a compatible ParaView SDK, not standalone VTK.
// Integrate these wrappers into the real vtkop; no application pipeline here.
#include "VtkopContract.h"
#include <vtkAMRDualContour.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkCompositeDataIterator.h>
#include <vtkDataObject.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkDummyController.h>
#include <vtkMultiProcessController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkNonOverlappingAMR.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUniformGrid.h>
#include <vtkVersionMacros.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace postop::vtkop_contract {
struct UniformGridObject {
    vtkSmartPointer<vtkUniformGrid> data = vtkSmartPointer<vtkUniformGrid>::New(); // Создание объекта сетки VTK
};
struct AmrObject {
    vtkSmartPointer<vtkNonOverlappingAMR> data =
        vtkSmartPointer<vtkNonOverlappingAMR>::New();
    UniformGrid block; // block сохраняет владеющий указатель на входную сетку, чтобы она не была удалена раньше времени.
};
struct ContourObject {
    vtkSmartPointer<vtkAMRDualContour> filter = vtkSmartPointer<vtkAMRDualContour>::New(); // Создание объектв VTK-фильтра для DC
    Amr input;
    bool updated = false; // Служебный флаг выполнения работы алгоритма
};
struct PolyDataObject {
    vtkSmartPointer<vtkPolyData> data; // Создание объектв для хранения данных VTK
    std::vector<std::array<std::size_t, 3>> faces; // faces - массив индексов треугольников faces[0]={1,2,3} и тд
    bool triangulated = false; // Служебный флаг выполнения работы алгоритма
};
// Безымянное пространство имён, в котором находится вспомогательный код (проверка, работа с ошибками, работа с форматами данных)
namespace {
template<class T> T& checked(const std::shared_ptr<T>& object) { // Проверка, что указатель на объект не пуст
    if (!object) throw std::invalid_argument("Null vtkop handle");
    return *object;
}
// Перевод размеров для VTK
vtkIdType vtkSize(std::size_t value) {
    if (static_cast<std::uintmax_t>(value) >
        static_cast<std::uintmax_t>((std::numeric_limits<vtkIdType>::max)()))
        throw std::runtime_error("Size exceeds vtkIdType range");
    return static_cast<vtkIdType>(value);
}
// Обратное преобразование размеров
std::size_t nativeSize(vtkIdType value) {
    if (value < 0 || static_cast<std::uintmax_t>(value) >
        static_cast<std::uintmax_t>((std::numeric_limits<std::size_t>::max)()))
        throw std::runtime_error("Invalid VTK size or index");
    return static_cast<std::size_t>(value);
}
class Errors final : public vtkCommand { // Запоминает ошибку, чтобы после вызова фильтра сообщить о ней
public:
    static Errors* New() { return new Errors; }
    bool failed = false;
    void Execute(vtkObject*, unsigned long, void*) override { failed = true; }
};
struct ObserveErrors { // Автоматически подключает и отключает обработчик ошибок
    vtkObject* object;
    unsigned long tag;
    ObserveErrors(vtkObject* o, vtkCommand* errors)
        : object(o), tag(o->AddObserver(vtkCommand::ErrorEvent, errors)) {}
    ~ObserveErrors() { object->RemoveObserver(tag); }
};
struct ControllerScope { // Временно меняет глобальный контроллер и гарантирует восстановление предыдущего
    vtkSmartPointer<vtkMultiProcessController> previous;
    vtkNew<vtkDummyController> serial;
    ControllerScope() : previous(vtkMultiProcessController::GetGlobalController()) {
        vtkMultiProcessController::SetGlobalController(serial);
    }
    ~ControllerScope() { vtkMultiProcessController::SetGlobalController(previous); }
};
} // namespace

// Создание регулярной сетки
UniformGrid createUniformGrid(const std::array<double, 3>& origin,
                              const std::array<double, 3>& spacing,
                              const std::array<int, 3>& cells,
                              std::size_t valueCount) {
    std::array<int, 3> dimensions{};
    std::size_t cellCount = 1;
    std::size_t points = 1;
    for (std::size_t d = 0; d < 3; ++d) {
        // Не менее 3 ячеек (4 узлов) по оси; начало и шаг конечны, шаг положителен.
        if (cells[d] < 3 || cells[d] == (std::numeric_limits<int>::max)() ||
            !std::isfinite(origin[d]) || 
            !std::isfinite(spacing[d]) || spacing[d] <= 0)
            throw std::invalid_argument("Invalid regular grid");
        dimensions[d] = cells[d] + 1;
        const auto cellSize = static_cast<std::size_t>(cells[d]);
        const auto n = static_cast<std::size_t>(dimensions[d]);
        // Защита от переполнения произведения
        if (points > (std::numeric_limits<std::size_t>::max)() / n ||
            cellCount > (std::numeric_limits<std::size_t>::max)() / cellSize)
            throw std::invalid_argument("Grid size overflow");
        points *= n;
        cellCount *= cellSize;
        const double upper = origin[d] + spacing[d] * (dimensions[d] - 1); // координата верхней границы области
        // Проверка вычисленной границы
        if (!std::isfinite(upper) || upper <= origin[d])
            throw std::invalid_argument("Invalid grid bounds");
    }
    vtkSize(points);
    if (valueCount != cellCount)
        throw std::invalid_argument("VF count does not match the grid");

    // Создание сетки с заданными параметрами
    auto result = std::make_shared<UniformGridObject>();
    result->data->SetOrigin(origin[0], origin[1], origin[2]);
    result->data->SetSpacing(spacing[0], spacing[1], spacing[2]);
    result->data->SetDimensions(dimensions[0], dimensions[1], dimensions[2]);
    return result;
}
// Заполнение ячеек значениями VF
void setCellDensities(const UniformGrid& grid,
                    const double* values, std::size_t count) {
    auto& g = checked(grid);
    if (!values ||
        count != nativeSize(g.data->GetNumberOfCells()))
        throw std::invalid_argument("Invalid cell scalar array");
    vtkNew<vtkDoubleArray> array;
    array->SetNumberOfComponents(1); // Одна компонента — одно число на ячейку.
    array->SetNumberOfTuples(vtkSize(count)); // количество компонент(tuples) = кол-во КЭ
    for (std::size_t i = 0; i < count; ++i) {
        if (!std::isfinite(values[i]) || values[i] < 0 || values[i] > 1)
            throw std::invalid_argument("VF must be finite and in [0,1]");
        array->SetValue(static_cast<vtkIdType>(i), values[i]);
    }
    g.data->GetCellData()->SetScalars(array);
}
// Создания AMR-контейнера со значениями объемных долей
Amr createSingleBlockAmr(const UniformGrid& grid) {
    checked(grid);
    auto result = std::make_shared<AmrObject>();
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9,6,0)
    result->data->Initialize(std::vector<unsigned int>{1});
#else
    int blocks[1] = {1};
    result->data->Initialize(1, blocks);
#endif
    result->block = grid;
    result->data->SetDataSet(0, 0, grid->data);
    return result;
}
// Создание объекта фильтра DC
Contour createAmrDualContour() { return std::make_shared<ContourObject>(); }
// Передача ячеек с виртуалной плотностью
void setDensityInput(const Contour& contour, const Amr& input) {
    auto& c = checked(contour); // исправить, нельзя одну букву
    auto& amr = checked(input);
    auto& block = checked(amr.block);
    if (!block.data->GetCellData()->GetScalars())
        throw std::invalid_argument("Input cell densities are missing");
    c.updated = false;
    c.input = input;
    c.filter->SetInputDataObject(amr.data);
    c.filter->SetInputArrayToProcess(0, 0, 0,
        vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS);
}
// Задание уровня Iso-поверхности
void setIsovalue(const Contour& contour, double value) {
    auto& c = checked(contour);
    if (!std::isfinite(value) || value < 0 || value > 1)
        throw std::invalid_argument("Invalid VF threshold");
    c.updated = false;
    c.filter->SetIsoValue(value);
}
// Вкл/выкл закрытие поверхности на границе расчетной области
void setCapping(const Contour& contour, bool enabled) {
    auto& c = checked(contour); c.updated = false;
    c.filter->SetEnableCapping(enabled);
}
// Вкл/выкл объединение общих точек выходной сетки
void setMergePoints(const Contour& contour, bool enabled) {
    auto& c = checked(contour); c.updated = false;
    c.filter->SetEnableMergePoints(enabled);
}
// Запуск алгоритма
void update(const Contour& contour) {
    auto& c = checked(contour);
    c.updated = false;
    if (!c.input || !c.input->block ||
        !c.input->block->data->GetCellData()->GetScalars())
        throw std::runtime_error("DC input or cell densities are missing");
    // Serializes THIS implementation only. The real vtkop must coordinate
    // this scope with other host users of the global VTK controller.
    static std::mutex controllerMutex;
    const std::lock_guard<std::mutex> lock(controllerMutex);
    ControllerScope controller;
    c.filter->SetController(controller.serial);
    c.filter->SetEnableMultiProcessCommunication(1);
    vtkNew<Errors> errors;
    ObserveErrors observer(c.filter, errors);
    c.filter->Update();
    if (errors->failed) throw std::runtime_error("vtkAMRDualContour failed");
    if (!vtkMultiBlockDataSet::SafeDownCast(c.filter->GetOutputDataObject(0)))
        throw std::runtime_error("Unexpected DC output type");
    c.updated = true;
}
// Контейнер с результатами DC (список полигональных блоков результата работы алгоритма DC)
std::vector<PolyData> polyDataBlocks(const Contour& contour) {
    auto& c = checked(contour);
    if (!c.updated) throw std::runtime_error("DC must be updated before reading output");
    auto* output = vtkMultiBlockDataSet::SafeDownCast(c.filter->GetOutputDataObject(0));
    if (!output) throw std::runtime_error("Missing DC output");
    std::vector<PolyData> result;
    auto iterator = vtkSmartPointer<vtkCompositeDataIterator>::Take(output->NewIterator());
    iterator->SkipEmptyNodesOn();
    for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal(); iterator->GoToNextItem()) {
        auto* poly = vtkPolyData::SafeDownCast(iterator->GetCurrentDataObject());
        if (!poly) continue;
        auto item = std::make_shared<PolyDataObject>();
        item->data = poly;
        result.push_back(item);
    }
    return result;
}
// Триангулирование результатов работы алгоритма (полигональная модель -> триангулирование представление)
PolyData triangulate(const PolyData& input) {
    auto& source = checked(input);
    vtkNew<vtkTriangleFilter> filter;
    vtkNew<Errors> errors;
    ObserveErrors observer(filter, errors);
    filter->SetInputData(source.data);
    filter->PassLinesOff();
    filter->PassVertsOff();
    filter->Update();
    if (errors->failed || !filter->GetOutput())
        throw std::runtime_error("VTK triangulation failed");
    auto result = std::make_shared<PolyDataObject>();
    result->data = filter->GetOutput();
    result->faces.reserve(nativeSize(result->data->GetNumberOfPolys()));
    auto* cells = result->data->GetPolys();
    vtkIdType count = 0;
    const vtkIdType* ids = nullptr;
    cells->InitTraversal();
    while (cells->GetNextCell(count, ids)) {
        if (count != 3) throw std::runtime_error("Non-triangle after triangulation");
        std::array<std::size_t, 3> face{};
        for (int j = 0; j < 3; ++j) {
            if (ids[j] < 0 || ids[j] >= result->data->GetNumberOfPoints())
                throw std::runtime_error("Invalid VTK face index");
            face[j] = nativeSize(ids[j]);
        }
        result->faces.push_back(face);
    }
    result->triangulated = true;
    return result;
}
std::size_t pointCount(const PolyData& data) {
    return nativeSize(checked(data).data->GetNumberOfPoints());
}
std::size_t triangleCount(const PolyData& data) {
    auto& p = checked(data);
    if (!p.triangulated) throw std::runtime_error("Triangulate before reading faces");
    return p.faces.size();
}
std::array<double, 3> point(const PolyData& data, std::size_t index) {
    auto& p = checked(data);
    if (index >= pointCount(data)) throw std::out_of_range("VTK point index");
    std::array<double, 3> result{};
    p.data->GetPoint(vtkSize(index), result.data());
    for (double coordinate : result)
        if (!std::isfinite(coordinate))
            throw std::runtime_error("Non-finite surface coordinate");
    return result;
}
std::array<std::size_t, 3> triangle(const PolyData& data, std::size_t index) {
    auto& p = checked(data);
    if (!p.triangulated) throw std::runtime_error("Triangulate before reading faces");
    return p.faces.at(index);
}
} // namespace postop::vtkop_contract
