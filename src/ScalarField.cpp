#include "ScalarField.h"

namespace geometry
{
    ScalarField::ScalarField(std::size_t size)
        : values_(size, 0.0)
    {
    }

    void ScalarField::setValue(
        std::size_t index, double value)
    {
        values_.at(index) = value;
    }

    double ScalarField::getValue(
        std::size_t index) const
    {
        return values_.at(index);
    }

    std::size_t ScalarField::size() const
    {
        return values_.size();
    }
}
