#include <vector>
#include <cstddef>

namespace geometry
{
    class ScalarField
    {
    public:
        explicit ScalarField(std::size_t size);

        void setValue(std::size_t index, double value);
        double getValue(std::size_t index) const;

        std::size_t size() const;

    private:
        std::vector<double> values_;
    };
}