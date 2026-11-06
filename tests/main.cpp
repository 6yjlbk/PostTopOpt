#include <iostream>
#include "ScalarField.h"
#include <Eigen/Dense>

int main()
{
    geometry::ScalarField field(10);

    field.setValue(0, 0.5);
    field.setValue(1, 1.0);

    std::cout << "Field size: "
              << field.size() << '\n';

    std::cout << "Value 0: "
              << field.getValue(0) << '\n';

    std::cout << "Value 1: "
              << field.getValue(1) << '\n';

    Eigen::Vector3d v(1.0, 2.0, 3.0);

    std::cout << "Vector: "
          << v.transpose() << '\n';

    std::cout << "Length: "
          << v.norm() << '\n';

    return 0;
}
