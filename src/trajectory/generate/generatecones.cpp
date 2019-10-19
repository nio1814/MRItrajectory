extern "C"
{
#include "cones.h"
}

int main()
{
  generateCones(28, 14, nullptr, 2, 2, 32, 1, InterConeCompensation::Compensation1, 5e-3, 4e-6, 28, 4, 15000, WaveformStorageType::STORE_BASIS);
}
