extern "C"
{
#include "cones.h"
}

int main()
{
  const Cones* cones = generateCones(28, 14, nullptr, 2, 2, 0, 1, InterConeCompensation::Compensation1, 5e-3, 4e-6, 28, 4, 15000, WaveformStorageType::STORE_BASIS);
  fprintf(stdout, "%d", cones->numCones);
}
