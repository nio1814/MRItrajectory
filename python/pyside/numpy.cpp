#include "numpy.h"

#include <algorithm>


PyObject* mapNumpyArray(void *data, std::vector<npy_intp> dimensions)
{
  initializeNumpy();

  std::vector<npy_intp> strides = {};
  npy_intp stride = 1;
  for (const npy_intp dimension : dimensions) {
    strides.push_back(stride);
    stride *= dimension;
  }

  std::reverse(std::begin(dimensions), std::end(dimensions));

  return PyArray_NewFromDescr(&PyArray_Type, PyArray_DescrFromType(NPY_FLOAT32), dimensions.size(), dimensions.data(), NULL, data, NPY_ARRAY_WRITEABLE, NULL);
}
