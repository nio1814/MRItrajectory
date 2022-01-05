#ifndef NUMPY_H
#define NUMPY_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>

/*!
 * \brief Set up code for using the NumPy API.
 */
inline void* initializeNumpy()
{
  Py_Initialize();
  if (PyArray_API == NULL) {
    import_array();
  }

  return NULL;
}

PyObject* mapNumpyArray(void *data, std::vector<npy_intp> dimensions);

#endif // NUMPY_H
