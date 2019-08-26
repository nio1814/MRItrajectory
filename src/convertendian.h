#ifndef CONVERTENDIAN_H
#define CONVERTENDIAN_H

#include <stdio.h>

enum Endian{LittleEndian, BigEndian};

void swapArrayEndian(void* buffer, int length, int ptsize);

int needEndianSwap(enum Endian endianDesired);

/*!
 * \brief Write an array that may or may not be initialized to file.
 * \param array The array, which may be null.
 * \param numPoints The number of elements in the array.
 * \param pointSize The number of bytes of each element.
 * \param file  The file buffer to write to.
 */
void writeArray(void* array, size_t numPoints, int pointSize, FILE* file);

void readArray(void** array, int pointSize, FILE* file, enum Endian endian);

#endif // ENDIAN_H
