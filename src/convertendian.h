#ifndef CONVERTENDIAN_H
#define CONVERTENDIAN_H

#include <stdio.h>

enum Endian{LittleEndian, BigEndian};

void swapArrayEndian(void* buffer, int length, int ptsize);

int needEndianSwap(enum Endian endianDesired);

void writeArray(void* array, unsigned long points, int pointSize, FILE* file);

void readArray(void** array, int pointSize, FILE* file, enum Endian endian);

#endif // ENDIAN_H
