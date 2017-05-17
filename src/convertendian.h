#ifndef CONVERTENDIAN_H
#define CONVERTENDIAN_H

enum Endian{LittleEndian, BigEndian};

void swapArrayEndian(void* buffer, int length, int ptsize);

int needEndianSwap(enum Endian endianDesired);

#endif // ENDIAN_H
