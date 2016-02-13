#ifndef CRC32_H_
#define CRC32_H_

// define endianess and some integer data types
#if defined(_MSC_VER) || defined(__MINGW32__)
  typedef unsigned __int8  uint8_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int32  int32_t;

  #define __LITTLE_ENDIAN 1234
  #define __BIG_ENDIAN    4321
  #define __BYTE_ORDER    __LITTLE_ENDIAN

  #include <xmmintrin.h>
  #ifdef __MINGW32__
    #define PREFETCH(location) __builtin_prefetch(location)
  #else
    #define PREFETCH(location) _mm_prefetch(location, _MM_HINT_T0)
  #endif
#else
  // uint8_t, uint32_t, in32_t
  #include <stdint.h>
  // defines __BYTE_ORDER as __LITTLE_ENDIAN or __BIG_ENDIAN
  #include <sys/param.h>

  #ifdef __GNUC__
    #define PREFETCH(location) __builtin_prefetch(location)
  #else
    #define PREFETCH(location) ;
  #endif
#endif

/// Slicing-By-16
#define MaxSlice    16

uint32_t crc32_fast(const void* data, size_t length, uint32_t previousCrc32);

#endif/*CRC32_H_*/
