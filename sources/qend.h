#ifndef ENDIAN_H__89B8A0C5_E27F_44B1_B84D_C5CC43ED3D61
#define ENDIAN_H__89B8A0C5_E27F_44B1_B84D_C5CC43ED3D61
#include <stddef.h>
int is_bigendian(void);
void swap_endian(char *buf, size_t wordsize, size_t n_words);
#endif/*ENDIAN_H__89B8A0C5_E27F_44B1_B84D_C5CC43ED3D61*/
