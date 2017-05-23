#include <endian.h>

int 
is_bigendian(void)
{
    int test_bigendian = 1;
    return (0 == *(char *)&test_bigendian);
}


void
swap_endian(char *buf, size_t wordsize, size_t n_words)
{
#define CHAR_SWAP(a,b) do { char x ; x=(a) ; (a)=(b) ; (b)=x ; } while(0)
    switch(wordsize) {
    case 2:
        for( ; n_words-- ; buf += wordsize) {
            CHAR_SWAP(buf[0], buf[1]);
        }
        break;
    case 4:
        for( ; n_words-- ; buf += wordsize) {
            CHAR_SWAP(buf[0], buf[3]);
            CHAR_SWAP(buf[1], buf[2]);
        } 
        break;
    case 8:
        for( ; n_words-- ; buf += wordsize) {
            CHAR_SWAP(buf[0], buf[7]);
            CHAR_SWAP(buf[1], buf[6]);
            CHAR_SWAP(buf[2], buf[5]);
            CHAR_SWAP(buf[3], buf[4]);
        }
        break;
    case 16:
        for( ; n_words-- ; buf += wordsize) {
            CHAR_SWAP(buf[0], buf[15]);
            CHAR_SWAP(buf[1], buf[14]);
            CHAR_SWAP(buf[2], buf[13]);
            CHAR_SWAP(buf[3], buf[12]);
            CHAR_SWAP(buf[4], buf[11]);
            CHAR_SWAP(buf[5], buf[10]);
            CHAR_SWAP(buf[6], buf[ 9]);
            CHAR_SWAP(buf[7], buf[ 8]);
        }
        break;
    case 1: 
    default:
        return;
    }
#undef CHAR_SWAP
}

