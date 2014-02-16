#ifndef MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9
#define MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9

typedef struct sha256_context SHA256_Context;
typedef struct { unsigned char v[32]; } SHA256_Sum;

SHA256_Context *sha256_create(lua_State *L);
void sha256_destroy(SHA256_Context *p);
void sha256_update(SHA256_Context *c, const void *ptr, unsigned int size);
void sha256_sum(SHA256_Sum *r, SHA256_Context *c);
void sha256_reset(SHA256_Context *c);
int  sha256_cmp(const SHA256_Sum *a, const SHA256_Sum *b);

void sha256_sum_add_string(SHA256_Context *p, const char *ptr, unsigned int count);
void sha256_sum_add_ints(SHA256_Context *p, const int *ptr, unsigned int count);
void sha256_sum_add_doubles(SHA256_Context *p, const double *ptr, unsigned int count);
void sha256_sum_add_floats(SHA256_Context *p, const float *ptr, unsigned int count);

void sha256_sum_string(SHA256_Sum *r, const char *ptr, unsigned int count);

void sha256_sum_clear(SHA256_Sum *r);

#endif /* !defined(MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9) */
