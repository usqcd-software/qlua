#ifndef MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9
#define MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9

typedef struct sha256_context SHA256_Context;
typedef struct { unsigned char v[32]; } SHA256_Sum;

SHA256_Context *sha256_create(lua_State *L);
void sha256_destroy(SHA256_Context *p);
void sha256_update(SHA256_Context *c, void *ptr, unsigned int size);
void sha256_sum(SHA256_Sum *r, SHA256_Context *c);
void sha256_reset(SHA256_Context *c);
int  sha256_cmp(const SHA256_Sum *a, const SHA256_Sum *b);

#endif /* !defined(MARK_7096CD73_AE49_436B_9F7E_38F2EE882AE9) */
