#include "qlua.h"                                                /* DEPS */
#include "sha256.h"                                              /* DEPS */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#define BLOCK_SIZE 64

struct sha256_context {
  lua_State   *L;
  uint32_t     h[8];
  uint8_t      x[BLOCK_SIZE];
  unsigned int x_len;
  uint64_t     len;
};

void
sha256_reset(SHA256_Context *c)
{
  c->h[0] = 0x6A09E667;
  c->h[1] = 0xBB67AE85;
  c->h[2] = 0x3C6EF372;
  c->h[3] = 0xA54FF53A;
  c->h[4] = 0x510E527F;
  c->h[5] = 0x9B05688C;
  c->h[6] = 0x1F83D9AB;
  c->h[7] = 0x5BE0CD19;
  c->x_len = 0;
  c->len = 0;
  c->L = NULL;
}

SHA256_Context *
sha256_create(lua_State *L)
{
  SHA256_Context *c = (SHA256_Context *)qlua_malloc(L, sizeof (SHA256_Context));
  sha256_reset(c);
  c->L = L;

  return c;
}

void
sha256_destroy(SHA256_Context *c)
{
  lua_State *L = c->L;
  memset(c, 0, sizeof (SHA256_Context));
  qlua_free(L, c);
}

static unsigned int
update_block(SHA256_Context *digest, uint8_t *p, unsigned int len)
{
  static const uint32_t magik[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
  };
  unsigned int lx;
  uint32_t w[64];
  uint32_t h0 = digest->h[0];
  uint32_t h1 = digest->h[1];
  uint32_t h2 = digest->h[2];
  uint32_t h3 = digest->h[3];
  uint32_t h4 = digest->h[4];
  uint32_t h5 = digest->h[5];
  uint32_t h6 = digest->h[6];
  uint32_t h7 = digest->h[7];

  for (lx = 0; lx + BLOCK_SIZE <= len; lx += BLOCK_SIZE, p += BLOCK_SIZE) {
    int i, j;
    uint32_t a, b, c, d, e, f, g, h;
    for (i = j = 0; i < 16; i++, j+= 4) {
      w[i] = ((uint32_t)p[j]) << 24 |
             ((uint32_t)p[j+1]) << 16 |
             ((uint32_t)p[j+2]) << 8 |
             ((uint32_t)p[j+3]);
    }
    for (i = 16; i < 64; i++) {
      uint32_t t1 = (w[i-2]>>17 | w[i-2]<<(32-17)) ^
                    (w[i-2]>>19 | w[i-2]<<(32-19)) ^
                    (w[i-2] >> 10);
      uint32_t t2 = (w[i-15]>>7 | w[i-15]<<(32-7)) ^
                    (w[i-15]>>18 | w[i-15]<<(32-18)) ^
                    (w[i-15] >> 3);
      w[i] = t1 + w[i-7] + t2 + w[i-16];
    }
    a = h0; b = h1; c = h2; d = h3; e = h4; f = h5; g = h6; h = h7;
    for (i = 0; i < 64; i++) {
      uint32_t t1 = h + ((e>>6 | e<<(32-6)) ^
                         (e>>11 | e<<(32-11)) ^
                         (e>>25 | e<<(32-25))) +
                  ((e & f) ^ (~e & g)) +
                   magik[i] + w[i];
      uint32_t t2 = ((a>>2 | a<<(32-2)) ^ (a>>13 | a<<(32-13)) ^ (a>>22 | a<<(32-22))) +
                    ((a & b) ^ (a & c) ^ (b & c));
      h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h0 += a; h1 += b; h2 += c; h3 += d; h4 += e; h5 += f; h6 += g; h7 += h;
  }
  digest->h[0] = h0;
  digest->h[1] = h1;
  digest->h[2] = h2;
  digest->h[3] = h3;
  digest->h[4] = h4;
  digest->h[5] = h5;
  digest->h[6] = h6;
  digest->h[7] = h7;

  return lx;
}

void
sha256_update(SHA256_Context *c, const void *p, unsigned int len)
{
  uint8_t *pp = (uint8_t *)p;
  unsigned int lx;

  c->len += len;
  if (c->x_len > 0) {
    lx = len;
    if (lx > BLOCK_SIZE - c->x_len)
      lx = BLOCK_SIZE - c->x_len;
    memcpy(c->x + c->x_len, pp, lx);
    c->x_len += lx;
    if (c->x_len == BLOCK_SIZE) {
      update_block(c, c->x, BLOCK_SIZE);
      c->x_len = 0;
    }
    pp += lx;
    len -= lx;
  }
  lx = update_block(c, pp, len);
  if (lx < len) {
    memcpy(c->x, pp + lx, len - lx);
    c->x_len = len - lx;
  }
  return;
}

void
sha256_sum(SHA256_Sum *r, SHA256_Context *c)
{
  SHA256_Context d = *c;
  uint8_t tmp[64];
  uint64_t len = d.len << 3;
  int i, j;

  memset(tmp, 0, sizeof (tmp));
  tmp[0] = 0x80;

  if (d.len % 64 < 56)
    sha256_update(&d, tmp, 56 - d.len % 64);
  else
    sha256_update(&d, tmp, 64 + 56 - d.len % 64);

  for (i = 0; i < 8; i++)
    tmp[i] = (uint8_t)(len >> (56 - 8 * i));

  sha256_update(&d, tmp, 8);

  assert(d.x_len == 0);
  for (j = i = 0; i < 8; i++) {
    r->v[j++] = (uint8_t)(d.h[i] >> 24);
    r->v[j++] = (uint8_t)(d.h[i] >> 16);
    r->v[j++] = (uint8_t)(d.h[i] >> 8);
    r->v[j++] = (uint8_t)(d.h[i]);
  }
}

int
sha256_cmp(const SHA256_Sum *a, const SHA256_Sum *b)
{
  return memcmp(a->v, b->v, 32);
}

void
sha256_sum_string(SHA256_Sum *r, const char *ptr, unsigned int count)
{
  SHA256_Context ctx;

  sha256_reset(&ctx);
  sha256_update(&ctx, ptr, count);
  sha256_sum(r, &ctx);
  memset(&ctx, 0, sizeof (SHA256_Context));
}

void
sha256_sum_add_ints(SHA256_Context *ctx, const int *ptr, unsigned int count)
{
  int i;

  for (i = 0; i < count; i++) {
    unsigned long long u = ptr[i];
    unsigned char c[sizeof (unsigned long long)];
    int j;

    for (j = 0; j < sizeof (unsigned long long); j++) {
      c[j] = (unsigned char)u;
      u >>= 8;
    }
    sha256_update(ctx, c, sizeof (unsigned long long));
  }
}

void
sha256_sum_add_doubles(SHA256_Context *ctx, const double *ptr, unsigned int count)
{
  int i;
  QLUA_ASSERT(sizeof (double) <= sizeof (unsigned long long));

  for (i = 0; i < count; i++) {
    union {
      double d;
      unsigned long long u;
    } v;
    unsigned long long u;
    unsigned char c[sizeof (unsigned long long)];
    int j;

    v.d = ptr[i];
    u = v.u;
    for (j = 0; j < sizeof (unsigned long long); j++) {
      c[j] = (unsigned char)u;
      u >>= 8;
    }
    sha256_update(ctx, c, sizeof (unsigned long long));
  }
}

void
sha256_sum_clear(SHA256_Sum *r)
{
  memset(r, 0, sizeof (SHA256_Sum));
}
