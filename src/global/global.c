#include "global.h"

#include <time.h>
#include <sys/time.h> // for seeding random
#include <unistd.h> // getpid()

#include "ctx_output.h" // ctx_output_init()
#include "misc/jenkins.h" // hash functions

#define rotl32(h,r) ((h)<<(r)|(h)>>(32-(r)))

static inline uint32_t get_rand_seed()
{
  struct timeval now;
  gettimeofday(&now, NULL);

  uint32_t h = rand();
  h = strhash_fast_mix(h, rotl32((uint32_t)now.tv_sec,  h & 31));
  h = strhash_fast_mix(h, rotl32((uint32_t)now.tv_usec, h & 31));
  h = strhash_fast_mix(h, (uint32_t)getpid());
  return h;
}

void seed_random()
{
  uint32_t seed = get_rand_seed();
  srand(seed);
  srand48(~seed);
}

void cortex_init()
{
  seed_random();
  // Cannot use die/warn/message/timestamp until we have completed setup
  ctx_output_init();
  // Now safe to use die/warn/message/timestamp methods
  // since mutex and cmdcode have been set
}

void cortex_destroy()
{
  ctx_output_destroy();
}
