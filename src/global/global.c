#include "global.h"

#include <time.h>
#include <sys/time.h> // for seeding random
#include <unistd.h> // getpid()

#include "ctx_output.h" // ctx_output_init()
#include "misc/jenkins.h" // hash functions

void seed_random()
{
  struct timeval now;
  gettimeofday(&now, NULL);

  uint32_t h;
  h = strhash_fast_mix(0, (uint32_t)now.tv_sec);
  h = strhash_fast_mix(h, (uint32_t)now.tv_usec);
  h = strhash_fast_mix(h, (uint32_t)getpid());

  srand(h);
  srand48(~h);
}

void cortex_init()
{
  // Cannot use die/warn/message/timestamp until we have completed setup
  ctx_output_init();
  // Now safe to use die/warn/message/timestamp methods
  // since mutex and cmdcode have been set
  seed_random();
}

void cortex_destroy()
{
  ctx_output_destroy();
}
