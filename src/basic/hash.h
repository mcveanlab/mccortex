#ifndef HASH_H_
#define HASH_H_

// Hash functions
#if defined(USE_CITY_HASH)
  // Use Google's CityHash
  #include "misc/city.h"
  #define HASH_NAME_STR "CityHash32"
  #define ctx_hash32(src,n,rehash) ((uint32_t)CityHash64WithSeed((char*)(src), (n), (rehash)))
  #define ctx_hash64(src,n,rehash) CityHash64WithSeed((char*)(src), (n), (rehash))
#elif defined(USE_XXHASH)
  // Use xxHash
  #include "xxHash/xxhash.h"
  #define HASH_NAME_STR "xxHash32"
  #define ctx_hash32(src,n,rehash) XXH32((src), (n), (rehash))
  #define ctx_hash64(src,n,rehash) XXH64((src), (n), (rehash))
#else
  // Use Bob Jenkin's lookup3
  #include "misc/lookup3.h"
  #define HASH_NAME_STR "Lookup3"
  #define ctx_hash32(src,n,rehash) lk3_hashlittle((src), (n), (rehash))

static inline uint64_t ctx_hash64(void *ptr, size_t n, uint64_t init)
{
  uint32_t a = init>>32, b = init;
  lk3_hashlittle2(ptr, n, &a, &b); // note: `a` slightly better mixed than `b`
  return (((uint64_t)b<<32) | a);
}

#endif

#endif /* HASH_H_ */
