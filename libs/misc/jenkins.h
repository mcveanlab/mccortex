#ifndef JENKINS_H_
#define JENKINS_H_

// 5 ops per byte
static inline uint32_t jenkins_mix(uint32_t h, uint8_t x) {
  h += x; h += (h<<10); h ^= (h>>6); return h;
}

static inline uint32_t jenkins_finish(uint32_t h) {
  h += (h<<3); h ^= (h>>11); h += (h<<15); return h;
}

// 5*bytes+6 ops [32bit => 26, 64 => 46]
static inline uint32_t jenkins_one_at_a_time_hash(const uint8_t *key, size_t len)
{
  uint32_t hash, i;
  for(hash = i = 0; i < len; ++i) hash = jenkins_mix(hash, key[i]);
  return jenkins_finish(hash);
}

// 2 ops per byte
#define strhash_fast_mix(h,x) ({ (h) = (h) * 37 + (x); (h); })

#endif /* JENKINS_H_ */
