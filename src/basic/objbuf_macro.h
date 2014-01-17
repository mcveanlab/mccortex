#ifndef OBJBUF_H_
#define OBJBUF_H_

//
// objbuf.h
// Define a buffer with functions to alloc, dealloc, resize, add, append, reset
//
// Example:
//
//   create_objbuf(charbuf,String,char)
//
// Creates:
//
//   typedef struct {
//     char *data;
//     size_t len, capacity;
//   } String;
//
//   charbuf_alloc(String *buf, size_t capacity)
//   charbuf_dealloc(String *buf)
//   charbuf_ensure_capacity(String *buf, size_t capacity)
//   charbuf_add(String *buf, char obj)
//   charbuf_append(String *buf, char *obj, size_t n)
//   charbuf_reset(String *buf)
//

// Round a number up to the nearest number that is a power of two
#ifndef roundup2pow
  #define roundup2pow(x) (1UL << (64 - leading_zeros(x)))
#endif

#define create_objbuf(FUNC,buf_t,obj_t)                                        \
                                                                               \
typedef struct {                                                               \
  obj_t *data;                                                                 \
  size_t len, capacity;                                                        \
} buf_t;                                                                       \
                                                                               \
static inline void FUNC ## _alloc(buf_t *buf, size_t capacity) {               \
  buf->capacity = roundup2pow(capacity);                                       \
  buf->data = malloc2(sizeof(obj_t) * buf->capacity);                          \
  buf->len = 0;                                                                \
}                                                                              \
                                                                               \
static inline void FUNC ## _dealloc(buf_t *buf) {                              \
  free(buf->data);                                                             \
}                                                                              \
                                                                               \
static inline void FUNC ## _ensure_capacity(buf_t *buf, size_t cap) {          \
  if(cap > buf->capacity) {                                                    \
    buf->capacity = roundup2pow(cap);                                          \
    buf->data = realloc2(buf->data, sizeof(obj_t) * buf->capacity);            \
  }                                                                            \
}                                                                              \
                                                                               \
static inline void FUNC ## _add(buf_t *buf, obj_t obj) {                       \
  FUNC ## _ensure_capacity(buf, buf->len+1);                                   \
  buf->data[buf->len++] = obj;                                                 \
}                                                                              \
                                                                               \
static inline void FUNC ## _append(buf_t *buf, obj_t *obj, size_t n) {         \
  FUNC ## _ensure_capacity(buf, buf->len+n);                                   \
  memcpy(buf->data+buf->len, obj, n*sizeof(obj_t));                            \
  buf->len += n;                                                               \
}                                                                              \
                                                                               \
static inline void FUNC ## _reset(buf_t *buf) { buf->len = 0; }                \


#endif /* OBJBUF_H_ */
