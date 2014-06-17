#ifndef CTX_ASSERT_H_
#define CTX_ASSERT_H_

#if defined(CTXCHECKS) && CTXCHECKS == 0
#  undef CTXCHECKS
#endif

#ifdef NDEBUG
#  define ASSERTSTR "ASSERTS=OFF"
#else
#  define ASSERTSTR "ASSERTS=ON"
#endif
#ifdef CTXCHECKS
#  define CHECKSTR "CHECKS=ON"
#else
#  define CHECKSTR "CHECKS=OFF"
#endif

//
// Internal Integrity Checks: ctx_check(), ctx_assert(), ctx_assume()
//
// ctx_assert() behaves like assert()
// ctx_assume() behaves like assert() when NDEBUG=1,
//              otherwise interpret as a truth statement for the compiler
// ctx_check()  behaves like assert() when CTXCHECKS=1 but for heavy checks
//              without CTXCHECKS, does nothing
//
//            | NDEBUG=1 |      Default
//----------------------------------------------
// ctx_assert | nothing  |  fast check + abort()
// ctx_assume | optimise |  fast check + abort()
//
//            | Default  |     CTXCHECKS=1
//----------------------------------------------
// ctx_check  | nothing  |  slow check + abort()

void ctx_assertf(const char *file, const char *func, int line,
                 const char *assert, const char *fmt, ...)
__attribute__((format(printf, 5, 6)))
__attribute__((noreturn));

void ctx_assertf_no_abort(const char *file, const char *func, int line,
                          const char *assert, const char *fmt, ...)
__attribute__((format(printf, 5, 6)));


// ctx_check(): slow heavy checks, deep structure traversal
#ifdef CTXCHECKS
  #define ctx_check2(x,msg,...) ctx_assert2(x,msg,##__VA_ARGS__)
#else
  #define ctx_check2(x,msg,...) do {} while(0)
#endif

// ctx_assert(): quick value tests
#ifdef NDEBUG
  #define ctx_assert2(x,msg,...) do {} while(0)
#else
  #define ctx_assert2(x,msg,...) \
((x) ? (void)0 : ctx_assertf(__FILE__,__func__,__LINE__,QUOTE_VALUE(x),msg,##__VA_ARGS__))
#endif

// ctx_assume(): tells the compiler a condition that always holds
#ifdef NDEBUG
  #define ctx_assume2(x,msg,...) do { if(x) (void)0; else __builtin_unreachable(); } while(0)
#else
  #define ctx_assume2(x,msg,...) ctx_assert2(x,msg,##__VA_ARGS__)
#endif

// Check is turned on with CTXCHECKS=1 -> heavy lifting involved
// assert -> no action if NDEBUG=1
// assume -> declares !x impossible (helps with optimisations)
#define ctx_check(x)  ctx_check2(x,NULL)
#define ctx_assert(x) ctx_assert2(x,NULL)
#define ctx_assume(x) ctx_assume2(x,NULL)

// Return false if a condition fails, rather than aborting
// Note: doesn't depend on CTXCHECKS
#define ctx_assert_ret2(x,msg,...) do { if(!(x)) {                             \
  ctx_assertf_no_abort(__FILE__,__func__,__LINE__,QUOTE_VALUE(x),msg,##__VA_ARGS__);\
  return false;                                                                \
}} while(0)

#define ctx_assert_ret(x) ctx_assert_ret2(x,NULL)


#endif /* CTX_ASSERT_H_ */
