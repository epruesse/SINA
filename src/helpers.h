#define unlikely(x) (__builtin_expect ((x), 0))
#define likely(x) (__builtin_expect ((x), 1))
