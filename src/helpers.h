#ifndef _HELPERS_H_
#define _HELPERS_H_

#define unlikely(x) (__builtin_expect ((x), 0))
#define likely(x) (__builtin_expect ((x), 1))

template<class T>
auto operator<<(std::ostream& out, const T& t) -> decltype(t.print_to(out)) {
  // auto f()->decltype() used so that template only applies to
  // things with proper print_to method
  return t.print_to(out);
}

#endif // _HELPERS_H_
