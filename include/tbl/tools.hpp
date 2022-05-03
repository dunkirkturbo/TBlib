#ifndef TBLIB_TOOLS_HPP
#define TBLIB_TOOLS_HPP

#include <cstdlib>
#include <utility>

namespace tbl {

template<class T, size_t align, class... Args>
T* alloc_aligned(size_t n, Args&& ...args) {
    T* ret = reinterpret_cast<T*>(_aligned_malloc(n * sizeof(T), align));
    for (size_t i = 0; i < n; i++) {
        new (&ret[i]) T(std::forward<Args>(args)...);
    }
    return ret;
}

template<class T>
void free_aligned(size_t n, T* p) {
    for (size_t i = 0; i < n; i++) {
        p[i].~T();
    }
    _aligned_free(p);
}

}

#endif //TBLIB_TOOLS_HPP
