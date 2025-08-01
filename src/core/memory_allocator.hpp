/*
================================================================================
||                                                                            ||
||                              RadioVectorLibrary                            ||
||                                   (RVL)                                    ||
||                                                                            ||
||                                                                            ||
||                         Developer: Benjamin Cance, KC8BWS                  ||
||                         Email: kc8bws@kc8bws.com                           ||
||                         MIT License                                        ||
||                                                                            ||
||                                                                            ||
================================================================================
*/

#ifndef RVL_CORE_MEMORY_ALLOCATOR_HPP
#define RVL_CORE_MEMORY_ALLOCATOR_HPP

#include <memory>
#include <cstdlib>
#include <cstddef>
#include <limits>
#include <vector>
#include <new>

#ifdef _WIN32
    #include <malloc.h>
#endif

namespace rvl {
namespace core {
namespace memory {

constexpr size_t cache_line_size = 64;
constexpr size_t simd_alignment = 64;  

inline void* aligned_alloc(size_t size, size_t alignment) {
    #ifdef _WIN32
        return _aligned_malloc(size, alignment);
    #else
        void* ptr = nullptr;
        if (posix_memalign(&ptr, alignment, size) != 0) {
            throw std::bad_alloc();
        }
        return ptr;
    #endif
}

inline void aligned_free(void* ptr) {
    #ifdef _WIN32
        _aligned_free(ptr);
    #else
        free(ptr);
    #endif
}

template<typename T>
T* aligned_new(size_t count = 1, size_t alignment = alignof(T)) {
    const size_t size = sizeof(T) * count;
    void* ptr = aligned_alloc(size, std::max(alignment, alignof(T)));
    return static_cast<T*>(ptr);
}

template<typename T>
void aligned_delete(T* ptr) {
    if (ptr) {
        ptr->~T();
        aligned_free(ptr);
    }
}

template<typename T>
struct aligned_deleter {
    void operator()(T* ptr) const {
        if (ptr) {
            ptr->~T();
            aligned_free(ptr);
        }
    }
};

template<typename T>
using aligned_unique_ptr = std::unique_ptr<T, aligned_deleter<T>>;

template<typename T, size_t Alignment = simd_alignment>
class simd_allocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template<typename U>
    struct rebind {
        using other = simd_allocator<U, Alignment>;
    };

    simd_allocator() = default;

    template<typename U>
    simd_allocator(const simd_allocator<U, Alignment>&) noexcept {}

    T* allocate(size_type n) {
        if (n > std::numeric_limits<size_type>::max() / sizeof(T)) {
            throw std::bad_alloc();
        }

        void* ptr = aligned_alloc(n * sizeof(T), Alignment);
        if (!ptr) {
            throw std::bad_alloc();
        }

        return static_cast<T*>(ptr);
    }

    void deallocate(T* ptr, size_type) noexcept {
        aligned_free(ptr);
    }

    template<typename U, size_t A>
    bool operator==(const simd_allocator<U, A>&) const noexcept {
        return Alignment == A;
    }

    template<typename U, size_t A>
    bool operator!=(const simd_allocator<U, A>&) const noexcept {
        return Alignment != A;
    }
};

template<typename T>
using simd_vector = std::vector<T, simd_allocator<T>>;

} 
} 
} 

#endif 