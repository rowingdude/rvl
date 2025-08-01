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

#ifndef RVL_CORE_ARENA_HPP
#define RVL_CORE_ARENA_HPP

#include "memory_allocator.hpp"
#include <utility>
#include <cstring>

namespace rvl {
namespace core {
namespace memory {

class arena {
public:
    explicit arena(size_t initial_size = 1024 * 1024) {
        if (initial_size > 0) {
            data_ = static_cast<char*>(aligned_alloc(initial_size, cache_line_size));
            if (!data_) {
                throw std::bad_alloc();
            }
            capacity_ = initial_size;
        }
    }

    ~arena() {
        if (data_) {
            aligned_free(data_);
        }
    }

    arena(const arena&) = delete;
    arena& operator=(const arena&) = delete;

    arena(arena&& other) noexcept
        : data_(other.data_), capacity_(other.capacity_), used_(other.used_) {
        other.data_ = nullptr;
        other.capacity_ = 0;
        other.used_ = 0;
    }

    arena& operator=(arena&& other) noexcept {
        if (this != &other) {
            if (data_) {
                aligned_free(data_);
            }
            data_ = other.data_;
            capacity_ = other.capacity_;
            used_ = other.used_;
            other.data_ = nullptr;
            other.capacity_ = 0;
            other.used_ = 0;
        }
        return *this;
    }

    void* allocate(size_t size, size_t alignment = alignof(std::max_align_t)) {
        size_t aligned_used = (used_ + alignment - 1) & ~(alignment - 1);
        size_t required = aligned_used + size;

        if (required > capacity_) {
            grow(required);
            aligned_used = (used_ + alignment - 1) & ~(alignment - 1);
        }

        void* ptr = data_ + aligned_used;
        used_ = aligned_used + size;
        return ptr;
    }

    template<typename T, typename... Args>
    T* construct(Args&&... args) {
        void* ptr = allocate(sizeof(T), alignof(T));
        return new(ptr) T(std::forward<Args>(args)...);
    }

    void reset() noexcept {
        used_ = 0;
    }

    size_t used() const noexcept { return used_; }
    size_t capacity() const noexcept { return capacity_; }

private:
    char* data_ = nullptr;
    size_t capacity_ = 0;
    size_t used_ = 0;

    void grow(size_t min_size) {
        size_t new_capacity = capacity_ * 2;
        if (new_capacity < min_size) {
            new_capacity = min_size;
        }

        char* new_data = static_cast<char*>(aligned_alloc(new_capacity, cache_line_size));
        if (!new_data) {
            throw std::bad_alloc();
        }

        if (data_ && used_ > 0) {
            std::memcpy(new_data, data_, used_);
            aligned_free(data_);
        }

        data_ = new_data;
        capacity_ = new_capacity;
    }
};

} 
} 
} 

#endif 