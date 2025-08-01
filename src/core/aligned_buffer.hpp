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

#ifndef RVL_CORE_ALIGNED_BUFFER_HPP
#define RVL_CORE_ALIGNED_BUFFER_HPP

#include "memory_allocator.hpp"
#include <algorithm>
#include <stdexcept>

namespace rvl {
namespace core {
namespace memory {

template<typename T, size_t Alignment = 64>
class aligned_buffer {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;

    aligned_buffer() = default;

    explicit aligned_buffer(size_t size) {
        resize(size);
    }

    aligned_buffer(aligned_buffer&& other) noexcept
        : data_(other.data_), size_(other.size_), capacity_(other.capacity_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    aligned_buffer(const aligned_buffer&) = delete;
    aligned_buffer& operator=(const aligned_buffer&) = delete;

    aligned_buffer& operator=(aligned_buffer&& other) noexcept {
        if (this != &other) {
            deallocate();
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    ~aligned_buffer() {
        deallocate();
    }

    void resize(size_t new_size) {
        if (new_size > capacity_) {
            reserve(new_size);
        }
        size_ = new_size;
    }

    void reserve(size_t new_capacity) {
        if (new_capacity > capacity_) {
            T* new_data = aligned_new<T>(new_capacity, Alignment);
            if (data_) {
                std::copy(data_, data_ + size_, new_data);
                aligned_free(data_);
            }
            data_ = new_data;
            capacity_ = new_capacity;
        }
    }

    void shrink_to_fit() {
        if (capacity_ > size_) {
            T* new_data = nullptr;
            if (size_ > 0) {
                new_data = aligned_new<T>(size_, Alignment);
                std::copy(data_, data_ + size_, new_data);
            }
            aligned_free(data_);
            data_ = new_data;
            capacity_ = size_;
        }
    }

    T* data() noexcept { return data_; }
    const T* data() const noexcept { return data_; }

    size_t size() const noexcept { return size_; }
    size_t capacity() const noexcept { return capacity_; }

    T& operator[](size_t i) noexcept { return data_[i]; }
    const T& operator[](size_t i) const noexcept { return data_[i]; }

    T& at(size_t i) {
        if (i >= size_) {
            throw std::out_of_range("aligned_buffer::at");
        }
        return data_[i];
    }

    const T& at(size_t i) const {
        if (i >= size_) {
            throw std::out_of_range("aligned_buffer::at");
        }
        return data_[i];
    }

private:
    T* data_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;

    void deallocate() {
        if (data_) {
            aligned_free(data_);
            data_ = nullptr;
            size_ = 0;
            capacity_ = 0;
        }
    }
};

} 
} 
} 

#endif 