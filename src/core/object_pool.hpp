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

#ifndef RVL_CORE_OBJECT_POOL_HPP
#define RVL_CORE_OBJECT_POOL_HPP

#include "memory_allocator.hpp"
#include <vector>
#include <memory>
#include <bitset>

namespace rvl {
namespace core {
namespace memory {

template<typename T, size_t BlockSize = 1024>
class object_pool {
public:
    object_pool() = default;

    ~object_pool() {
        for (auto& block : blocks_) {
            for (size_t i = 0; i < BlockSize; ++i) {
                if (block->used[i]) {
                    T* ptr = reinterpret_cast<T*>(&block->storage[i * sizeof(T)]);
                    ptr->~T();
                }
            }
        }
    }

    object_pool(const object_pool&) = delete;
    object_pool& operator=(const object_pool&) = delete;

    template<typename... Args>
    T* allocate(Args&&... args) {
        if (!available_.empty()) {
            T* ptr = available_.back();
            available_.pop_back();
            new(ptr) T(std::forward<Args>(args)...);
            ++allocated_;
            return ptr;
        }

        for (auto& block : blocks_) {
            if (block->used_count < BlockSize) {
                T* ptr = allocate_from_block(*block);
                if (ptr) {
                    new(ptr) T(std::forward<Args>(args)...);
                    ++allocated_;
                    return ptr;
                }
            }
        }

        blocks_.emplace_back(std::make_unique<block>());
        T* ptr = allocate_from_block(*blocks_.back());
        new(ptr) T(std::forward<Args>(args)...);
        ++allocated_;
        return ptr;
    }

    void deallocate(T* ptr) {
        if (!ptr) return;

        ptr->~T();
        available_.push_back(ptr);
        --allocated_;

        for (auto& block : blocks_) {
            char* block_start = reinterpret_cast<char*>(&block->storage[0]);
            char* block_end = block_start + sizeof(block->storage);
            char* ptr_char = reinterpret_cast<char*>(ptr);

            if (ptr_char >= block_start && ptr_char < block_end) {
                size_t index = (ptr_char - block_start) / sizeof(T);
                block->used[index] = false;
                --block->used_count;
                break;
            }
        }
    }

    size_t allocated_count() const noexcept { return allocated_; }
    size_t available_count() const noexcept { return available_.size(); }

private:
    struct block {
        alignas(T) char storage[sizeof(T) * BlockSize];
        std::bitset<BlockSize> used;
        size_t used_count = 0;
    };

    std::vector<std::unique_ptr<block>> blocks_;
    std::vector<T*> available_;
    size_t allocated_ = 0;

    T* allocate_from_block(block& b) {
        for (size_t i = 0; i < BlockSize; ++i) {
            if (!b.used[i]) {
                b.used[i] = true;
                ++b.used_count;
                return reinterpret_cast<T*>(&b.storage[i * sizeof(T)]);
            }
        }
        return nullptr;
    }
};

} 
} 
} 

#endif 