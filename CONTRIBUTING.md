# Contributing to RadioVectorLib

## Development Setup

### Prerequisites
- C++20 compatible compiler (GCC 10+, Clang 11+, MSVC 2019+)
- CMake 3.20+
- Git

### Building
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

### Running Tests
```bash
ctest --verbose
```

## Code Style Guidelines

### General Rules
1. Always consider vectorization opportunities
2. Use aligned allocations and cache-friendly layouts
3. Prefer references and views over copies
4. All resources must be managed automatically

### Naming Conventions
- Please adhere to using `PascalCase` for top level architectural components (structs, classes)
- `snake_case` is used for functions and file names
- All constants shoudl use `UPPER_SNAKE_CASE` or `snake_case` for template constants
- All member variables should use `variable_name_`
- All things Templates should use `CamelCase` or single capital letters

### File Organization
```cpp
// 1. Header guard
#ifndef RVL_MODULE_FILENAME_HPP
#define RVL_MODULE_FILENAME_HPP

// 2. System includes
#include <vector>
#include <complex>

// 3. Project includes
#include "core/memory_allocator.hpp"

// 4. Forward declarations
namespace rvl {
namespace module {
class ForwardDeclared;
}
}

// 5. Implementation
namespace rvl {
namespace module {

// Classes and functions

} // namespace module
} // namespace rvl

#endif // RVL_MODULE_FILENAME_HPP
```

### SIMD Guidelines

1. Always use Structure of Arrays (SoA) for vectorizable data
```cpp
// Good
struct ComplexArraySoA {
    aligned_buffer<float> real;
    aligned_buffer<float> imag;
};

// Bad
struct ComplexArrayAoS {
    std::vector<std::complex<float>> data;
};
```

2. Ensure proper alignment for SIMD operations
```cpp
simd_vector<float> data(1024);  // 64-byte aligned
```

3. Process data in SIMD-width chunks
```cpp
constexpr size_t simd_width = 8;
for (size_t i = 0; i < n; i += simd_width) {

}
```

### Memory Management

1. **Smart Pointers**: Use `unique_ptr` by default, `shared_ptr` only when necessary
2. **Custom Allocators**: Use provided aligned allocators for SIMD data
3. **Object Pools**: Use for frequently allocated/deallocated objects
4. **Arena Allocators**: Use for temporary calculations

### Error Handling

1. Use provided error types from `core/error.hpp`
2. Check inputs at API boundaries
3. Use `RVL_ASSERT` for debug-only checks
4. Use `RVL_REQUIRE` for release-mode checks

### Performance Considerations

1. **Minimize Allocations**: Pre-allocate buffers and reuse
2. **Cache Efficiency**: Keep working set small, access memory sequentially
3. **Branch Prediction**: Avoid unpredictable branches in hot paths
4. **Inline Functions**: Use `inline` for small, frequently called functions

### Testing Requirements

1. **Unit Tests**: Minimum 90% code coverage
2. **Performance Tests**: Benchmark against scalar baseline
3. **Accuracy Tests**: Verify numerical accuracy
4. **SIMD Tests**: Test all SIMD paths (SSE, AVX, AVX512, NEON)

### Documentation

1. **Doxygen Comments**: All public APIs must be documented
2. **Implementation Notes**: Complex algorithms need explanatory comments
3. **References**: Cite papers/standards for implemented algorithms

### Commit Guidelines

1. **Message Format**: `module: Brief description`
2. **Atomic Commits**: One logical change per commit
3. **Testing**: All commits must pass tests
4. **Sign-off**: Use `git commit -s` to sign commits

## Adding New Features

### 1. Propagation Model Example
```cpp
namespace rvl {
namespace propagation {
namespace new_model {

template<typename T>
class new_propagation_model {
public:
    using value_type = T;
    using vector_type = core::memory::simd_vector<T>;

    void calculate_batch(const vector_type& distances,
                        const vector_type& frequencies,
                        vector_type& path_loss);

    T calculate_scalar(T distance, T frequency) const;
};

} // namespace new_model
} // namespace propagation
} // namespace rvl
```

### 2. Implementation Checklist
- [ ] Header-only implementation in appropriate module
- [ ] SIMD optimized version
- [ ] Scalar fallback
- [ ] Unit tests with accuracy verification
- [ ] Performance benchmark
- [ ] Documentation with references
- [ ] Example usage

## Review Process

1. **Self Review**: Run tests, benchmarks, and linters
2. **Code Review**: All code must be reviewed before merge
3. **Performance Review**: Benchmarks must show improvement or no regression
4. **Documentation Review**: Ensure docs are complete and accurate