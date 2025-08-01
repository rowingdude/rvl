# RadioVectorLib (RVL)

High-performance radio propagation library with SIMD/vector optimizations.

## About

While studying for my Amateur Extra exam, I realized I needed to learn some more math, to aid in this, I pulled up the EZNEC software so I could start getting an idea of the relationships these equations I was working with had with one another. I then realize the EZNEC and derivative software was written in archaric C (nothing wrong this this) and would not fully make use of the new CPU features.. this project was born from that goal. 

I am releasing this library and all future versions of it under the MIT license, I will maintain this library under active development until all of the optimization intrinsics are incorporated, thereafter, it will be developed further a new C++ features relavent to it become adopted by the standard. This library is also being written primarily as a "header" library, so it's easier to port. I've tried hard to leave out anything that would cause breaking changes to operating system ports.

## Features

- Uses widely tested CPU feature detection to enable compile-time optimizations such as vector AVX and SSE
- Maintains strict concern separation with concepts being housed in discrete file pairs
- All equations use `core::memory::simd_vector<T>` for efficient vectorized calculations
- Zero-copy, move-semantic buffer classes with aligned memory allocation
- Complete foundation for antenna design programs with current distribution, array factors, impedance analysis, and efficiency calculations
- Ionospheric ray tracing and VHF/UHF tropospheric scatter modeling
- Header-only library designed for easy porting across operating systems

## Concept Table of Contents

### Core Infrastructure
- [Constants](src/core/constants.hpp) - Physical and mathematical constants
- [Error Handling](src/core/error.hpp) - Exception classes with source location
- [Memory Management](src/core/memory_allocator.hpp) - SIMD-aligned allocators and buffers
- [Logging](src/core/logging.hpp) - Thread-safe logging infrastructure
- [Units](src/core/units.hpp) - Unit conversion utilities

### Antenna Models

#### Design Foundations
- [Current Distribution](src/antenna/current_distribution.hpp) - Current distribution for arbitrary wire antennas (Method of Moments foundation)
  - Sinusoidal and King's theory current models
  - Center-fed and off-center-fed dipole calculations
  - Monopole current distribution analysis
  - Current null detection and analysis
- [Array Factor](src/antenna/array_factor.hpp) - Array factor calculations for linear, planar, and circular arrays
  - Linear, planar, and circular array geometries
  - Dolph-Chebyshev amplitude synthesis for sidelobe control
  - Beamforming and steering calculations
  - Grating lobe analysis and prevention
  - Array directivity and beamwidth calculations
- [Impedance vs Frequency](src/antenna/impedance_frequency.hpp) - Bandwidth analysis, King's theory, resonance prediction
  - King's three-term impedance theory
  - Arbitrary length dipole impedance calculations
  - Loop antenna impedance modeling
  - Bandwidth analysis with VSWR thresholds
  - Q-factor and resonance frequency prediction
- [Radiation Efficiency](src/antenna/radiation_efficiency.hpp) - Comprehensive efficiency analysis with conductor, dielectric, and mismatch losses
  - Ohmic loss calculations with skin effect modeling
  - Dielectric loss analysis for various materials
  - Ground loss effects for finite ground planes
  - Mismatch efficiency from VSWR measurements
  - Material comparison (copper, aluminum, steel, etc.)
  - Aperture efficiency for reflector antennas

#### Antenna Types
- [Dipole](src/antenna/dipole_resonant_frequency.hpp) - Half-wave dipole calculations
- [Helical Antenna](src/antenna/helical_antenna.hpp) - Helical antenna gain and geometry
- [Microstrip Patch](src/antenna/microstrip_patch_antenna.hpp) - Microstrip antenna design
- [Near Field Distance](src/antenna/near_field_distance.hpp) - Near/far field boundaries

### Propagation Models

#### Basic Propagation
- [Free Space Path Loss](src/propagation/free_space_path_loss.hpp) - FSPL calculations
- [Friis Transmission](src/propagation/friis_transmission_equation.hpp) - Link budget basics
- [Friis Polarization](src/propagation/friis_with_polarization_mismatch.hpp) - Polarization losses
- [Ground Wave Basic](src/propagation/ground_wave_propagation.hpp) - Basic ground wave

#### HF Ionospheric Propagation
- [Ionospheric Ray Tracing](src/propagation/ionospheric_ray_tracing.hpp) - **Complete HF propagation prediction system**
  - Hamilton's equations solver for accurate ray trajectories
  - MUF (Maximum Usable Frequency) calculations for any path
  - Multi-hop analysis for long-distance HF communication
  - Chapman layer ionospheric profile modeling
  - Real-time ionospheric condition integration
  - Ham radio band propagation prediction (80m-10m)
- [Ionospheric Refractive Index](src/propagation/ionospheric_refractive_index.hpp) - N² calculations with collision frequency
- [Ray Trajectory](src/propagation/ray_trajectory.hpp) - Hamilton's equations for ray tracing
- [Group Phase Path](src/propagation/group_phase_path.hpp) - Group and phase path integration
- [Ionospheric Scintillation](src/propagation/ionospheric_scintillation.hpp) - Rytov variance and S4 scintillation indices

#### VHF/UHF Tropospheric Propagation
- [Tropospheric Scatter](src/propagation/tropospheric_scatter.hpp) - **Professional beyond-horizon modeling**
  - Volume scattering theory with Cn² turbulence parameters
  - Path loss analysis with detailed component breakdown
  - Availability statistics with log-normal fading models
  - Frequency and space diversity calculations
  - Optimal frequency selection for required availability
  - Real-world ham radio scenarios (VHF/UHF weak signal, contest operation)
  - Geographic coordinate path calculation
  - Atmospheric condition modeling (seasonal variations)
- [Tropospheric Ducting](src/propagation/tropospheric_ducting.hpp) - Ducting and evaporation ducts
- [Atmospheric Attenuation](src/propagation/atmospheric_attenuation.hpp) - Gas, rain, cloud losses

#### Ground Wave and Terrain
- [Ground Wave](src/propagation/ground_wave.hpp) - Sommerfeld ground wave solution
- [Terrain Diffraction](src/propagation/terrain_diffraction.hpp) - Knife-edge and multiple obstacles
- [Hata Model](src/propagation/hata_model.hpp) - Urban/suburban/rural path loss

### Feedline Models
- [Smith Chart](src/feedline/smith_chart.hpp)
  - Impedance to reflection coefficient conversions
  - Transmission line impedance transformations
  - L-network design for impedance matching
  - Stub matching network calculations
  - Quarter-wave transformer design
  - VSWR circle and constant resistance/reactance circles
  - Batch processing for multiple impedances
- [Coaxial Cable](src/feedline/coaxial_cable_impedance.hpp) - Coax characteristic impedance
- [VSWR](src/feedline/vswr_calculations.hpp) - Standing wave ratio calculations
- [Twin Lead](src/feedline/twin_lead_impedance.hpp) - Parallel transmission lines
- [Waveguide](src/feedline/waveguide_rectangular.hpp) - Rectangular waveguide parameters

### RF System Parameters
- [S-Parameters](src/rf_systems/s_parameter_conversions.hpp) - Network parameter conversions
- [Intermodulation](src/rf_systems/intermodulation_distortion.hpp) - IP3 and IMD calculations
- [Phase Noise](src/rf_systems/phase_noise.hpp) - Phase noise and jitter analysis


## Quick Start

```bash
# Clone the repository
git clone https://github.com/rowingdude/rvl.git
cd rvl

# Build and test individual components
g++ -std=c++20 -O2 -I. test_antenna_design_foundations.cpp src/antenna/*.cpp src/feedline/*.cpp -o test_antenna
g++ -std=c++20 -O2 -I. test_ionospheric_ray_tracing.cpp src/propagation/*.cpp -o test_ionospheric
g++ -std=c++20 -O2 -I. test_tropospheric_scatter.cpp src/propagation/*.cpp -o test_tropospheric

# Run tests
./test_antenna
./test_ionospheric
./test_tropospheric
```

### Example Usage

To build a prototype 4el Yagi, you could:

```cpp
#include "src/antenna/array_factor.hpp"
#include "src/propagation/ionospheric_ray_tracing.hpp"
#include "src/propagation/tropospheric_scatter.hpp"

using namespace rvl::antenna;
using namespace rvl::propagation;

array_factor_d::linear_array_params yagi;
yagi.num_elements = 4;
yagi.element_spacing_m = 0.35;  // 0.35λ spacing
yagi.frequency_hz = 144.2e6;    // 2m band

auto directivity = array_factor_d::calculate_array_directivity_db(yagi);
std::cout << "Yagi directivity: " << directivity << " dBi" << std::endl;

auto ionosphere = ionospheric_ray_tracing_d::create_chapman_profile(12e6, 320.0, 45.0);
auto muf_analysis = ionospheric_ray_tracing_d::calculate_muf(
    tx_position, rx_position, ionosphere);
std::cout << "MUF: " << muf_analysis.maximum_usable_frequency_hz / 1e6 << " MHz" << std::endl;

auto geometry = tropospheric_scatter_d::create_path_geometry(
    40.7, -74.0, 100.0,  // NYC coordinates, 100m height
    42.4, -71.1, 100.0   // Boston coordinates
);
auto atmosphere = tropospheric_scatter_d::create_standard_atmosphere();
auto scatter_analysis = tropospheric_scatter_d::calculate_enhanced_scatter_analysis(
    432e6, geometry, atmosphere, 12.0, 12.0);  // 432 MHz with 12 dBi antennas
```


## Dependencies

**Required:**
- C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- Standard library support for `<complex>`, `<vector>`, `<cmath>`

**Optional:**
- CMake 3.20+ (for build system)
- OpenMP (for parallelization)
- Intel TBB (for advanced threading)
- Google Test (for unit tests)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

## License
MIT
