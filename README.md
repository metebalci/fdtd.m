# fdtd.m

This repo contains FDTD EM Solver Implementations I created using Octave while I was following the lectures of [Dr. Raymond Rumpf](https://empossible.net/).

The implementations are as follows:

# fdtd_1d

1D FDTD implementation. This is a complete (but without dispersion correction and lossy material support) implementation in 1-dimension:

- Implements a perfect boundary condition (no reflections)
- Uses a TF/SF source at z=2
- Measures reflectance at z=1
- Measures transmittance at the last z
- There are 10 (spacer_region variable in fdtd_init function) spacer regions (vacuum) before and after the actual structure given.

During the simulation, a plot showing Ey and Hx waves are updated in real-time, the layers of the structure are also shown in gray colors in this plot. At the end of the simulation, fourier (frequency) analysis results are also plotted.

There is also a simple 3D plot functionality but it is commented out.

Below are plots while and at the end of running fdtd_demo3:

![Waves](fdtd_1d_demo3_t.png)

![Fourier Analysis](fdtd_1d_demo3_f.png)

## Differences from the lecture

There are two optimizations in the implementation comparing to the lectures:

- The inner E and H update iterations (in space) are replaced by vectoral implementations. This is obviously faster.

- Fourier kernels (K) are pre-computed but also multiplied at each time iteration (Kc) in order to not run power operation (^T) as in the lectures.

- Fourier calculation iterations (in frequency) are also replaced by vectoral implementations.

- The frequency domain plots are given as s-parameters, in dB, rather than plain power ratios of reflection/transmittance to source.

## Demos

Three demos (also based on the examples from the lectures) are provided:

- fdtd_demo1: single region
- fdtd_demo2: anti-reflection layers
- fdtd_demo3: bragg gratings

## How to use

- Run fdtd_1d first.
- Then either run fdtd_demo functions directly or fdtd_run function yourself.

fdtd_run parameters are:

- layers: Nx3 matrix for N layers, each row contains one layer with [thickness dielectric_constant relative_permeability]
- length_unit: the unit used in the layers structure, e.g. if you are using mm then provide 1e-3 here
- f_max: maximum frequency to simulate, used for preparing the source and fourier analysis
