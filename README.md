# XCUMP

CUDA-Compatible C/C++ Multiprecision Library in VERY early development

the goal of XCUMP is to create a well-documented, easy-to-use, CUDA-Compatible Multiprecision Library that significantly outperforms every other CPU-based Multiprecision Library on most nvidia GPUs.

there will be a CPU-part of the library with arithmetic algorithms optimized for the CPU, and a CUDA-part. The conversion
between a object on the CPU and one on the GPU will be very fast. due to limitations on the GPU, it can only use 32-bit
limb size. This is not an issue on the CPU, but 32-bit limbs will still be used to make conversion to GPU O(1). because of
this, the CPU part of the library will be at least ~50% slower than other CPU-based MP libraries.

this library will be in a single header. It will support the use of templates for imortant parameters known at compile time.
there will be a C binding (without templates and such), and a C++ dll for people who dont care about templates and want faster
compile time.

it will include a Schönhagen-Strassen multiplcation algorithm. This is by far the most complicated part of the library
and the most complicated algorithm i have ever written. The most difficult part of SS (NTT and CRT) is already pretty much finished. Optimization is another thing though.

current todo:

todo:

- basic schoolbook addition: done
- basic schoolbook subtraction: -
- limb division: done
- schoolbook multiplicatio: done
- NTTCRT multiplication: almost done
- Schönhagen-Strassen Multiplication: in progress
- faster division algorithms
- optimize NTT
- optimize CRT
