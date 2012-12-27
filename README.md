StencilProbe: A Microbenchmark for Stencil Applications
=====
The Stencil Probe is a small, self-contained serial microbenchmark that we developed as a tool to explore the behavior of grid-based computations. As such it is suitable for experimentation on architectures in varying states of implementation -- from production CPUs to cycle-accurate simulators. By modifying the operations in the inner loop of the benchmark, the Stencil Probe can effectively mimic the kernels of applications that use stencils on regular grids. In this way, we can easily simulate the memory access patterns and performance of large applications, as well as use the Stencil Probe as a testbed for potential optimizations, without having to port or modify the entire application.

See a longer description at http://people.csail.mit.edu/skamil/projects/stencilprobe
