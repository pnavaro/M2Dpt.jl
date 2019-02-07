# M2Dpt.jl

M2Dpt: 	2-D pseudo-transient routines, part of the [M2Di software distribution](https://bitbucket.org/lraess/m2di/).

Julia version of the matlab script `TM2Dpt_GJI.m`

https://doi.org/10.1093/gji/ggy434.

- MATLAB 2018a (probably uses more than one thread)
```
219 seconds
```

- Julia vectorized code (`test` directory)

```
475.035169 seconds (70.78 M allocations: 1005.967 GiB, 14.47% gc time)
```

- Julia code with loops (`test` directory).

```
297.820268 seconds (2.81 M allocations: 407.538 MiB, 0.05% gc time)
```

- Fortran reference program

In fortran directory type `make` and run the executable program `m2dpt`:

```
227 seconds
```
