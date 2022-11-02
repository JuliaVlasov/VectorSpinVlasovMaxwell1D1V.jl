# VectorSpinVlasovMaxwell1D1V.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaVlasov.github.io/VectorSpinVlasovMaxwell1D1V.jl/dev/)
[![Build Status](https://github.com/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl)

```bash
git clone https://github.com/JuliaVlasov/VectorSpinVlasovMaxwell1D1V.jl.git
cd VectorSpinVlasovMaxwell1D1V.jl
julia -O3 --check-bounds=no --project examples/inplace.jl
```
This first implementation is translated from a MATLAB version. 

There is two implementations of the same algorithm. If you want to test the other one

```bash
julia -O3 --check-bounds=no --project examples/operators.jl
```

In this version, we use Julia Structs to reduce memory print and improve the readibility.
