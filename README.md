# ComplexOptInterface

**This package is still in early development, feedback is welcome!**

| **Build Status** | **Social** |
|:----------------:|:----------:|
| [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] |
| [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/en/a/af/Discourse_logo.png" width="64">][discourse-url] |

An extension of MathOptInterface to complex sets.

## How to use

Add to the JuMP model the bridges of this package with `add_all_bridges` and then you can create complex equality constraints and complex hermitian matrices:
```julia
using JuMP
import ComplexOptInterface
const COI = ComplexOptInterface

model = Model()
COI.add_all_bridges(model)
@variable(model, x in COI.ComplexPlane(), start = 5 + 6im, lower_bound = 1 + 2im, upper_bound = 3 + 4im)
@variable(model, Q[1:2, 1:2] in COI.HermitianPSDCone())
@constraint(model, (1 + 2im) * x + Q[1, 1] + Q[2, 2] * im == 1 + 2im)
```

## Design considerations

There are two types of complex sets:
1) Some sets have complex entries, i.e. `MOI.EqualTo{Complex{Float64}}`,
2) Some sets have real entries even if they model a complex mathematical set, i.e.
   `HermitianPositiveSemidefiniteConeTriangle` where the imaginary part of the
   off-diagonal entries are stored in separate output indices.
   For instance, `[x y+zim; y-zim w]` is vectorized as
   `[x, y, w, z]`.

Sets with complex entries **should not** be used with `VariableIndex` or `VectorOfVariables` constraints
as the variables in MathOptInterface are assumed to be real.
Indeed, for instance doing
```julia
x, cx = MOI.add_constrained_variable(model, MOI.EqualTo(1 + 2im))
```
fallbacks to
```julia
x = MOI.add_variable(model)
cx = MOI.add_constraint(model, x, MOI.EqualTo(1 + 2im))
```
In the first line, the solvers create a real variable.
Moreover, in the bridge from `MOI.ScalarAffineFunction{Complex{T}}`-in-`EqualTo{Complex{T}}`
to `ScalarAffineFunction{T}`-in-`EqualTo{T}`, we assume that the variables are real
to split the equality in real and imaginary part.
So in conclusion, complex variables may violate assumptions in MOI so use it at your own risks.
Note that even if you add a `MOI.ScalarAffineFunction{Complex{T}}`-in-`EqualTo{Complex{T}}` constraint,
it may be bridged by a slack bridge into a variable constrained in `EqualTo{Complex{T}}`.
Luckily, this bridge might not be selected but remain causious with sets with complex entries.

[build-img]: https://travis-ci.com/jump-dev/ComplexOptInterface.jl.svg?branch=master
[build-url]: https://travis-ci.com/jump-dev/ComplexOptInterface.jl
[codecov-img]: http://codecov.io/github/jump-dev/ComplexOptInterface.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/jump-dev/ComplexOptInterface.jl?branch=master

[gitter-url]: https://gitter.im/JuliaOpt/JuMP-dev?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/JuMP-dev.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt
