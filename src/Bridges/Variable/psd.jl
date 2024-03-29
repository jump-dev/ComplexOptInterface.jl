const EQ{T} = MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}}

"""
Let `H = HermitianToSymmetricPSDBridge(n)`, `S = PositiveSemidefiniteConeTriangle(2n)`
and `P = S ∩ con_11_22 ∩ con_12_21 ∩ con12diag`.
Suppose for simplicity that the elements of P are ordered as:
```
\\ 1 |\\ 2
 \\  | 3
  \\ |4_\\
   \\  5
    \\
     \\
```
We have `P = A * H` where
```
    [I  0]
    [0  I]
A = [0  0]
    [0 -I]
    [I  0]
```
Therefore, `H* = A* * P*` where
```
     [I 0 0  0 I]
A* = [0 I 0 -I 0]
```
Moreover, as `(S ∩ T)* = S* + T*` for cones `S` and `T`, we have
```
P* = S* + con_11_22* + con_12_21* + con12diag*
```
the dual vector of `P*` is the dual vector of `S*` for which we add in the corresponding
entries the dual of the three constraints, multiplied by the coefficients for the `EqualTo` constraints.
Note that these contribution cancel out when we multiply them by `A*`:
A* * (S* + con_11_22* + con_12_21* + con12diag*) = A* * S*
so we can just ignore them.
"""
struct HermitianToSymmetricPSDBridge{T} <: MOIB.Variable.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    psd_constraint::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    }
    con_11_22::Vector{EQ{T}}
    con12diag::Vector{EQ{T}}
    con_12_21::Vector{EQ{T}}
end

function MOIB.Variable.bridge_constrained_variable(
    ::Type{HermitianToSymmetricPSDBridge{T}},
    model::MOI.ModelLike,
    set::COI.HermitianPositiveSemidefiniteConeTriangle,
) where {T}
    n = set.side_dimension
    variables, psd_constraint = MOI.add_constrained_variables(
        model,
        MOI.PositiveSemidefiniteConeTriangle(2n),
    )

    k11 = 0
    k12 = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    k21 = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(2n)) + 1
    k22 = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    X11() = variables[k11]
    X12() = variables[k12]
    function X21(i, j)
        I = j
        J = n + i
        k21 = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(J - 1)) + I
        return variables[k21]
    end
    X22() = variables[k22]
    con_11_22 = EQ{T}[]
    con12diag = EQ{T}[]
    con_12_21 = EQ{T}[]
    for j in 1:n
        k22 += n
        for i in 1:j
            k11 += 1
            k12 += 1
            k22 += 1
            push!(
                con_11_22,
                MOI.add_constraint(
                    model,
                    MOI.Utilities.operate(-, T, X11(), X22()),
                    MOI.EqualTo(zero(T)),
                ),
            )
            if i == j
                push!(
                    con12diag,
                    MOI.add_constraint(
                        model,
                        convert(MOI.ScalarAffineFunction{T}, X12()),
                        MOI.EqualTo(zero(T)),
                    ),
                )
            else
                push!(
                    con_12_21,
                    MOI.add_constraint(
                        model,
                        MOI.Utilities.operate(+, T, X21(i, j), X12()),
                        MOI.EqualTo(zero(T)),
                    ),
                )
            end
        end
        k12 += n
    end

    return HermitianToSymmetricPSDBridge(
        variables,
        psd_constraint,
        con_11_22,
        con12diag,
        con_12_21,
    )
end

function MOIB.Variable.supports_constrained_variable(
    ::Type{<:HermitianToSymmetricPSDBridge},
    ::Type{COI.HermitianPositiveSemidefiniteConeTriangle},
)
    return true
end
function MOIB.added_constrained_variable_types(
    ::Type{<:HermitianToSymmetricPSDBridge},
)
    return [(MOI.PositiveSemidefiniteConeTriangle,)]
end
function MOIB.added_constraint_types(
    ::Type{HermitianToSymmetricPSDBridge{T}},
) where {T}
    return [(MOI.ScalarAffineFunction{T}, MOI.EqualTo{T})]
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::HermitianToSymmetricPSDBridge, ::MOI.NumberOfVariables)
    return length(bridge.variables)
end
function MOI.get(
    bridge::HermitianToSymmetricPSDBridge,
    ::MOI.ListOfVariableIndices,
)
    return bridge.variables
end
function MOI.get(
    bridge::HermitianToSymmetricPSDBridge,
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    },
)
    return 1
end
function MOI.get(
    bridge::HermitianToSymmetricPSDBridge,
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    },
)
    return [bridge.psd_constraint]
end
function MOI.get(
    bridge::HermitianToSymmetricPSDBridge{T},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
) where {T}
    return length(bridge.con_11_22) +
           length(bridge.con12diag) +
           length(bridge.con_12_21)
end
function MOI.get(
    bridge::HermitianToSymmetricPSDBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
) where {T}
    return [bridge.con_11_22; bridge.con12diag; bridge.con_12_21]
end

# References
function MOI.delete(model::MOI.ModelLike, bridge::HermitianToSymmetricPSDBridge)
    for ci in bridge.con_11_22
        MOI.delete(model, ci)
    end
    for ci in bridge.con12diag
        MOI.delete(model, ci)
    end
    for ci in bridge.con_12_21
        MOI.delete(model, ci)
    end
    return MOI.delete(model, bridge.variables)
end

# Attributes, Bridge acting as a constraint

function MOI.get(
    model::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::HermitianToSymmetricPSDBridge,
)
    return COI.HermitianPositiveSemidefiniteConeTriangle(
        length(bridge.con12diag),
    )
end

function _matrix_indices(k)
    # If `k` is a diagonal index, `s(k)` is odd and 1 + 8k is a perfect square.
    n = 1 + 8k
    s = isqrt(n)
    if s^2 == n
        j = div(s, 2)
    else
        # Otherwise, if it is after the diagonal index `k` but before the diagonal
        # index `k'` with `s(k') = s(k) + 2`, we have `s(k) <= s < s(k) + 2`.
        # By shifting by `+1` before `div`, we make sure to have the right column.
        j = div(s + 1, 2)
    end
    i = k - MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(j - 1))
    return i, j
end

function _variable_map(idx::MOIB.IndexInVector, n)
    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    if idx.value <= N
        return idx.value
    else
        i, j = _matrix_indices(idx.value - N)
        return N +
               j * n +
               MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(j)) +
               i
    end
end
function _variable(bridge::HermitianToSymmetricPSDBridge, i::MOIB.IndexInVector)
    return bridge.variables[_variable_map(i, length(bridge.con12diag))]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintPrimal,
    bridge::HermitianToSymmetricPSDBridge{T},
) where {T}
    values = MOI.get(model, attr, bridge.psd_constraint)
    M = MOI.dimension(MOI.get(model, MOI.ConstraintSet(), bridge))
    n = length(bridge.con12diag)
    return [values[_variable_map(MOIB.IndexInVector(i), n)] for i in 1:M]
end

# See docstring of bridge for why we ignore the dual of the constraints
# `con_11_22`, `con_12_21` and `con12diag`.
function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::HermitianToSymmetricPSDBridge{T},
) where {T}
    dual = MOI.get(model, attr, bridge.psd_constraint)
    M = MOI.dimension(MOI.get(model, MOI.ConstraintSet(), bridge))
    mapped = zeros(T, M)
    n = length(bridge.con12diag)
    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    k11 = 0
    k12 = N
    k21 = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(2n)) + 1
    k22 = N
    k = 0
    for j in 1:n
        k21 -= n + 1 - j
        k22 += n
        for i in 1:j
            k11 += 1
            k12 += 1
            k21 -= 1
            k22 += 1
            mapped[k11] += dual[k11]
            mapped[k11] += dual[k22]
            if i != j
                k += 1
                mapped[N+k] += dual[k12]
                mapped[N+k] -= dual[k21]
            end
        end
        k12 += n
        k21 -= n - j
    end
    return mapped
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::HermitianToSymmetricPSDBridge{T},
    i::MOIB.IndexInVector,
) where {T}
    return value = MOI.get(model, attr, _variable(bridge, i))
end

function MOIB.bridged_function(
    bridge::HermitianToSymmetricPSDBridge{T},
    i::MOIB.IndexInVector,
) where {T}
    func = _variable(bridge, i)
    return convert(MOI.ScalarAffineFunction{T}, func)
end
function MOIB.Variable.unbridged_map(
    bridge::HermitianToSymmetricPSDBridge{T},
    vi::MOI.VariableIndex,
    i::MOIB.IndexInVector,
) where {T}
    func = convert(MOI.ScalarAffineFunction{T}, vi)
    return (_variable(bridge, i) => func,)
end
