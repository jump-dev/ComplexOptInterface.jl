function operate_coefficient(f, T::Type, term::MOI.ScalarAffineTerm)
    MOI.ScalarAffineTerm(f(term.coefficient), term.variable_index)
end
function operate_coefficient(f, T::Type, term::MOI.VectorAffineTerm)
    return MOI.VectorAffineTerm(term.output_index, operate_coefficient(f, T, term.scalar_term))
end
function operate_coefficients(f, T::Type, func::MOI.VectorAffineFunction)
    return MOI.VectorAffineFunction(
        [operate_coefficient(f, T, term) for term in func.terms],
        map(f, func.constants)
    )
end

similar_type(::Type{<:MOI.VectorAffineFunction}, T::Type) = MOI.VectorAffineFunction{T}

struct SplitZeroBridge{T, F<:MOI.Utilities.TypedLike{T}, G<:MOI.Utilities.TypedLike{Complex{T}}} <: MOI.Bridges.Constraint.AbstractBridge
    dimension::Int
    constraint::MOI.ConstraintIndex{F, MOI.Zeros}
    real_indices::Vector{Int}
    imag_indices::Vector{Int}
end
function _nonzero_indices(func::MOI.AbstractVectorFunction)
    return [i for (i, scalar_func) in enumerate(MOIU.scalarize(func)) if !iszero(scalar_func)]
end
function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SplitZeroBridge{T, F, G}}, model::MOI.ModelLike,
    f::G,
    set::MOI.Zeros
) where {T, F, G}
    real_part = operate_coefficients(real, T, f)
    imag_part = operate_coefficients(imag, T, f)
    real_indices = _nonzero_indices(real_part)
    imag_indices = _nonzero_indices(imag_part)
    func = MOIU.operate(
        vcat, T,
        MOIU.eachscalar(real_part)[real_indices],
        MOIU.eachscalar(imag_part)[imag_indices]
    )
    constraint = MOI.add_constraint(model, func, MOI.Zeros(length(real_indices) + length(imag_indices)))
    return SplitZeroBridge{T, F, G}(MOI.dimension(set), constraint, real_indices, imag_indices)
end

function MOI.supports_constraint(
    ::Type{SplitZeroBridge{T}}, ::Type{<:MOI.Utilities.TypedLike{Complex{T}}},
    ::Type{MOI.Zeros}) where T
    return true
end
MOIB.added_constrained_variable_types(::Type{<:SplitZeroBridge}) = Tuple{DataType}[]
function MOIB.added_constraint_types(::Type{SplitZeroBridge{T, F, G}}) where {T, F, G}
    return Tuple{DataType, DataType}[(F, MOI.Zeros)]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SplitZeroBridge{T}}, G::Type{<:MOI.Utilities.TypedLike},
    ::Type{MOI.Zeros}) where T
    return SplitZeroBridge{T, similar_type(G, T), G}
end

# Attributes, Bridge acting as a model
function MOI.get(::SplitZeroBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.Zeros}) where {T, F}
    return 1
end
function MOI.get(bridge::SplitZeroBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.Zeros}) where {T, F}
    return [bridge.constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SplitZeroBridge)
    MOI.delete(model, bridge.constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.supports(
    ::MOI.ModelLike,
    ::Union{MOI.ConstraintPrimalStart, MOI.ConstraintDualStart},
    ::Type{<:SplitZeroBridge})

    return true
end
function MOI.get(model::MOI.ModelLike, attr::Union{MOI.ConstraintPrimal, MOI.ConstraintPrimalStart, MOI.ConstraintDual, MOI.ConstraintDualStart},
                 bridge::SplitZeroBridge)
    values = MOI.get(model, attr, bridge.constraint)
    output = zeros(Complex{eltype(values)}, bridge.dimension)
    for (i, idx) in enumerate(bridge.real_indices)
        output[idx] = values[i]
    end
    for (i, idx) in enumerate(bridge.imag_indices)
        output[idx] = values[length(bridge.real_indices) + i] * im
    end
    return output
end
function MOI.set(model::MOI.ModelLike, attr::Union{MOI.ConstraintPrimalStart, MOI.ConstraintDualStart},
                 bridge::SplitZeroBridge{T}, value) where T
    input = Vector{T}(undef, length(bridge.real_indices) + length(bridge.imag_indices))
    for (i, idx) in enumerate(bridge.real_indices)
        input[i] = real(value[idx])
    end
    for (i, idx) in enumerate(bridge.imag_indices)
        input[length(bridge.real_indices) + i] = imag(value[idx])
    end
    MOI.set(model, attr, bridge.constraint, input)
end
