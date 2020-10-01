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
    real::MOI.ConstraintIndex{F, MOI.Zeros}
    imag::MOI.ConstraintIndex{F, MOI.Zeros}
end
function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SplitZeroBridge{T, F, G}}, model::MOI.ModelLike,
    f::G,
    set::MOI.Zeros
) where {T, F, G}
    real_con = MOI.add_constraint(model, operate_coefficients(real, T, f), set)
    imag_con = MOI.add_constraint(model, operate_coefficients(imag, T, f), set)
    return SplitZeroBridge{T, F, G}(real_con, imag_con)
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
    return 2
end
function MOI.get(bridge::SplitZeroBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.Zeros}) where {T, F}
    return [bridge.real, bridge.imag]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SplitZeroBridge)
    MOI.delete(model, bridge.imag)
    MOI.delete(model, bridge.real)
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
    return MOI.get(model, attr, bridge.real) + im * MOI.get(model, attr, bridge.imag)
end
function MOI.set(model::MOI.ModelLike, attr::Union{MOI.ConstraintPrimalStart, MOI.ConstraintDualStart},
                 bridge::SplitZeroBridge, value)
    MOI.set(model, attr, bridge.real, map(real, value))
    MOI.set(model, attr, bridge.imag, map(imag, value))
end
