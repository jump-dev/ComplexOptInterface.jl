struct SplitEqualToBridge{
    T,
    F<:MOI.Utilities.TypedScalarLike{T},
    G<:MOI.Utilities.TypedScalarLike{Complex{T}},
} <: MOI.Bridges.Constraint.AbstractBridge
    real_constraint::Union{Nothing,MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
    imag_constraint::Union{Nothing,MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
end
function _add_constraint_if_nonzero(model, func, set)
    if iszero(func) && iszero(MOI.constant(set))
        return nothing
    else
        return MOI.add_constraint(model, func, set)
    end
end
function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SplitEqualToBridge{T,F,G}},
    model::MOI.ModelLike,
    func::G,
    set::MOI.EqualTo,
) where {T,F,G}
    real_func = real(func)
    imag_func = MOI.Utilities.operate(imag, T, func)
    real_set = MOI.EqualTo(real(MOI.constant(set)))
    imag_set = MOI.EqualTo(imag(MOI.constant(set)))
    real_constraint = _add_constraint_if_nonzero(model, real_func, real_set)
    imag_constraint = _add_constraint_if_nonzero(model, imag_func, imag_set)
    return SplitEqualToBridge{T,F,G}(real_constraint, imag_constraint)
end

# We don't support `MOI.VariableIndex` as it would be a self-loop in the bridge graph
function MOI.supports_constraint(
    ::Type{SplitEqualToBridge{T}},
    ::Type{<:MOI.Utilities.TypedLike{Complex{T}}},
    ::Type{<:MOI.EqualTo},
) where {T}
    return true
end
MOIB.added_constrained_variable_types(::Type{<:SplitEqualToBridge}) = Tuple{DataType}[]
function MOIB.added_constraint_types(::Type{SplitEqualToBridge{T,F,G}}) where {T,F,G}
    return Tuple{DataType,DataType}[(F, MOI.EqualTo{T})]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SplitEqualToBridge{T}},
    G::Type{<:MOI.Utilities.TypedLike},
    ::Type{<:MOI.EqualTo},
) where {T}
    F = MA.promote_operation(imag, G)
    return SplitEqualToBridge{T,F,G}
end

# Attributes, Bridge acting as a model
function MOI.get(
    bridge::SplitEqualToBridge{T,F},
    ::MOI.NumberOfConstraints{F,MOI.EqualTo{T}},
) where {T,F}
    return !isnothing(bridge.real_constraint) + !isnothing(bridge.imag_constraint)
end
function MOI.get(
    bridge::SplitEqualToBridge{T,F},
    ::MOI.ListOfConstraintIndices{F,MOI.EqualTo{T}},
) where {T,F}
    list = MOI.ConstraintIndex{F,MOI.EqualTo{T}}[]
    if !isnothing(bridge.real_constraint)
        push!(list, bridge.real_constraint)
    end
    if !isnothing(bridge.imag_constraint)
        push!(list, bridge.imag_constraint)
    end
    return list
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SplitEqualToBridge)
    if !isnothing(bridge.real_constraint)
        MOI.delete(model, bridge.real_constraint)
    end
    if !isnothing(bridge.imag_constraint)
        MOI.delete(model, bridge.imag_constraint)
    end
end

# Attributes, Bridge acting as a constraint
function MOI.supports(
    ::MOI.ModelLike,
    ::Union{MOI.ConstraintPrimalStart,MOI.ConstraintDualStart},
    ::Type{<:SplitEqualToBridge},
)
    return true
end
function MOI.get(
    model::MOI.ModelLike,
    attr::Union{
        MOI.ConstraintPrimal,
        MOI.ConstraintPrimalStart,
        MOI.ConstraintDual,
        MOI.ConstraintDualStart,
    },
    bridge::SplitEqualToBridge{T},
) where {T}
    if isnothing(bridge.real_constraint)
        real_value = zero(T)
    else
        real_value = MOI.get(model, attr, bridge.real_constraint)
    end
    if isnothing(bridge.imag_constraint)
        imag_value = zero(T)
    else
        imag_value = MOI.get(model, attr, bridge.imag_constraint)
    end
    return real_value + imag_value * im
end
function MOI.set(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintPrimalStart,MOI.ConstraintDualStart},
    bridge::SplitEqualToBridge{T},
    value,
) where {T}
    if !isnothing(bridge.real_constraint)
        MOI.set(model, attr, bridge.real_constraint, real(value))
    end
    if !isnothing(bridge.imag_constraint)
        MOI.set(model, attr, bridge.real_constraint, imag(value))
    end
    return
end
