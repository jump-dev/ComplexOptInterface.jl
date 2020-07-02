module ComplexOptInterface

import MutableArithmetics
const MA = MutableArithmetics

import MathOptInterface
const MOI = MathOptInterface

struct HermitianPositiveSemidefiniteConeTriangle <: MOI.AbstractVectorSet
    side_dimension::Int
end
function MOI.dimension(set::HermitianPositiveSemidefiniteConeTriangle)
    return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)) +
        MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension - 1))
end

function MOI.Utilities.set_dot(x::AbstractVector, y::AbstractVector,
                               set::HermitianPositiveSemidefiniteConeTriangle)
    sym = MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)
    result = MOI.Utilities.set_dot(x, y, sym)
    for k in (MOI.dimension(sym) + 1):MOI.dimension(set)
        result = MA.add_mul!(result, 2, x[k], y[k])
    end
    return result
end

include("Bridges/Bridges.jl")

end # module
