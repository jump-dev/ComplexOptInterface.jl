module ComplexOptInterface

import MathOptInterface
const MOI = MathOptInterface

struct PositiveSemidefiniteConeTriangle <: MOI.AbstractSymmetricMatrixSetTriangle
    side_dimension::Int
end

struct HermitianPositiveSemidefiniteConeTriangle <: MOI.AbstractVectorSet
    side_dimension::Int
end
function MOI.dimension(set::HermitianPositiveSemidefiniteConeTriangle)
    return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)) +
        MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension - 1))
end


include("Bridges/Bridges.jl")

end # module
