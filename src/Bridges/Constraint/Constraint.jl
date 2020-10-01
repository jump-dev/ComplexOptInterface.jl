module Constraint

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import ComplexOptInterface
const COI = ComplexOptInterface

include("split_zero.jl")

end
