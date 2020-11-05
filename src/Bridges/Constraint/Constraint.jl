module Constraint

import MutableArithmetics
const MA = MutableArithmetics

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import ComplexOptInterface
const COI = ComplexOptInterface

include("split_zero.jl")

end
