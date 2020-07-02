using Test

import MathOptInterface
const MOI = MathOptInterface

import ComplexOptInterface
const COI = ComplexOptInterface

"""
    projection_test(optimizer, config)

Test computation of the projection of a hermitian matrix to the cone
of hermitian positive semidefinite matrices.

[ 1    -1+im]   [ 1-im  ?]   [√3    ]   [ 1+im  -1+√3]
[-1-im -1   ] = [-1+√3  ?] * [   -√3] * [  ?      ?  ] / (6 - 2√3)
So it's projection to the Hermitian PSD cone is
                 [1-im]  [ 1+im  1-√3]
√3 / (6 - 2√3) * [1-√3]
which is
                 [1-im]  [ 1+im  1-√3]
(1 + √3) / 4   * [1-√3]
which is
[(1+√3)/2  -1/2+im/2]
[-1/2+im/2 (-1+√3)/2]
"""
function projection_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    set = COI.HermitianPositiveSemidefiniteConeTriangle(2)
    x, cx = MOI.add_constrained_variables(optimizer, set)
    fx = MOI.SingleVariable.(x)
    x11 = fx[1:3]
    x12 = fx[4]
    t = MOI.add_variable(optimizer)
    ft = MOI.SingleVariable(t)
    MOI.add_constraint(optimizer, MOI.Utilities.operate(vcat, Float64, ft, fx[1] - 1.0, √2 * (fx[2] + 1.0), fx[3] + 1.0, √2 * (fx[4] - 1.0)),
                       MOI.SecondOrderCone(5))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(optimizer, MOI.ObjectiveFunction{typeof(ft)}(), ft)
    MOI.optimize!(optimizer)
    primal = [(1 + √3) / 2, -1/2, (-1 + √3)/2, 1/2]
    dual = [(3 - √3) / 6, 0.2886751198, (3 + √3) / 6, -0.2886751197]
    @test MOI.Utilities.set_dot(primal, dual, set) ≈ 0 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ primal atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ primal atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ dual atol=atol rtol=rtol
end

import CSDP
@testset "CSDP" begin
    bridged = MOI.instantiate(CSDP.Optimizer, with_bridge_type=Float64)
    MOI.Bridges.add_bridge(bridged, COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64})
    test(bridged, MOI.Test.TestConfig(atol=1e-4, rtol=1e-4))
end
