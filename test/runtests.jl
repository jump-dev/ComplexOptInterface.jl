using Test

import MathOptInterface
const MOI = MathOptInterface

import ComplexOptInterface
const COI = ComplexOptInterface

"""
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
function test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, COI.HermitianPositiveSemidefiniteConeTriangle(2))
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
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ primal atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ primal atol=atol rtol=rtol
    dual = MOI.get(optimizer, MOI.ConstraintDual(), cx)
    @test dual[1] ≈ (3 - √3) / 6 atol=atol rtol=rtol
    @test dual[2] ≈ 0.2886751197468812 atol=atol rtol=rtol
    @test dual[3] ≈ (3 + √3) / 6 atol=atol rtol=rtol
    @test dual[4] ≈ 0.0 atol=atol rtol=rtol
end

import CSDP
@testset "CSDP" begin
    bridged = MOI.instantiate(CSDP.Optimizer, with_bridge_type=Float64)
    MOI.Bridges.add_bridge(bridged, COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64})
    test(bridged, MOI.Test.TestConfig(atol=1e-4, rtol=1e-4))
end
