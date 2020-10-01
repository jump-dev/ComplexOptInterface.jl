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

function zero_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.Nonnegatives(2))
    fx = MOI.SingleVariable.(x)
    func = (1.0 + 2.0im) * fx[1] + (1.0 - 1.0im) * fx[2] + (-1.0 + -1.0im)
    c = MOI.add_constraint(
        optimizer,
        MOI.Utilities.operate(vcat, Complex{Float64}, func),
        MOI.Zeros(1)
    )
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ [2/3, 1/3] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ [2/3, 1/3] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(2) atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c) ≈ [0.0 + 0.0im] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), c) ≈ [0.0 + 0.0im] atol=atol rtol=rtol
end

import CSDP
@testset "CSDP" begin
    config = MOI.Test.TestConfig(atol=1e-4, rtol=1e-4)
    bridged = MOI.instantiate(CSDP.Optimizer, with_bridge_type=Float64)
    MOI.Bridges.add_bridge(bridged, COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64})
    projection_test(bridged, config)
    bridged = MOI.Bridges.LazyBridgeOptimizer(MOI.Utilities.CachingOptimizer(MOI.Utilities.Model{Float64}(), CSDP.Optimizer()))
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.ScalarizeBridge{Float64})
    MOI.Bridges.add_bridge(bridged, COI.Bridges.Constraint.SplitZeroBridge{Float64})
    zero_test(bridged, config)
end
