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
    x11 = x[1:3]
    x12 = x[4]
    t = MOI.add_variable(optimizer)
    MOI.add_constraint(
        optimizer,
        MOI.Utilities.operate(
            vcat,
            Float64,
            t,
            x[1] - 1.0,
            √2 * (x[2] + 1.0),
            x[3] + 1.0,
            √2 * (x[4] - 1.0),
        ),
        MOI.SecondOrderCone(5),
    )
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(optimizer, MOI.ObjectiveFunction{typeof(t)}(), t)
    MOI.optimize!(optimizer)
    primal = [(1 + √3) / 2, -1 / 2, (-1 + √3) / 2, 1 / 2]
    dual = [(3 - √3) / 6, 0.2886751198, (3 + √3) / 6, -0.2886751197]
    @test MOI.Utilities.set_dot(primal, dual, set) ≈ 0 atol = atol rtol = rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ primal atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ primal atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ dual atol = atol rtol =
        rtol
end

function hermitian_psd_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    set = COI.HermitianPositiveSemidefiniteConeTriangle(3)
    x, cx = MOI.add_constrained_variables(optimizer, set)
    MOI.add_constraints(
        optimizer,
        x,
        MOI.EqualTo.([1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0]),
    )
    MOI.optimize!(optimizer)
    primal = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0]
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ primal atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ primal atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(9) atol = atol rtol =
        rtol
end

function equalto_1_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.Nonnegatives(2))
    func = (1.0 + 2.0im) * x[1] + (1.0 - 1.0im) * x[2]
    c = MOI.add_constraint(optimizer, func, MOI.EqualTo(1.0 + 1.0im))
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ [2 / 3, 1 / 3] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ [2 / 3, 1 / 3] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(2) atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c) ≈ 1.0 + 1.0im atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), c) ≈ 0.0 + 0.0im atol = atol rtol =
        rtol
end

function equalto_2_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.Nonnegatives(1))
    func = (1.0 + 0.0im) * x[1] + 1.0im * x[1] - (1.0 + 0.0im) * x[1]
    c = MOI.add_constraint(optimizer, func, MOI.EqualTo(2.0im))
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ [2.0] atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ [2.0] atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(1) atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c) ≈ 2.0im atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), c) ≈ 0.0 + 0.0im atol = atol rtol =
        rtol
end

function zero_1_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.Nonnegatives(2))
    func = (1.0 + 2.0im) * x[1] + (1.0 - 1.0im) * x[2] + (-1.0 + -1.0im)
    c = MOI.add_constraint(
        optimizer,
        MOI.Utilities.operate(vcat, Complex{Float64}, func),
        MOI.Zeros(1),
    )
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ [2 / 3, 1 / 3] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ [2 / 3, 1 / 3] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(2) atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c) ≈ [0.0 + 0.0im] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), c) ≈ [0.0 + 0.0im] atol =
        atol rtol = rtol
end

function zero_2_test(optimizer, config)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.Nonnegatives(1))
    func = (1.0 + 0.0im) * x[1] + 1.0im * x[1] - 2.0im - (1.0 + 0.0im) * x[1]
    c = MOI.add_constraint(
        optimizer,
        MOI.Utilities.operate(vcat, Complex{Float64}, func),
        MOI.Zeros(1),
    )
    MOI.optimize!(optimizer)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ [2.0] atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ [2.0] atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ zeros(1) atol = atol rtol =
        rtol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), c) ≈ [0.0 + 0.0im] atol =
        atol rtol = rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), c) ≈ [0.0 + 0.0im] atol =
        atol rtol = rtol
end

import CSDP
@testset "CSDP" begin
    config = MOI.Test.Config(atol = 1e-4, rtol = 1e-4)
    bridged = MOI.instantiate(CSDP.Optimizer, with_bridge_type = Float64)
    MOI.Bridges.add_bridge(
        bridged,
        COI.Bridges.Variable.HermitianToSymmetricPSDBridge{Float64},
    )
    projection_test(bridged, config)
    hermitian_psd_test(bridged, config)

    bridged = MOI.Bridges.LazyBridgeOptimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.Model{Float64}(),
            CSDP.Optimizer(),
        ),
    )
    MOI.Bridges.add_bridge(
        bridged,
        COI.Bridges.Constraint.SplitEqualToBridge{Float64},
    )
    equalto_1_test(bridged, config)
    cis = MOI.get(
        bridged.model,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    @test length(cis) == 2
    equalto_2_test(bridged, config)
    cis = MOI.get(
        bridged.model,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    @test length(cis) == 1

    bridged = MOI.Bridges.LazyBridgeOptimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.Model{Float64}(),
            CSDP.Optimizer(),
        ),
    )
    MOI.Bridges.add_bridge(
        bridged,
        MOI.Bridges.Constraint.ScalarizeBridge{Float64},
    )
    MOI.Bridges.add_bridge(
        bridged,
        COI.Bridges.Constraint.SplitZeroBridge{Float64},
    )
    zero_1_test(bridged, config)
    cis = MOI.get(
        bridged.model,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    @test length(cis) == 2
    zero_2_test(bridged, config)
    cis = MOI.get(
        bridged.model,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
    @test length(cis) == 1
end
