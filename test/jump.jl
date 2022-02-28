module TestJuMP

using JuMP
import GLPK
import ComplexOptInterface
const COI = ComplexOptInterface
using Test

function test_readme_example()
    model = Model()
    COI.add_all_bridges(model)
    @variable(model, x[1:2, 1:2] in COI.HermitianPSDCone())
    @test num_variables(model) == 4
    v = all_variables(model)
    @test x[1, 1] == 1v[1]
    @test x[1, 2] == v[2] + v[4] * im
    @test x[2, 2] == 1v[3]
    # FIXME needs https://github.com/jump-dev/JuMP.jl/pull/2899
    #@test x[2, 1] == conj(x[1, 2])
    @constraint(model, x[1, 1] + x[2, 2] * im == 1 + 2im)
end

function test_simple_equality()
    model = Model(GLPK.Optimizer)
    add_bridge(model, COI.Bridges.Constraint.SplitEqualToBridge)

    @variable(model, 0 <= x[1:3] <= 1)

    @objective(model, Max, x[1] + 2x[2] - x[3])

    # FIXME needs https://github.com/jump-dev/JuMP.jl/pull/2895
    #@constraint(model, (1 + im) * x[1] + (1 + im) * x[2] + (1 - 3im) * x[3] == 1.0)
    aff = (1 + im) * x[1] + (1 + im) * x[2] + (1 - 3im) * x[3]
    @constraint(model, aff == 1.0)

    # Optimize the model:

    optimize!(model)

    # Check the termination and primal status to see if we have a solution:

    println("Termination status : ", termination_status(model))
    println("Primal status      : ", primal_status(model))

    # Print the solution:

    obj_value = objective_value(model)
    x_values = value.(x)

    println("Objective value : ", obj_value)
    println("x values        : ", x_values)

    @test obj_value ≈ 1.25
    @test x_values ≈ [0.0, 0.75, 0.25]
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestJuMP.runtests()
