import LinearAlgebra
import JuMP

function add_all_bridges(model::JuMP.Model)
    JuMP.add_bridge(model, Bridges.Variable.HermitianToSymmetricPSDBridge)
    JuMP.add_bridge(model, Bridges.Constraint.SplitEqualToBridge)
    JuMP.add_bridge(model, Bridges.Constraint.SplitZeroBridge)
end

struct HermitianPSDCone end

struct HermitianMatrixShape <: JuMP.AbstractShape
    side_dimension::Int
end

function JuMP.vectorize(matrix::Matrix, ::HermitianMatrixShape)
    n = LinearAlgebra.checksquare(matrix)
    return vcat(
        JuMP.vectorize(_real.(matrix), JuMP.SymmetricMatrixShape(n)),
        JuMP.vectorize(_imag.(matrix[1:(end-1),2:end]), JuMP.SymmetricMatrixShape(n - 1)),
    )
end

function JuMP.reshape_vector(v::Vector{T}, shape::HermitianMatrixShape) where T
    NewType = MA.promote_operation(MA.add_mul, T, Complex{Bool}, T)
    n = shape.side_dimension
    matrix = Matrix{NewType}(undef, n, n)
    real_k = 0
    imag_k = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    for j in 1:n
        for i in 1:(j-1)
            real_k += 1
            imag_k += 1
            matrix[i, j] = v[real_k] + im * v[imag_k]
            matrix[j, i] = v[real_k] - im * v[imag_k]
        end
        real_k += 1
        matrix[j, j] = v[real_k]
    end
    return matrix
end

function _mapinfo(f::Function, v::JuMP.ScalarVariable)
    info = v.info
    return JuMP.ScalarVariable(JuMP.VariableInfo(
        info.has_lb,
        f(info.lower_bound),
        info.has_ub,
        f(info.upper_bound),
        info.has_fix,
        f(info.fixed_value),
        info.has_start,
        f(info.start),
        info.binary,
        info.integer,
    ))
end

_real(s::String) = s
_imag(s::String) = s

_real(v::JuMP.ScalarVariable) = _mapinfo(real, v)
_imag(v::JuMP.ScalarVariable) = _mapinfo(imag, v)
_conj(v::JuMP.ScalarVariable) = _mapinfo(conj, v)
function _isreal(v::JuMP.ScalarVariable)
    info = v.info
    return isreal(info.lower_bound) && isreal(info.upper_bound) && isreal(info.fixed_value) && isreal(info.start)
end

function _vectorize_variables(_error::Function, matrix::Matrix)
    n = LinearAlgebra.checksquare(matrix)
    for j in 1:n
        if !_isreal(matrix[j, j])
            _error(
                "Non-real bounds or starting values for diagonal of hermitian variable.",
            )
        end
        for i in 1:j
            if matrix[i, j] != _conj(matrix[j, i])
                _error(
                    "Non-conjugate bounds, integrality or starting values for hermitian variable.",
                )
            end
        end
    end
    JuMP.vectorize(matrix, HermitianMatrixShape(n))
end

function JuMP.build_variable(
    _error::Function,
    variables::Matrix{<:JuMP.AbstractVariable},
    ::HermitianPSDCone,
)
    n = JuMP._square_side(_error, variables)
    set = HermitianPositiveSemidefiniteConeTriangle(n)
    shape = HermitianMatrixShape(n)
    return JuMP.VariablesConstrainedOnCreation(
        _vectorize_variables(_error, variables),
        set,
        shape,
    )
end

