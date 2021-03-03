module PGField

export GFElem, PGInf, order

struct GFElem{p}
    value::Int
    GFElem{p}(value::Integer) where p = new{p}(mod(value, p))
end
struct PGInf{p} end
Base.zero(::Union{GFElem{p},PGInf{p}}) where p = GFElem{p}(0)
Base.one(::Union{GFElem{p},PGInf{p}}) where p = GFElem{p}(1)
Base.iszero(x::GFElem) = iszero(x.value)
Base.iszero(::PGInf) = false

# Addition
Base.:+(x::GFElem{p},y::GFElem{p}) where p = GFElem{p}(widen(x.value)+widen(y.value))
Base.:+(x::GFElem{p},∞::PGInf{p}) where p = ∞
Base.:+(∞::PGInf{p},x::GFElem{p}) where p = x+∞
Base.:+(∞::PGInf{p},::PGInf{p}) where p = zero(∞)

# Multiplication
Base.:*(x::GFElem{p},y::GFElem{p}) where p = GFElem{p}(widemul(x.value,y.value))
Base.:*(x::GFElem{p},∞::PGInf{p}) where p = iszero(x) ? one(x) : ∞
Base.:*(∞::PGInf{p},x::GFElem{p}) where p = x*∞
Base.:*(∞::PGInf{p},::PGInf{p}) where p = ∞  # is this right?

# Negation and inverse
Base.:-(x::GFElem{p}) where p = GFElem{p}(-x.value)
Base.:-(∞::PGInf) = ∞
function Base.inv(x::GFElem{p}) where p  # based on https://github.com/Nemocas/AbstractAlgebra.jl/blob/master/src/julia/GF.jl#L270-L277
    iszero(x) && return PGInf{p}()
    g, s, t = gcdx(x.value, p)
    g != 1 && error("Characteristic not prime in GF($p).")
    return GFElem{p}(s)
end
Base.inv(∞::PGInf) = zero(∞)

# Powers
Base.:^(x::GFElem{p},n) where p = GFElem{p}(BigInt(x.value)^n)

# Order of element
order(x::GFElem{p}) where p = length(unique(x^i for i in 1:p))

# Broadcasting
Base.broadcastable(x::GFElem) = Ref(x)
Base.broadcastable(x::PGInf) = Ref(x)

end # module
