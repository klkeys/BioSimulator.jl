"""
```
PropensityVector(m::Integer)
```

An m-vector that caches reaction propensities and maintains their running sum.

### Internals
`cache`: The cached reaction propensities.
`intensity`: The running sum of cached propensities.
`error_bound`: A number that quantifies the error resulting from updating sums with Kahan summation. Once this value exceeds a threshold, the running sum is explicitly recalculated.
"""
mutable struct PropensityVector{T <: AbstractFloat} <: AbstractVector{T}
  cache       :: Vector{T}
  intensity   :: T
  error_bound :: T

  PropensityVector{T}(m::Integer) where T = new(zeros(T, m), zero(T), zero(T))
end

PVec{T} = PropensityVector{T}

Base.eltype(x::PVec{T}) where T = T

##### PVec interface #####
intensity(x::PVec) = x.intensity
isstable{T}(x::PVec{T})  = (x.error_bound <= eps(T) * x.intensity)

@fastmath function update_errorbound!{T}(x::PVec{T}, xi::T, i::Integer)
  x.error_bound = x.error_bound + eps(T) * (x.intensity + x[i] + xi)
end

@fastmath function update_intensity!{T}(x::PVec{T}, xi::T, i::Integer)
  # if x.intensity - x[i] + xi < 0
  #   error("""
  #   Intensity update failed!
  #   old a₀ = $(x.intensity)
  #   new a₀ = $(x.intensity - x[i] + xi)
  #
  #   old x[$i] = $(x[i]) / $(typeof(x[i]))
  #   new x[$i] = $(xi) / $(typeof(xi))
  #
  #   error_bound = $(x.error_bound)
  #   is_stable?  = $(isstable(x))
  #   """)
  # end
  x.intensity = x.intensity - x[i] + xi
end

##### AbstractArray interface #####

Base.size(x::PVec)             = size(x.cache)
Base.size(x::PVec, d::Integer) = size(x.cache, d)

Base.getindex(x::PVec, i::Integer) = x.cache[i]

Base.setindex!(x::PVec, xi, i::Integer) = (x.cache[i] = xi)

Base.IndexStyle(x::PVec) = Base.IndexStyle(x.cache)
