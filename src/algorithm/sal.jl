mutable struct SAL <: TauLeapMethod
  # parameters
  end_time :: Float64
  ϵ :: Float64
  δ :: Float64
  β :: Float64

  # state variables
  t              :: Float64
  dxdt           :: Vector{Float64}
  drdt           :: Vector{Float64}
  events         :: Vector{Int}

  # statistics
  stats_tracked :: Bool
  stats :: Dict{Symbol,Int}

  function SAL(end_time::AbstractFloat, ϵ, δ, β, stats_tracked)
    new(end_time, ϵ, δ, β,
      0.0, Float64[], Float64[], Int[],
      stats_tracked,
      Dict{Symbol,Int}(
        :negative_excursions => 0,
        :contractions => 0,
        :leaping_steps => 0,
        :gillespie_steps => 0
    ))
  end
end

set_time!(algorithm::SAL, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

##### accessors #####
get_derivatives(x::SAL) = x.dxdt, x.drdt

function init!(x::SAL, Xt, r)
  c = length(Xt)
  d = size(stoichiometry(r), 2)

  setfield!(x, :dxdt,   zeros(Float64, c))
  setfield!(x, :drdt,   zeros(Float64, d))
  setfield!(x, :events, zeros(Int,     d))

  return nothing
end

function reset!(algorithm::SAL, Xt, r)
  dxdt = algorithm.dxdt
  drdt = algorithm.drdt
  algorithm.t = zero(algorithm.t)

  update_all_propensities!(r, Xt)
  mean_derivatives!(dxdt, r)
  time_derivatives!(drdt, Xt, r, dxdt)

  return nothing
end

function step!(algorithm::SAL, Xt, r)
  # unpack propensities
  a = propensities(r)

  # unpack SAL variables
  ϵ = algorithm.ϵ
  β = algorithm.β
  δ = algorithm.δ
  dxdt, drdt = get_derivatives(algorithm)
  events = algorithm.events

  if intensity(a) > 0
    τ = tau_leap(r, drdt, ϵ)
    τ = min(τ, end_time(algorithm) - get_time(algorithm))

    # if τ is too small, do a Gillespie update
    if τ < δ / intensity(a)
      if algorithm.stats_tracked
        algorithm.stats[:gillespie_steps] += 1
      end
      τ = randexp() / intensity(a)
      set_time!(algorithm, τ)

      if !done(algorithm)
        μ = select_reaction(a)
        fire_reaction!(Xt, r, μ)
        update_propensities!(r, Xt, μ)
      end
    # otherwise, proceed with SAL
    else
      if algorithm.stats_tracked
        algorithm.stats[:leaping_steps] += 1
      end
      τ = sal_update!(algorithm, Xt, r, τ)
      update_all_propensities!(r, Xt)
      set_time!(algorithm, τ)
    end

    # update SAL variables
    mean_derivatives!(dxdt, r)
    time_derivatives!(drdt, Xt, r, dxdt)
  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  return nothing
end

function sal_update!(algorithm, Xt, r, τ)
  β = algorithm.β
  _, drdt = get_derivatives(algorithm)
  events = algorithm.events

  generate_events!(events, r, τ, drdt)

  isbadleap = is_badleap(Xt, r, events)
  
  if isbadleap && algorithm.stats_tracked
    algorithm.stats[:negative_excursions] += 1
  end

  while isbadleap
    if algorithm.stats_tracked
      algorithm.stats[:contractions] += 1
    end
    contract!(events, β)
    τ = τ * β
    isbadleap = is_badleap(Xt, r, events)
  end

  fire_reactions!(Xt, r, events)

  return τ
end

function mean_derivatives!(dxdt, r::DenseReactionSystem)
  V = stoichiometry(r)
  a = propensities(r)

  for k in eachindex(dxdt)
    dxdt[k] = 0.0
    for j in eachindex(a)
      dxdt[k] = dxdt[k] + a[j] * V[k, j]
    end
  end

  return dxdt
end

function mean_derivatives!(dxdt, r::SparseReactionSystem)
  V = stoichiometry(r)
  a = propensities(r)

  Vj = nonzeros(V)
  ix = rowvals(V)
  c  = size(V, 2)

  fill!(dxdt, 0.0)

  for j in 1:c
    for k in nzrange(V, j)
      i = ix[k]
      dxdt[i] = dxdt[i] + a[j] * Vj[k]
    end
  end

  return dxdt
end

function time_derivatives!(drdt, Xt, r::DenseReactionSystem, dxdt)
  for i in eachindex(drdt)
    drdt[i]  = 0.0
    for k in eachindex(Xt)
      ∂r∂x_k = compute_mass_action_deriv(Xt, r, i, k)
      drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
    end
  end
  return drdt
end

function time_derivatives!(drdt, Xt, r::SparseReactionSystem, dxdt)
  U  = coefficients(r)
  ix = rowvals(U)

  for i in eachindex(drdt)
    drdt[i]  = 0.0
    for k in nzrange(U, i)
      ∂r∂x_k = compute_mass_action_deriv(Xt, r, i, ix[k])
      drdt[i] = drdt[i] + ∂r∂x_k * dxdt[ix[k]]
    end
  end

  return drdt
end

function tau_leap(r::AbstractReactionSystem, drdt::Vector, ϵ::AbstractFloat)
  τ = Inf
  a = propensities(r)
  k = scaled_rates(r) # should we be using the unscaled rate instead?

  for j in eachindex(a)
    A = ϵ * max(a[j], k[j])
    B = abs(drdt[j])

    τ = min(τ, A / B)
  end

  return τ
end

function generate_events!(events, r, τ, drdt)
  a = propensities(r)

  for j in eachindex(a)
    λ = τ * a[j] + 0.5 * τ * τ * drdt[j]

    events[j] = poisrand(max(λ, 0))
  end

  return nothing
end

function is_badleap(Xt, r::DenseReactionSystem, events)
  V = stoichiometry(r)

  for i in eachindex(Xt)
    xi = Xt[i]
    for j in eachindex(events)
      xi = xi + events[j] * V[i, j]
      if xi < 0 return true end
    end
  end

  return false
end

function is_badleap(Xt, r::SparseReactionSystem, events)
  V  = stoichiometry(r)
  Vj = nonzeros(V)
  ix = rowvals(V)
  c  = size(V, 2)

  for j in 1:c
    xi = 0
    for k in nzrange(V, j)
      xi = Xt[ix[k]] + events[j] * Vj[k]

      if xi < 0 return true end
    end
  end

  return false
end

function contract!(events, β)
  for i in eachindex(events)
    k = 0
    for j in 1:events[i]
      k = rand() < β ? k + 1 : k
    end
    events[i] = k
  end

  return nothing
end
