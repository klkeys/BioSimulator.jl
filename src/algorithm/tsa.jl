mutable struct TSA <: ExactMethod
  # parameters
  end_time :: Float64
  l_step :: Int

  # state
  t :: Float64
  p_cache :: Vector{Float64}

  # statistics

  function TSA(end_time, l)
    new(end_time, l, 0.0, zeros(l),
      Dict{Symbol,Int}(
        :tsa_steps => 0
      )
      )
  end
end

TSA(end_time; l :: Int = 10) = TSA(end_time, l)

set_time!(algorithm :: TSA, τ :: AbstractFloat) = (algorithm.t = algorithm.t + τ)

##### implementation #####

function step!(algorithm :: TSA, Xt :: Vector, r :: AbstractReactionSystem)
  l_step  = algorithm.l_step
  p_cache = algorithm.p_cache

  a = propensities(r)

  if intensity(a) > 0
    # simulate a total of l_step events
    for l in 1:l_step
      # store the intensity
      p_cache[l] = intensity(a)

      # select and fire a reaction
      μ = select_reaction(a)
      fire_reaction!(Xt, r, μ)

      # update the propensities
      update_propensities!(a, r, Xt, μ)
    end
    # impute_time()
  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(Error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  algorithm.stats[:tsa_steps] += 1

  return nothing
end