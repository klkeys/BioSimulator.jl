abstract Coupling

immutable Tight <: Coupling end
immutable Loose <: Coupling end

export Tight, Loose

type ODM <: Algorithm
  itr::Int

  tf::Float64
  dt::Float64

  c::Coupling
  pre_steps::Int
  samples::Int

  t::Float64
  intensity::Float64
  steps::Int

  g::DiGraph

  function ODM(itr, tf, dt; kwargs...)
    args    = Dict{Symbol,Any}(kwargs)
    c       = get(args, :c,       Tight())
    steps   = get(args, :steps,   100)
    samples = get(args, :samples, 1)

    new(itr, tf, dt, c, steps, samples, 0.0, 0.0, 0, DiGraph())
  end
end

function init(alg::ODM, rxns, spcs, initial, params)
  alg.g = typeof(alg.c) <: Loose ? init_dep_graph(rxns) : DiGraph()

  # Presimulate to sort reactions according to multiscale property.
  # This will modify spcs and rxns
  init_odm!(spcs, rxns, params, initial, alg.pre_steps, alg.samples)
  return;
end

function reset(alg::ODM, rxns, spcs, params)
  alg.t     = 0.0
  alg.steps = 0
  alg.intensity = compute_propensities!(rxns, spcs, params)

  return;
end

function step(alg::ODM, rxns, spcs, params)
  τ, intensity = odm_update!(alg.c, spcs, rxns, params, alg.intensity, alg.g)
  alg.t = alg.t + τ
  alg.steps = alg.steps + 1
  alg.intensity = intensity
  return;
end

function odm_update!(c::Coupling, spcs::Vector{Int}, rxns::Vector{Reaction}, param, intensity::Float64, g::LightGraphs.DiGraph)
	τ = rand(Exponential(1 / intensity))
	jump = intensity * rand()
	μ = sample(rxns, jump)
	μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")

	intensity = update_propensities!(c, rxns, spcs, param, g)
	return τ, intensity
end

function update_propensities!(c::Loose, rxns, spcs, param)
	dependents = neighbors(g, μ)

	@inbounds for α in dependents
		intensity = intensity - rxns[α].propensity
		propensity!(rxns[α], spcs, param);
		intensity = intensity + rxns[α].propensity
	end
	return intensity
end

function update_propensities!(c::Tight, rxns, spcs, param, g)
	return compute_propensities!(rxns, spcs, param)
end

function presimulate!(spcs, rxns, params, initial, n, itr)
	events = zeros(Float64, length(rxns))

	for i = 1:itr
		copy!(spcs, initial)
		for k = 1:n
			intensity = compute_propensities!(rxns, spcs, params)
			τ = rand(Exponential(1 / intensity))
			u = rand()
		  jump = intensity * u
		  j = sample(rxns, jump)
		  j > 0 ? update!(spcs, rxns[j]) : error("No reaction occurred!")
			events[j] = events[j] + 1
		end
	end

	for i in eachindex(rxns)
		rxns[i].propensity = events[i]
	end

	return rxns
end

function init_odm!(spcs, rxns, params, initial, n, itr)
	rxns = presimulate!(spcs, rxns, params, initial, n, itr)
	sort!(rxns, alg=Base.MergeSort, lt=isless, rev=true)
	return rxns
end

function isless(x::Reaction, y::Reaction)
	return x.propensity < y.propensity
end
