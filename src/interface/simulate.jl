"""
```
simulate(model::Network, algname::T, [output_type::Val{S} = Val(:fixed)]; time::Float64=1.0, epochs::Int=1, trials::Int=1, track_stats::Bool=false, kwargs...) where {T,S}
```

Simulate a `model` with `algname`. The simulation routine will run until the termination `time` for the given number of `trials`.

### Arguments

- `model`: The `Network` to simulate.
- `algname`: The name of the algorithm to carry out the simulation. One of `Direct()`, `FirstReaction()`, `NextReaction()`, `OptimizedDirect()`, `TauLeaping()`, or `StepAnticipation()`.
- `output_type = Val(:fixed)`: The type of output to record. The `Val(:full)` option records every simulation step, whereas the `Val(:fixed)` option records a fixed number of `epochs`.

### Optional Arguments

- `time = 1.0`: The amount of time to simulate the model, in units of the model.
- `epochs = 1`: The number of times to sample the vector of counts.
- `trials = 1`: The number of independent realizations to generate.
- `track_stats = false`: An option to toggle tracking algorithm statistics (e.g. number of steps).
- `kwargs`: Additional keyword arguments specific to each algorithm. Reference a specific `algname` for details.

"""
function simulate(model::Network, algname::T, output_type::Val{S}=Val(:fixed);
  time::Float64=1.0,
  epochs::Int=1,
  trials::Int=1,
  track_stats::Bool=false,
  kwargs...) where {T,S}

  # build algorithm
  algorithm = build_algorithm(algname, time, track_stats; kwargs...)

  # extract model information
  c = n_species(model)
  d = n_reactions(model)

  species   = species_list(model)
  reactions = reaction_list(model)

  # create simulation data structures
  x0, rxn, id, id2ind = make_datastructs(species, reactions, c, d)

  # get output type
  output = build_output(output_type, c, epochs, trials, time)

  # initialize
  xt = copy(x0)
  init!(algorithm, xt, rxn)

  # delegate trials
  simulate_wrapper!(output, xt, x0, algorithm, rxn)

  result = SimulationSummary(model, algname, algorithm, time, epochs, trials, id2ind, output; kwargs...)

  return result
end

function build_output(::Val{:fixed}, nspecies, epochs, ntrials, tfinal)
  n = epochs + 1
  tdata = collect(linspace(0.0, tfinal, n))
  output = RegularEnsemble(ntrials, nspecies, epochs)
  for i in eachindex(output)
    @inbounds copy!(output[i].tdata, tdata)
  end
  return output
end

function build_output(::Val{:full}, nspecies, epochs, ntrials, tfinal)
  return Ensemble(ntrials)
end

function make_datastructs(species, reactions, c, d)
  # state vector
  x0, id, id2ind = make_species_vector(species)

  # reactions
  if d <= 8
    rxn = DenseReactionSystem(reactions, id2ind, c, d)
  else
    rxn = SparseReactionSystem(reactions, id2ind, c, d)
  end

  return x0, rxn, id, id2ind
end

function simulate_wrapper!(output, xt, x0, alg, rxn)
  N = Threads.nthreads()

  if N == 1
    simulate_chunk!(output, xt, x0, alg, rxn, 1:length(output))
  else
    Threads.@threads for i in 1:N
      len = div(length(output), N)
      domain = ((i-1)*len+1):i*len
      simulate_chunk!(output, deepcopy(xt), x0, deepcopy(alg), deepcopy(rxn), domain)
    end
  end
end

function simulate_chunk!(output, Xt, X0, algorithm, reactions, trial_set)
  a = propensities(reactions)
  for trial in trial_set
    copy!(Xt, X0)
    reset!(algorithm, Xt, reactions)

    xw = output[trial]

    simulate!(xw, Xt, algorithm, reactions)
  end

  return output
end

function simulate!(xw :: RegularPath, Xt, algorithm, reactions)
  epoch = 1
  while !done(algorithm)
    epoch = update!(xw, algorithm.t, Xt, epoch)
    step!(algorithm, Xt, reactions)
  end
  update!(xw, algorithm.t, Xt, epoch)

  return xw
end

function simulate!(xw :: SamplePath, Xt, algorithm, reactions)
  while !done(algorithm)
    update!(xw, algorithm.t, Xt)
    step!(algorithm, Xt, reactions)
  end

  return xw
end
