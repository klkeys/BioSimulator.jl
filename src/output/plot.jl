@recipe function f(tr :: Trajectory)
  id = tr.id
  t  = tr.t_index
  x  = tr.val

  # global attributes
  legend     --> true
  grid       --> false
  xguide     --> "time"
  yguide     --> "population size"
  xlims      --> (t[1], t[end])
  seriestype --> :steppre
  label      --> id

  t, x
end

@recipe function f(mt :: MeanTrajectory)
  id = mt.id
  t  = mt.t_index
  μ  = mt.mean_val
  σ  = mt.std_val

  # global attributes
  legend     --> true
  grid       --> false
  xguide     --> "time"
  yguide     --> "mean population size"
  xlims      --> (t[1], t[end])
  fillalpha  --> 0.3
  seriestype --> :path
  label      --> id
  ribbon     --> σ

  t, μ
end

@recipe function f(hg :: Histogram)
  id = hg.id
  t  = hg.t
  x  = hg.val

  legend     --> true
  grid       --> false
  xguide     --> "population size"
  yguide     --> "frequency"
  seriestype --> :histogram

  x
end

immutable PhaseTrajectory
  ids :: Vector{String}
  x   :: Vector{Int}
  y   :: Vector{Int}
end

function PhaseTrajectory(output :: SimData, species1, species2, trial)
  t_index = output.t_index

  x = reshape(output[species1][:, trial], length(t_index))
  y = reshape(output[species2][:, trial], length(t_index))

  return PhaseTrajectory([species1, species2], x, y)
end

@recipe function f(pt :: PhaseTrajectory)
  ids = pt.ids
  x   = pt.x
  y   = pt.y

  # global attributes
  legend     --> false
  grid       --> false
  xguide     --> "$(ids[1]) population size"
  yguide     --> "$(ids[2]) population size"
  seriestype --> :path
  label      --> ids

  x, y
end
