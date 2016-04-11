#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Compare different ways of solving a stock problem :
# Min   E [\sum_{t=1}^TF c_t u_t]
# s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
#         0 <= s_t <= 1
#         u_min <= u_t <= u_max
#         u_t choosen knowing xi_1 .. xi_t
#############################################################################
push!(LOAD_PATH, "../src")
using StochDynamicProgramming, JuMP, Clp, Distributions

const XI_MAX = 0.3
const XI_MIN = 0
const N_XI = 10
const N_STAGES = 5
const COSTS = rand(N_STAGES)
const x0 = 0.5
const MAX_ITER = 100 # maximum iteration of SDDP
const SENSIBILITY = 0
const FORWARD_PASS_NUMBER = 2 # maximum iteration of SDDP

proba = 1/N_XI*ones(N_XI) # uniform probabilities
xi_support = collect(linspace(XI_MIN,XI_MAX,N_XI))
xi = [ [proba'; xi_support'] for t in 1:N_STAGES-1]

solver = ClpSolver();

m = SPModel("Linear");

@addState(m, 0 <= x[1:N_STAGES] <= 1)

@addControl(m, 0 <= u[1:N_STAGES-1] <= 0.5)

@addNoise(m, xi)
@setStochObjective(m, Min, sum{COSTS[i]*u[i], i = 1:N_STAGES-1})
@addDynamic(m, x[i+1] = x[i] + u[i] - xi[i])

@addConstraintsdp(m, x[1]==x0)

solveInterface(m, "SDDP", solver, FORWARD_PASS_NUMBER, SENSIBILITY, MAX_ITER)

solveInterface(m, "SDP", "HD", [0.01], [0.01])




