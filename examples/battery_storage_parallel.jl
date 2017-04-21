#  Copyright 2017, V.Leclere, H.Gerard, F.Pacaud, T.Rigaut
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Solving an decision hazard battery storage problem using multiprocessed
# dynamic programming.
#
# to launch (with N processes)
#       $ julia -p N battery_storage_parallel.jl
#
#       PROBLEM DESCRIPTION
# For t = 1..T, we want to satisfy a deterministic energy demand d_t using:
# - a random renewable energy unit with production (at time t) W_t
# - a battery with stock (at time t) S_t
# - a quantity G_{t} bought on the market for a deterministic price c_t
# The decision variable is the quantity U_t of energy stored (or discharged) in
# the battery, is decided knowing the uncertainty W_0,...,W_t.
# We want to minimize the money spent on buying electricity (cost are deterministic) :
#               min     E[\sum_{t=0}^T c_t * G_{t}]
#                       d_t + U_t <= W_t + G_t                          (a)
#                       S_{t+1} = S_t + r (U_t)^+ + 1/r (U_t)^-         (b)
#                       Smin <= S_t <= Smax                             (c)
#                       G_t >= 0                                        (d)
#
# (a) : more production than demand (excess is wasted)
# (b) : dynamic of the stock with a charge coefficient r
# (c) : limits on the battery
# (d) : we don't sell electricity




import StochDynamicProgramming, Distributions
println("library loaded")

# We have to define the instance on all the workers (processes)
@everywhere begin

    ######## Stochastic Model  Parameters  ########
    srand(1234) #Seed is essential when defining random variables accross multiple workers
    const N_STAGES = 50
    const COSTS = rand(N_STAGES)
    const DEMAND = rand(N_STAGES)

    const CONTROL_MAX = 0.5
    const CONTROL_MIN = 0

    const STATE_MAX = 1
    const STATE_MIN = 0

    const XI_MAX = 0.3
    const XI_MIN = 0
    const N_XI = 10
    # initial stock
    const S0 = 0.5

    # charge and discharge efficiency parameters
    const rho_c = 0.98
    const rho_dc = 0.97

    # create law of noises
    const proba = 1/N_XI*ones(N_XI) # uniform probabilities
    const xi_support = collect(linspace(XI_MIN,XI_MAX,N_XI))
    const xi_law = StochDynamicProgramming.NoiseLaw(xi_support, proba)
    const xi_laws = StochDynamicProgramming.NoiseLaw[xi_law for t in 1:N_STAGES-1]

    # Define dynamic of the stock:
    function dynamic(t, x, u, xi)
    	return [ x[1] + 1/rho_dc * min(u[1],0) + rho_c * max(u[1],0) ]
    end

    # Define cost corresponding to each timestep:
    function cost_t(t, x, u, xi)
        return COSTS[t] * max(0, DEMAND[t] + u[1] - xi[1])
    end

    ######## Setting up the SPmodel
    s_bounds = [(STATE_MIN, STATE_MAX)]
    u_bounds = [(CONTROL_MIN, CONTROL_MAX)]
    spmodel = StochDynamicProgramming.StochDynModel(N_STAGES,
                                                    u_bounds,
                                                    [S0],
                                                    cost_t,
                                                    dynamic,
                                                    xi_laws,
                                                    xbounds = s_bounds)

    scenarios = StochDynamicProgramming.simulate_scenarios(xi_laws,1000)

    stateSteps = [0.01]
    controlSteps = [0.001]

    paramSDP = StochDynamicProgramming.ExhaustiveSdpParameters(stateSteps,
                                                    controlSteps,
                                                    infoStructure = "HD")
end

@time  Vs = StochDynamicProgramming.solve_dp(spmodel,paramSDP, 1)

lb_sdp = StochDynamicProgramming.get_bellman_value(spmodel,paramSDP,Vs)
println("Value obtained by SDP: "*string(lb_sdp))
@time costsdp, states, stocks = StochDynamicProgramming.forward_simulations(spmodel,paramSDP,Vs,scenarios)
println(mean(costsdp))

