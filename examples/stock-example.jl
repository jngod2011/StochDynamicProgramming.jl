#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#############################################################################
# Compare three ways to solve the following optimal inventory problem
# in stochastic optimal control (SOC) form, with state s_t, control u_t and noise xi_t
#
# Min   E [\sum_{t=1}^{TF} c_t u_t]
# s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
#         0 <= s_t <= 1
#         u_min <= u_t <= u_max
#         u_t chosen knowing xi_1 .. xi_t
#
#############################################################################

using StochDynamicProgramming, Clp #calling the library and a solver.

######## Switching on or off the three resolution methods
run_sddp = true # false if you don't want to run sddp
run_sdp  = true # false if you don't want to run sdp
run_ef   = true # false if you don't want to run extensive formulation
test_simulation = true # false if you don't want to simulate the strategies on scenarios

######## Optimization parameters  ########
# choose the LP solver used.
SOLVER = ClpSolver() 			   # require "using Clp"
#const SOLVER = CplexSolver(CPX_PARAM_SIMDISPLAY=0) # require "using CPLEX"

# convergence test
MAX_ITER = 10 # number of iterations of SDDP
step = 0.01   # discretization step of SDP

######## Stochastic Model  Parameters  ########
N_STAGES = 6              # number of stages of the SP problem
COSTS = [sin(3*t)-1 for t in 1:N_STAGES-1]
#const COSTS = rand(N_STAGES)    # randomly generating deterministic costs

CONTROL_MAX = 0.5         # bounds on the control
CONTROL_MIN = 0

XI_MAX = 0.3              # bounds on the noise
XI_MIN = 0
N_XI = 10                 # discretization of the noise

S0 = 0.5                  # initial stock

# Specify the marginal probability distributions of the noises xi_1 ... xi_{N_STAGES}
proba = 1/N_XI*ones(N_XI) # uniform probabilities
xi_support = collect(linspace(XI_MIN,XI_MAX,N_XI))
xi_law = NoiseLaw(xi_support, proba)
# By assumption, the process xi_1 ... xi_{N_STAGES -1} is made of independent random variables
xi_laws = NoiseLaw[xi_law for t in 1:N_STAGES-1]

# Define the dynamic of the (inventory) stock
function dynamic(t, x, u, xi)
    return [x[1] + u[1] - xi[1]]
end

# Define the instantaneous cost (corresponding to each stage)
function cost_t(t, x, u, w)
    return COSTS[t] * u[1]
end

######## Setting up the SOC problem
s_bounds = [(0, 1)] 			        # bounds on the state
u_bounds = [(CONTROL_MIN, CONTROL_MAX)] # bounds on the controls
spmodel = LinearSPModel(N_STAGES,u_bounds,S0,cost_t,dynamic,xi_laws) # constructing the spmodel
set_state_bounds(spmodel, s_bounds) 	# adding states bounds to the problem

######### Solving the problem via SDDP (Stochastic Dynamic Dual Programming)
if run_sddp
    tic()
    println("Starting resolution by SDDP (Stochastic Dynamic Dual Programming)")
    # 10 forward pass, stop at MAX_ITER
    paramSDDP = SDDPparameters(SOLVER,
                               passnumber=10,
                               max_iterations=MAX_ITER)
    sddp = solve_SDDP(spmodel, paramSDDP, 2) # display information every 2 iterations
    lb_sddp = StochDynamicProgramming.get_lower_bound(spmodel, paramSDDP, sddp.bellmanfunctions)
    println("Lower bound obtained by SDDP: "*string(round(lb_sddp,4)))
    toc(); println();
end

######### Solving the problem via SDP (Stochastic Dynamic Programming)
if run_sdp
    tic()
    println("Starting resolution by SDP (Stochastic Dynamic Programming)")
    stateSteps = [step] # discretization step of the state
    controlSteps = [step] # discretization step of the control
    println("The discretization step of the state and of the control is "*string(step))
    infoStruct = "HD" # noise at time t is known before making the decision at time t
#   infoStruct = "DH" # noise at time t is known after taking the decision at time t
    paramSDP = SDPparameters(spmodel, stateSteps, controlSteps, infoStruct)
    Vs = solve_dp(spmodel,paramSDP, 1)
    value_sdp = StochDynamicProgramming.get_bellman_value(spmodel,paramSDP,Vs)
    println("Value obtained by SDP: "*string(round(value_sdp,4)))
    toc(); println();
end

######### Solving the problem via Extensive Formulation
if run_ef
    tic()
    println("Starting resolution by EF (Extensive Formulation)")
    value_ef = extensive_formulation(spmodel, paramSDDP)[1]
    println("Value obtained by EF (true value of the problem): "*string(round(value_ef,4)))
    solve_sdp && println("Relative error of SDP value: "*string(100*round(abs(value_sdp/value_ef)-1,4))*"%")
    solve_sddp && println("Relative error of SDDP lower bound: "*string(100*round(abs(lb_sddp/value_ef)-1,4))*"%")
    toc(); println();
end

######### Evaluating and comparing the three resolution methods on 1000 scenarios
#srand(1234) # fix the random seed accross runs
if test_simulation
    println("Evaluating the resolution methods on 1000 scenarios (the lower the better)")
    scenarios = StochDynamicProgramming.simulate_scenarios(xi_laws,1000)
    if run_sddp
        costsddp, stocks = forward_simulations(spmodel, paramSDDP, pbs, scenarios)
        println("Value obtained by SDDP: "*string(round(mean(costsddp),4)))
    end
    if run_sdp
        costsdp, states, controls = sdp_forward_simulation(spmodel,paramSDP,scenarios,Vs)
        println("Value obtained by SDP: "*string(round(mean(costsdp),4)))
    end
end
