#  Copyright 2015, Vincent Leclere, Francois Pacaud and Henri Gerard
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#############################################################################
# Compare three ways to solve the following intraday step problem
# in stochastic optimal control (SOC) form, with state s_t, control u_t and noise xi_t
#
# Min   E [\sum_{t=1}^{TF} c_t u_t]
# s.t.    s_{t+1} = s_t + u_t, s_0 given
#         -q_t <= s_t <= smax + q_t
#         u_min <= u_t <= u_max
#         u_t chosen knowing c_1 .. c_t
#
#############################################################################

using StochDynamicProgramming, Clp #calling the library and a solver.

const SOLVER = ClpSolver() 			   # requires "using Clp"
#const SOLVER = CplexSolver(CPX_PARAM_SIMDISPLAY=0) # requires "using CPLEX"

const MAX_ITER = 2 # number of iterations of SDDP


######## Stochastic optimization problem parameters  ########
const N_STAGES = 24             # number of stages of the SP problem

const CONTROL_MAX = 5           # upper bound on the control
const CONTROL_MIN = -5           # lower bound on the control

const S0 = [0.]                 # initial stock
const SMAX = 50                 # step capacity

#### prices support

# reading data
f = open("hourlyprices.dat")
text = readlines(f);
n_data = length(text);
data = Array{Float64,1}(n_data);
for i in 1 : n_data
    data[i] =  parse(Float64,text[i][1:end-2])
end
close(f)

# sorting and reducing data
data = reshape(data, (24,96));
data = sort(data,2);
supports = Array{Float64,2}(24,3);
supports[:,1] = mean(data[:,1:30],2);
supports[:,2] = mean(data[:,31:62],2);
supports[:,3] = mean(data[:,63:end],2);

### Constructing Noises laws
c_laws = NoiseLaw[ NoiseLaw(supports[t,:], [1/3 1/3 1/3]) for t in 1:24 ] 

## Define the dynamic of the (inventory) stock
function dynamic_P(t, x, u, xi)
    return [x[1] + u[1]]
end

## Define the instantaneous cost (corresponding to each stage)
function cost_P(t, x, u, w)
    return w[1] * u[1]
end

######### Setting up the SOC problem
q = zeros(N_STAGES)

function Intraday_Primal(q)
    s_bounds = repmat([(-Inf,Inf)],1,N_STAGES)
    for t in 1:N_STAGES
        s_bounds[1,t]=(-q[t], SMAX - q[t])
    end
    u_bounds = [(CONTROL_MIN, CONTROL_MAX)]
    spmodel = LinearSPModel(N_STAGES,u_bounds,S0,cost_P,dynamic_P,c_laws)
    set_state_bounds(spmodel, s_bounds) 

    paramSDDP = SDDPparameters(SOLVER,
                               passnumber=10,
                               max_iterations=MAX_ITER)
    V, pbs = solve_SDDP(spmodel, paramSDDP, 5)
    lb_sddp = StochDynamicProgramming.get_lower_bound(spmodel, paramSDDP, V)
    return lb_sddp
end

#Intraday_Primal(q)

###########################################################################

# lambda = (alpha,beta,gamma,delta[1],delta[2])
#FIXME hacky use of noise for qt --> not working
function cost_D(t,z,lambda, w)
    res = -lambda[2]*q[t] +  lambda[3]*(q[t]-SMAX) - lambda[4]*CONTROL_MAX+lambda[5]*CONTROL_MIN
    if t >= 2 
        return res 
    else
        return res + z[1]*S0[1]
    end
end

function dynamic_D(t,z,lambda,ct)
    # z_{t+1} = alpha_t - beta_t + gamma_t
    return [lambda[1]-lambda[2]+lambda[3]]
end
## Contraintes : pour t intermédiaire
# alpha_t = z_{t-1} 
# alpha_t -delta_1 + delta_2 =c_t --> z[1] - lambda[4] + lambda[5] == ct[1]
# beta gamma delta >=0 --> Dans les bornes
#  pour t = T on ajoute (et K==0)
# alpha_T - beta_T + gamma_T =0 --> z[1]-lambda[2]+lambda[3] == 0
# pour t = 1 il faut libérer alpha, i.e. on a 
# lambda[1] - lambda[4] + lambda[5] == ct[1]


function cst(t,z,lambda,ct)
    if t == 24
        return [lambda[1]-z[1],lambda[1] - lambda[4] + lambda[5] - ct[1], lambda[1]-lambda[2]+lambda[3]]
    elseif t>= 2
        return [lambda[1]-z[1],lambda[1] - lambda[4] + lambda[5] - ct[1]]
    else 
        return [lambda[1] - lambda[4] + lambda[5] - ct[1]]
    end
end


function Intraday_Dual(q)
    z_bounds = repmat([(-100,100)],1,N_STAGES) #TODO should be Inf
    lambda_bounds = repmat([(-Inf,Inf);(0,Inf);(0,Inf);(0,Inf);(0,Inf)],1,N_STAGES)

    spmodel_d = LinearSPModel(N_STAGES,lambda_bounds,S0,cost_D,dynamic_D,c_laws, eqconstr=cst )
    set_state_bounds(spmodel_d, z_bounds) 
    paramSDDP = SDDPparameters(SOLVER,
                               passnumber=10,
                               max_iterations=MAX_ITER,
                               verbose=10)
    V, pbs = solve_SDDP(spmodel_d, paramSDDP, 5)
    lb_sddp = StochDynamicProgramming.get_lower_bound(spmodel_d, paramSDDP, V)
    return lb_sddp
end

Intraday_Dual(q)

			        # bounds on the state
#u_bounds = [(CONTROL_MIN, CONTROL_MAX)] # bounds on the controls
#spmodel = LinearSPModel(N_STAGES,u_bounds,S0,cost_t,dynamic,c_laws) # constructing the spmodel
#set_state_bounds(spmodel, s_bounds) 	# adding states bounds to the problem

########## Solving the problem via SDDP (Stochastic Dynamic Dual Programming)
#if run_sddp
#    tic()
#    println("Starting resolution by SDDP (Stochastic Dynamic Dual Programming)")
#    # 10 forward pass, stop at MAX_ITER
#    paramSDDP = SDDPparameters(SOLVER,
#                               passnumber=10,
#                               max_iterations=MAX_ITER)
#    V, pbs = solve_SDDP(spmodel, paramSDDP, 2) # display information every 2 iterations
#    lb_sddp = StochDynamicProgramming.get_lower_bound(spmodel, paramSDDP, V)
#    println("Lower bound obtained by SDDP: "*string(round(lb_sddp,4)))
#    toc(); println();
#end

