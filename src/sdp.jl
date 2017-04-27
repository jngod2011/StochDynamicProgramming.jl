#  Copyright 2017, V.Leclere, H.Gerard, F.Pacaud, T.Rigaut
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
#  Stochastic dynamic programming algorithm
#
#############################################################################

using ProgressMeter, Interpolations


"""
Compute interpolation of the value function at time t

# Arguments
* `model::SPmodel`:
* `dim_states::Int`:
    the number of state variables
* `v::Array`:
    the value function to interpolate
* `time::Int`:
    time at which we have to interpolate V

# Return
* Interpolation
    the interpolated value function (working as an array with float indexes)

"""
function value_bspline_interpolation( dim_states::Int, V::Union{SharedArray, Array}, time::Int)
    return interpolate(V[[Colon() for i in 1:dim_states]...,time], BSpline(Linear()), OnGrid())
end

function convex_inner_approximation( dim_states::Int, V::Union{SharedArray, Array}, time::Int)
    points = 0
    values = 0
    return points, values
end


"""
Compute the cartesian products of discretized state spaces

# Arguments
* `model::SPmodel`:
    the model of the problem
* `param::SdpParameters`:
    the parameters of the problem

# Return
* Iterators: product_states
    the cartesian product iterators for states

"""
function generate_state_grid(model::SPModel, param::SdpParameters, t::Int = -1 )

    xlim = model.xlim
    xdim = model.dimStates
    xsteps = param.stateSteps
    if (ndims(xlim)==2)
        if (ndims(xsteps)==1)
            error("the discretization of state grid should be indexed by time")
        end
        if t<0
            warn("Time step is not provided for state grid. Defaulted to 1")
            return collect(Base.product([xlim[i,1][1]:xsteps[i,1]:xlim[i,1][2] for i in 1:xdim]...))
        else
            return collect(Base.product([xlim[i,t][1]:xsteps[i,1]:xlim[i,t][2] for i in 1:xdim]...))
        end
    else
        return collect(Base.product([xlim[i][1]:xsteps[i]:xlim[i][2] for i in 1:xdim]...))
    end

end

"""
Compute the cartesian products of discretized control spaces or more complex space if provided

# Arguments
* `model::SPmodel`:
    the model of the problem
* `param::SdpParameters`:
    the parameters of the problem
* `t::Int`:
    time step of the value function computation
* `x::Array{Float64}`:
    the  current state explored

# Return
* Iterators: product_states and product_controls
    the cartesian product iterators for both states and controls

"""
function generate_control_grid(model::SPModel, param::SdpParameters,
                                t::Nullable{Int} = Nullable{Int}(),
                                x::Nullable{Array} = Nullable{Array}(),
                                w::Nullable{Array} = Nullable{Array}())

    if (isnull(param.buildSearchSpace))||(isnull(t))||(isnull(x))
        product_controls = Base.product([model.ulim[i][1]:param.controlSteps[i]:model.ulim[i][2] for i in 1:model.dimControls]...)
    else
        product_controls = param.buildSearchSpace(t, x, w)
    end

    return collect(product_controls)
end


"""
Dynamic programming algorithm to compute optimal value functions
by backward induction using bellman equation in the finite horizon case.
The information structure can be Decision Hazard (DH) or Hazard Decision (HD)

# Arguments
* `model::SPmodel`:
    the DPSPmodel of our problem
* `param::SdpParameters`:
    the parameters for the SDP algorithm
* `display::Int`:
    the output display or verbosity parameter

# Return
* `value_functions::Array`:
    the vector representing the value functions as functions of the state
    of the system at each time step

"""
function solve_dp(model::SPModel, param::SdpParameters, display=0::Int64)
    # Start of the algorithm
    V = compute_value_functions_grid(model, param, display)
    return V
end

function build_cost_function(costFunctions::Union{Nullable{Function},Function,Array{Function},PolyhedralFunction})
    if isa(costFunctions, Function)
        return costFunctions
    else
        return function cost(t,x,u,w)
            return maximum([aff_func(t,x,u,w) for aff_func in costFunctions])
        end
    end
end

function build_final_cost_function(costFunctions::Union{Function, Array{Function}, PolyhedralFunction})
    if isa(costFunctions, Function)
        return costFunctions
    elseif isa(costFunctions, Array{Function})
        return (x) -> maximum([aff_func(x) for aff_func in costFunctions])
    elseif isa(costFunctions, PolyhedralFunction)
        return (x) -> maximum([costFunctions.betas[k] + dot(costFunctions.lambdas[k,:],[x...]) for k in 1:costFunctions.numCuts])
    else
        return (x) -> 0
    end
end

function initialize_final_value!(final_cost, x_grid::Array, x_bounds,
                                x_steps, V::Union{Array, SharedArray})
    for x in x_grid
        ind_x = BellmanSolvers.index_from_variable(x, x_bounds, x_steps)
        V[ind_x..., end] = final_cost(x)
    end

end


function build_contraints_function(ineq_cons::Nullable{Function}, eq_cons::Nullable{Function})
    if !isnull(ineq_cons)&&!isnull(eq_cons)
        return (t,x,u,w) -> (find(abs.(get(eq_cons)(t,x,u,w)).>1e-10)==[])&&(find(get(ineq_cons)(t,x,u,w).>1e-10)==[])
    elseif !isnull(ineq_cons)
        return (t,x,u,w) -> (find(get(ineq_cons)(t,x,u,w).>1e-10)==[])
    elseif !isnull(eq_cons)
        return (t,x,u,w) -> (find(get(eq_cons)(t,x,u,w).!=0)==[])
    else
        return (t,x,u,w) -> true
    end
end

function build_marginal_law(model, param, t)
    law = model.noises

    if (param.expectation_computation=="MonteCarlo")
        return (param.monteCarloSize, [sampling(law,t) for i in 1:sampling_size],
                (1/sampling_size)*ones(sampling_size))
    else
        return (law[t].supportSize, law[t].support, law[t].proba)
    end
end

"""
Dynamic Programming algorithm to compute optimal value functions

# Parameters
* `model::StochDynProgModel`:
    the StochDynProgModel of the problem
* `param::SdpParameters`:
    the parameters for the algorithm
* `display::Int`:
    the output display or verbosity parameter

# Returns
* `value_functions::Array`:
    the vector representing the value functions as functions of the state
    of the system at each time step

"""
function compute_value_functions_grid(model::SPModel,
                                        param::SdpParameters,
                                        display=0::Int64)

    TF = model.stageNumber

    u_bounds = model.ulim
    x_bounds = model.xlim
    x_steps = param.stateSteps
    x_dim = model.dimStates

    dynamics = model.dynamics

    cost = build_cost_function(model.costFunctions)

    fin_cost = build_final_cost_function(model.finalCost)

    constraints = build_contraints_function(model.inequalityConstraints,
                                            model.equalityConstraints)

    law = model.noises

    build_Ux = param.buildSearchSpace

    #Compute cartesian product spaces
    product_states = generate_state_grid(model, param)


    V = SharedArray{Float64}(zeros(Float64, size(product_states)..., TF))

    #Compute final value functions
    initialize_final_value!(fin_cost, product_states, x_bounds, x_steps, V)

    if param.infoStructure == "HD"
        get_V_t_x = BellmanSolvers.exhaustive_search_hd
    else
        param.infoStructure == "DH" || warn("Information structure defaulted to DH")
        param.infoStructure = "DH"
        get_V_t_x = BellmanSolvers.exhaustive_search_dh
    end

    #Construct a progress meter
    p = 0
    if display > 0
        p = Progress((TF-1), 1)
        println("Starting resolution by SDP")
    end

    product_controls = generate_control_grid(model, param)
    # Loop over time:
    for t = (TF-1):-1:1

        display == 0 || next!(p)

        sampling_size, samples, probas = build_marginal_law(model, param, t)

        Vitp = value_bspline_interpolation(x_dim, V, t+1)

        function compute_V_t_x!(x)
            ind_x = BellmanSolvers.index_from_variable(x, x_bounds, x_steps)
            V[ind_x...,t] = get_V_t_x(sampling_size, samples, probas,
                                        u_bounds, x_bounds, x_steps, x_dim,
                                        product_controls, dynamics,
                                        constraints, cost, Vitp, t,
                                        x, build_Ux)[1]

        end

        pmap(compute_V_t_x!, product_states)

    end
    return V
end

"""
Get the optimal value of the problem from the optimal Bellman Function

# Arguments
* `model::SPmodel`:
    the DPSPmodel of our problem
* `param::SdpParameters`:
    the parameters for the SDP algorithm
* `V::Array{Float64}`:
    the Bellman Functions

# Return
* `V_x0::Float64`:

"""
function get_bellman_value(model::SPModel, param::SdpParameters,
                            V::Union{SharedArray, Array})
    ind_x0 = BellmanSolvers.real_index_from_variable(model.initialState, model.xlim, param.stateSteps)
    Vi = value_bspline_interpolation(model.dimStates, V, 1)
    return Vi[ind_x0...,1]
end


"""
Get the optimal control at time t knowing the state of the system in the decision
hazard case

# Arguments
* `model::SPmodel`:
    the DPSPmodel of our problem
* `param::SdpParameters`:
    the parameters for the SDP algorithm
* `V::Array{Float64}`:
    the Bellman Functions
* `t::int`:
    the time step
* `x::Array`:
    the state variable
* `w::Array`:
    the alea realization

# Return
* `V_x0::Float64`:

"""
function get_control(model::SPModel,param::SdpParameters,
                     V, t::Int64, x::Array, w::Union{Void, Array} = nothing)

    args = []
    optional_args = []

    dynamics = model.dynamics

    cost = build_cost_function(model.costFunctions)

    constraints = build_contraints_function(model.inequalityConstraints,
                                            model.equalityConstraints)

    if w==nothing

        get_u = BellmanSolvers.exhaustive_search_dh_get_u

        push!(args, build_marginal_law(model, param, t)...)

        push!(optional_args, param.buildSearchSpace)

    else

        get_u = BellmanSolvers.exhaustive_search_hd_get_u

        push!(optional_args, w, param.buildSearchSpace)

    end

    push!(args, model.ulim, model.xlim, param.stateSteps,
            model.dimStates, generate_control_grid(model, param),
            dynamics, constraints, cost,
            value_bspline_interpolation(model.dimStates, V, t+1), t, x)

    return get_u(args..., optional_args...)[1]
end


"""
Simulation of optimal control given an initial state and multiple scenarios

# Arguments
* `model::SPmodel`:
    the DPSPmodel of our problem
* `param::SdpParameters`:
    the parameters for the SDP algorithm
* `scenarios::Array`:
    the scenarios of uncertainties realizations we want to simulate on
* `X0::SdpParameters`:
    the initial state of the system
* `V::Array`:
    the vector representing the value functions as functions of the state
    of the system at each time step
* `display::Bool`:
    the output display or verbosity parameter

# Return
* `costs::Array{Float}`:
    the cost of the optimal control over the scenario provided
* `states::Array`:
    the state of the controlled system at each time step for each scenario
* `controls::Array`:
    the controls applied to the system at each time step for each scenario
"""
function forward_simulations(model::SPModel,
                            param::SdpParameters,
                            V::Union{SharedArray, Array},
                            scenarios::Array,
                            display=true::Bool)

    nb_scenarios = size(scenarios,2)

    TF = model.stageNumber
    law = model.noises
    x_dim = model.dimStates
    product_states = generate_state_grid(model, param)
    costs = SharedArray{Float64}(zeros(nb_scenarios))
    states = SharedArray{Float64}(zeros(TF,nb_scenarios,x_dim))
    controls = SharedArray{Float64}(zeros(TF-1,nb_scenarios,model.dimControls))

    dynamics = model.dynamics

    cost = build_cost_function(model.costFunctions)

    fin_cost = build_final_cost_function(model.finalCost)

    constraints = build_contraints_function(model.inequalityConstraints,
                                            model.equalityConstraints)

    args = [model.ulim, model.xlim, param.stateSteps, x_dim,
    generate_control_grid(model, param), dynamics, constraints,
    cost]

    X0 = model.initialState
    for s in 1:nb_scenarios
        states[1, s, :] = X0
    end

    best_control = tuple()

    info = param.infoStructure

    if  info == "DH"
        get_u = BellmanSolvers.exhaustive_search_dh_get_u
    elseif info == "HD"
        get_u = BellmanSolvers.exhaustive_search_hd_get_u
    else
        warn("Information structure should be DH or HD. Defaulted to DH")
        get_u = BellmanSolvers.exhaustive_search_dh_get_u
    end

    build_Ux = Nullable{Function}(param.buildSearchSpace)


    @sync @parallel for s in 1:nb_scenarios

        current_scen = scenarios[:,s,:]

        for t = 1:(TF-1)

            args_w = []

            x = states[t,s,:]
            w = current_scen[t,:]
            args_t = [value_bspline_interpolation(x_dim, V, t+1), t, x]

            if info == "DH"
                if (param.expectation_computation=="MonteCarlo")
                    sampling_size = param.monteCarloSize
                    push!(args_w,sampling_size,
                            [sampling(law,t) for i in 1:sampling_size],
                            (1./sampling_size)*ones(sampling_size))
                else
                    push!(args_w,law[t].supportSize, law[t].support, law[t].proba)
                end
            else
                push!(args_t, w)
            end

            best_control = get_u(args_w..., args..., args_t..., build_Ux)[1]

            if best_control == tuple()
                error("No u admissible")
            else
                controls[t,s,:] = [best_control...]
                states[t+1,s,:] = dynamics(t, x, best_control, w)
                costs[s] = costs[s] + cost(t, x, best_control, w)
            end
        end

    end

    for s in 1:nb_scenarios
        costs[s] = costs[s] + fin_cost(states[TF,s,:])
    end

    return costs, states, controls
end


