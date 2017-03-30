module SdpCoreFunctions
using Interpolations, Gallium

export index_from_variable, real_index_from_variable

"""
Convert the state and control float tuples (stored as arrays or tuples) of the
problem into int tuples that can be used as indexes for the discretized
value functions

# Parameters
* `variable::Array`:
    the vector variable we want to convert to an index (integer)
* `bounds::Array`:
    the lower bounds for each component of the variable
* `variable_steps::Array`:
    discretization step for each component

# Returns
* `index::Tuple{Integeres}`:
    the indexes of the variable

"""
function index_from_variable(variable, bounds::Array, variable_steps::Array)
    return tuple([ 1 + floor(Int64,(1e-10+( variable[i] - bounds[i][1] )/ variable_steps[i] )) for i in 1:length(variable)]...)
end


"""
Convert the state and control float tuples (stored as arrays or tuples) of the
problem into float tuples that can be used as indexes for the interpolated
value function

# Parameters
* `variable::Array`:
    the vector variable we want to convert to an index (integer)
* `bounds::Array`:
    the lower bounds for each component of the variable
* `variable_steps::Array`:
    discretization step for each component

# Returns
* `index::Tuple{Float64}`:
    the indexes of the variable

"""
function real_index_from_variable(variable, bounds::Array, variable_steps::Array)
    return tuple([1 + ( variable[i] - bounds[i][1] )/variable_steps[i] for i in 1:length(variable)]...)
end

"""
Check if next state x_{t+1} satisfies state bounds constraints

# Parameters
* `next_stae::Array`:
    the state we want to check
* `x_dim::Int`:
    the number of state variables
* `x_bounds::Array`:
    the state variables bounds

# Returns
* `index::Tuple{Float64}`:
    the indexes of the variable

"""
function is_next_state_feasible(next_state, x_dim, x_bounds)

    next_state_box_const = true

    for i in 1:x_dim
        next_state_box_const =  (next_state_box_const&&
                                (next_state[i]>=x_bounds[i][1]-1e-10)&&
                                (next_state[i]<=x_bounds[i][2]+1e-10))
    end

    return next_state_box_const
end

"""
Computes the value function at time t evaluated at state x in a decision
hazard setting

# Parameters
* `sampling_size::int`:
    the size of the uncertainty space
* `samples::Array`:
    the uncertainties realizations
* `probas::Array`:
    the probabilities of all the uncertainties realizations
* `u_bounds::Array`:
    the control variables bounds
* `x_bounds::Array`:
    the state variables bounds
* `x_steps::Array`:
    the state variables steps
* `x_dim::int`:
    the number of state variables
* `product_controls::Array`:
    the discretized control space
* `dynamics::Function`:
    the dynamic of the problem
* `constraints::Function`:
    the constraints of the problem
* `cost::Function`:
    the cost function of the problem
* `V::Array`:
    the value functions ofthe problem
* `Vitp::Interpolations`:
    the interpolated value function at time t+1
* `t::float`:
    the time step at which the value function is computed
* `x::Array or Tuple`:
    the state at which the value function needs to be evaluated

"""
function sdp_u_w_loop(sampling_size, samples, probas, u_bounds,
                                x_bounds, x_steps, x_dim, product_controls,
                                dynamics, constraints, cost, Vitp, t, x, u_space_builder = nothing)
    expected_V = Inf
    optimal_u = tuple()
    #Loop over controls
    if typeof(u_space_builder) != Void
        controls_search_space = u_space_builder(t,x)
    else
        controls_search_space = product_controls
    end

    for u in controls_search_space

        expected_V_u = 0.
        count_admissible_w = 0

        for w = 1:sampling_size
            w_sample = samples[:, w]
            proba = probas[w]
            next_state = dynamics(t, x, u, w_sample)

            #@conditional breakpoint(constraints, Tuple{Any, Any, Any, Any}) (length(u)==2)

            if constraints(t, x, u, w_sample)&&is_next_state_feasible(next_state, x_dim, x_bounds)

                count_admissible_w = count_admissible_w + proba

                ind_next_state = real_index_from_variable(next_state, x_bounds,
                                                            x_steps)
                next_V = Vitp[ind_next_state...]
                expected_V_u += proba*(cost(t, x, u, w_sample) + next_V)

            end
        end

        if (count_admissible_w>0)

            next_V = next_V / count_admissible_w

            if (expected_V_u < expected_V)

                expected_V = expected_V_u
                optimal_u = u

            end
        end
    end

    return expected_V, optimal_u
end


"""
Computes the value function at time t evaluated at state x in a hazard
decision setting

# Parameters
* `sampling_size::int`:
    the size of the uncertainty space
* `samples::Array`:
    the uncertainties realizations
* `probas::Array`:
    the probabilities of all the uncertainties realizations
* `u_bounds::Array`:
    the control variables bounds
* `x_bounds::Array`:
    the state variables bounds
* `x_steps::Array`:
    the state variables steps
* `x_dim::int`:
    the number of state variables
* `product_controls::Array`:
    the discretized control space
* `dynamics::Function`:
    the dynamic of the problem
* `constraints::Function`:
    the constraints of the problem
* `cost::Function`:
    the cost function of the problem
* `V::Array`:
    the value functions ofthe problem
* `Vitp::Interpolations`:
    the interpolated value function at time t+1
* `t::float`:
    the time step at which the value function is computed
* `x::Array or Tuple`:
    the state at which the value function needs to be evaluated

"""
function sdp_w_u_loop(sampling_size, samples, probas, u_bounds,
                                x_bounds, x_steps, x_dim, product_controls,
                                dynamics, constraints, cost, Vitp, t, x, u_space_builder = nothing)

    expected_V = 0.
    count_admissible_w = 0.

    #Compute expectation
    for w in 1:sampling_size
        admissible_u_w_count = 0
        w_sample = samples[:, w]
        proba = probas[w]

        uopt, best_V_x_w, admissible_u_w_count = sdp_hd_get_u(sampling_size, samples,
                                                            probas, u_bounds,
                                                            x_bounds, x_steps,
                                                            x_dim,
                                                            product_controls,
                                                            dynamics, constraints,
                                                            cost, Vitp, t, x, w_sample,
                                                            u_space_builder)


        expected_V += proba*best_V_x_w
        count_admissible_w += (admissible_u_w_count>0)*proba
    end
    if (count_admissible_w>0.)
        expected_V = expected_V / count_admissible_w
    end

    return expected_V
end

function sdp_hd_get_u(sampling_size, samples, probas, u_bounds,
                    x_bounds, x_steps, x_dim, product_controls,
                    dynamics, constraints, cost, Vitp, t, x, w,
                    u_space_builder = nothing)

    if typeof(u_space_builder) == Function
        controls_search_space = u_space_builder(t,x,w)
    else
        controls_search_space = product_controls
    end

    best_V_x_w = Inf
    next_V_x_w = Inf
    optimal_u = tuple()
    admissible_u_w_count = 0

    for u in product_controls

        next_state = dynamics(t, x, u, w)

        if constraints(t, x, u, w)&&is_next_state_feasible(next_state, x_dim, x_bounds)
            admissible_u_w_count  += 1
            ind_next_state = real_index_from_variable(next_state, x_bounds,
                                                        x_steps)
            next_V_x_w = cost(t, x, u, w) + Vitp[ind_next_state...]

            if (next_V_x_w < best_V_x_w)
                best_V_x_w = next_V_x_w
                optimal_u = u
            end

        end
    end

    return optimal_u, best_V_x_w, admissible_u_w_count
end

function sdp_dh_get_u(sampling_size, samples, probas, u_bounds,
                        x_bounds, x_steps, x_dim, product_controls,
                        dynamics, constraints, cost, Vitp, t, x, w, u_space_builder = nothing)

    return (sdp_u_w_loop(sampling_size, samples, probas, u_bounds,
                        x_bounds, x_steps, x_dim, product_controls,
                        dynamics, constraints, cost, Vitp, t, x, u_space_builder)[2])

end

end