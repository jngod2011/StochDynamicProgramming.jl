using Base.Meta

function SPModel(str_type)
    if str_type == "Linear"
        return(LinearDynamicLinearCostSPmodel())
    elseif str_type == "PiecewiseLinear"
        return(PiecewiseLinearCostSPmodel())
    else
        return(StochDynProgModel())
    end
end

macro addState(model, arg...)
    state_name = string(arg[1].args[3].args[1])
    n_stages = arg[1].args[3].args[2].args[2]
    quote
        push!($(esc(model)).state_names, $state_name);
        push!($(esc(model)).dynamic_expressions, "empty")
        push!($(esc(model)).xlim, ($(arg[1].args[1]),$(arg[1].args[5])));
        $(esc(model)).dimStates += 1;
        $(esc(model)).stageNumber = $(esc(n_stages));
    end
end

macro addControl(model, arg...)
    control_name = string(arg[1].args[3].args[1])
    quote
        push!($(esc(model)).control_names, $control_name);
        push!($(esc(model)).ulim, ($(arg[1].args[1]),$(arg[1].args[5])));
        $(esc(model)).dimControls += 1;
    end
end

macro addNoise(model, arg...)
    #TODO : Multidimensional case!!!
    aleas_name = string(arg[1])
    alea = arg[1]
    quote
        push!($(esc(model)).aleas_names, $(aleas_name));
        $(esc(model)).noises = NoiseLaw[NoiseLaw($(esc(alea))[t][2,:], $(esc(alea))[t][1,:]) for t in 1:$(esc(model)).stageNumber-1];
        $(esc(model)).dimNoises += 1;
    end
end

macro setStochObjective(model, ope, arg...)
    #TODO : almost everything, Max, finalCostFunction, prevent case where a variable name is at the end of another one (regExp?)...
    #Better way to handle this part?
    if ope == :Min
        ind = find(arg[1].args .== :sum)[1]
        costExpression = string(arg[1].args[ind+1])
        esc(quote
                costString = $(costExpression)
                for i in length($(model).state_names)
                    costString = replace(costString, string($(model).state_names[i],"[i]"),string("x"string([i])))
                end
                for i in length($(model).control_names)
                    costString = replace(costString, string($(model).control_names[i],"[i]"),string("u"string([i])))
                end
                for i in length($(model).aleas_names)
                    costString = replace(costString, string($(model).aleas_names[i],"[i]"),string("w"string([i])))
                end
                expr =  quote
                            function cost_t(i,x,u,w)
                                return($(parse(costString)))
                            end
                        end
                $(model).costFunctions = eval(expr)
                expr = 0;
                costString = 0;
            end)

    elseif ope == :Max

    else
        #@error("Missing min or max")
    end
end

macro addDynamic(model, arg...)
        #TODO : test in the multidimensional case, test execution time
        state_name = string(arg[1].args[1].args[1])
        dynamicExpression = string(arg[1].args[2])

            esc(quote
                    ind_dyn = find($(model).state_names .== string($(state_name)))[1]
                    dynamicString = $(dynamicExpression)
                    for i in length($(model).state_names)
                        dynamicString = replace(dynamicString, string($(model).state_names[i],"[i]"),string("x"string([i])))
                    end
                    for i in length($(model).control_names)
                        dynamicString = replace(dynamicString, string($(model).control_names[i],"[i]"),string("u"string([i])))
                    end
                    for i in length($(model).aleas_names)
                        dynamicString = replace(dynamicString, string($(model).aleas_names[i],"[i]"),string("w"string([i])))
                    end
                    dynamicString = parse(dynamicString)

                    $(model).dynamic_expressions[ind_dyn] = dynamicString

                    num_undefined_func = length(find($(model).dynamic_expressions .== "empty" ))

                    if num_undefined_func == 0

                        expr_tab = $(model).dynamic_expressions

                        expr =  quote
                                    function dyn_t(i,x,u,w)
                                        return([$(expr_tab...)])
                                    end
                                end
                        $(model).dynamics = eval(expr)
                    end

                    ind_dyn = 0;
                    dynamicString = 0;
                    num_undefined_func = 0;
                    expr = 0;
                end)
end

macro addConstraintsdp(model, arg...)
        #TODO : Change bad name (conflict with JuMP)
        if (match(r"\[1]",string(arg[1])) == nothing)
            1+1
        else
            esc(quote
                    $(model).initialState = [$(arg[1].args[3])]
                end
                )
        end
end

macro addCost(model, arg...)
    #For the piecewise case and maybe to replace setStochObjective or allow the user to choose
end

macro setRiskMeasure(model, arg...)
    #For the piecewise case and maybe to replace setStochObjective or allow the user to choose
end

function solveInterface(model::SPModel, strSDP, arg...)
    #Currently works in one particular case exclusively
    #Bad name but conflicts with the solve function of JuMP...
    if strSDP== "SDP"
        if arg[1]=="HD"
            paramSDP = SDPparameters(model, arg[2], arg[2], arg[1])
            Vs = sdp_optimize(model,paramSDP)
            lb_sdp = get_value(model,paramSDP,Vs)
            println("Value obtained by SDP: "*string(lb_sdp))
        end
    elseif strSDP=="SDDP"
        paramSDDP = SDDPparameters(arg...) # 10 forward pass, stop at MAX_ITER
        V, pbs = solve_SDDP(model, paramSDDP, 10) # display information every 10 iterations
        lb_sddp = StochDynamicProgramming.get_lower_bound(model, paramSDDP, V)
        println("Lower bound obtained by SDDP: "*string(lb_sddp))
    end
end
