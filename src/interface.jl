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

    ind = find(arg[1].args .== :sum)[1]
    costExpression = string(arg[1].args[ind+1])
    if ope == :Min
        esc(quote
                $(model).costFunctions =    @generated  function cost_t(i,x,u,w)
                                        return(StochDynamicProgramming.generateStochObjective($(model), $(costExpression)))
                                    end
            end)

    elseif ope == :Max

    else
        #@error("Missing min or max")
    end
end

function generateStochObjective(model, costExpr)
    costString = costExpr
    for i in 1:length(model.state_names)
        costString = replace(costString, string(model.state_names[i],"[i]"),string("x"string([i])))
    end
    for i in 1:length(model.control_names)
        costString = replace(costString, string(model.control_names[i],"[i]"),string("u"string([i])))
    end
    for i in 1:length(model.aleas_names)
        costString = replace(costString, string(model.aleas_names[i],"[i]"),string("w"string([i])))
    end

    return(parse(costString))
end

macro addDynamic(model, arg...)
        #TODO : test in the multidimensional case, test execution time
        state_name = string(arg[1].args[1].args[1])
        dynamicExpression = string(arg[1].args[2])

        esc(quote
                tabExpr = StochDynamicProgramming.generateDynamicEquation($(model), string($(state_name)), $(dynamicExpression))
                $(model).dynamics = @generated  function dyn_t(i,x,u,w)
                                        return(tabExpr)
                                    end
        end)
end

function generateDynamicEquation(model, dyn_state, dynamicExpression)

    ind_dyn = find(model.state_names .== dyn_state)[1]
    dynamicString = dynamicExpression
    for i in 1:length(model.state_names)
        dynamicString = replace(dynamicString, string(model.state_names[i],"[i]"),string("x"string([i])))
    end
    for i in 1:length(model.control_names)
        dynamicString = replace(dynamicString, string(model.control_names[i],"[i]"),string("u"string([i])))
    end
    for i in 1:length(model.aleas_names)
        dynamicString = replace(dynamicString, string(model.aleas_names[i],"[i]"),string("w"string([i])))
    end

    model.dynamic_expressions[ind_dyn] = dynamicString

    tabDynStr = "["

    tabDynStr = string(tabDynStr,model.dynamic_expressions[1])
    if (length(model.dynamic_expressions)>1)&(length(find(model.dynamic_expressions .== "empty"))==0)
        for str in model.dynamic_expressions[2:end]
            tabDynStr = string(tabDynStr,",",str)
        end
    end

    tabDynStr = string(tabDynStr,"]")

    return(parse(tabDynStr))

end

macro addDPConstraint(model, arg...)
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
