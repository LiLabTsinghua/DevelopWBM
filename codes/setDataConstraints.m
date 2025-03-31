function model = setDataConstraints(model,fluxes,expIDs,ecFlag,const_level,fixed)
%Biomass UB
GrowthRate = fluxes(end-1);
%value   = GrowthRate;
GrowthIndx = find(strcmpi(model.rxns,expIDs(end-1)));
%model   = setBounds(model,GrowthIndx,value,ecFlag,true);
if ecFlag
    Prot_biomass = fluxes(end); %Total protein content in biomass [g prot/g DW]
    Prot_pool    = fluxes(end); %Amount of enzymes available for biochemical reactions
    model        = rescaleBiomassProtein(model,Prot_pool,ecFlag);
end
%Glucose exchange bound
if const_level >0
    %value   = GrowthRate;
    %model   = setBounds(model,GrowthIndx,value,ecFlag,true);
    index   = strcmpi(expIDs,'MAR09034');
    rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
    value   = fluxes(index);
    model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
    %L-lactate exchange bound
    if const_level>1
        index   = strcmpi(expIDs,'MAR09135');
        rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
        value   = fluxes(index);
        model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
        %Threonine exchange bound
        if const_level>2
            index   = strcmpi(expIDs,'MAR09044');
            rxnIndx = find(strcmpi(model.rxns,expIDs(index)));
            value   = fluxes(index);
            model   = setBounds(model,rxnIndx,value,ecFlag,fixed);
        end
    end
end
end