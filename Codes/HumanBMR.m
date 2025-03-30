function Energy_kcal = HumanBMR(model, solution)

    S = model.S;
    F = max(-S,0);
    R = max(S,0);
    vf = max(solution.x,0);
    vr = max(-solution.x,0);
    production=[R F]*[vf; vr];
    consumption=[F R]*[vf; vr];
    % find all reactions in the model that involve atp
    atp = (find(~cellfun(@isempty,strfind(model.mets,'MAM01371'))));
    % sum of atp consumption in the flux distribution
    Sum_atp=sum(consumption(atp,1)); % (in mmol ATP/day/person)
    % compute the energy release in kJ
    Energy_kJ = Sum_atp/1000 * 64; % (in kJ/day/person)
    % compute the energy release in kcal, where 1 kJ = 0.239006 kcal
    % Energy_kcal = Energy_kJ*0.239006; % (in kcal/day/person)
    Energy_kcal = Energy_kJ*0.239006*2; % (in kcal/day/person)
end