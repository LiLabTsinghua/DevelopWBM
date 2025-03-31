function [JI, Intersection, Union] = JaccardIndex(model1,model2, includeGene)
% This function performs to get JaccardIndex to evaluate similarity of two
% models.
%
% Input:
%   model1/2        cell/organ-specific metabolic models
%
% Output:
%   JI              a value for JaccardIndex
%   intersection    Intersection of reactions in the two models
%   union           Union of reactions in the two models

if nargin < 3 || isempty(includeGene)
    includeGene = true;
end

if nargin < 2 
    EM = 'Only two models needed£¡';
    dispEM(EM);
end

if ~(isfield(model1, 'rxns') && isfield(model2, 'rxns'))
    disp('Something wrong with the models, which may not contain rxns');
end

if includeGene
    M1_index = ~cellfun('isempty', model1.grRules);
    M1 = model1.rxns(M1_index);
    M2_index = ~cellfun('isempty', model2.grRules);
    M2 = model2.rxns(M2_index);
else
    M1 = model1.rxns;
    M2 = model2.rxns;
end


Intersection = intersect(M1, M2);
Union = union(M1, M2);

JI = numel(Intersection)/numel(Union);

end



