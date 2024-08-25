function [JI, Intersection, Union] = JaccardIndex(model1,model2)
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

if nargin ~= 2 
    EM = 'Only two models needed£¡';
    dispEM(EM);
end

if ~(isfield(model1, 'rxns') && isfield(model2, 'rxns'))
    disp('Something wrong with the models, which may not contain rxns');
end

Intersection = intersect(model1.rxns, model2.rxns);
Union = union(model1.rxns, model2.rxns);

JI = numel(Intersection)/numel(Union);

end



