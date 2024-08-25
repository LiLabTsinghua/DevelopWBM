function cmodels = RemoveSexOrgans(models,Sex)
% This function performs to remove sexual organ GEMs from models.
%
% Input:
%   models        a struct of all organ-specific GEMs
%   Sex           "Male" or "Female"
%
% Output:
%   cmodels       a struct wuithout sexual organ GEMs

if nargin ~= 2 
    EM = 'Only need a struct of GEMs and gender info in this function!';
    dispEM(EM);
end

sexorgan = zeros(length(models), 1);

if isequal(Sex, 'Male')
    removeOrgans = ["Prostate", "Testis"];
    for i=1:length(models)
        model = models{i,1};
        if contains(model.id, removeOrgans)
           sexorgan(i,1) = 1;
        end
    end
elseif isequal(Sex, 'Female')
    removeOrgans = ["Breast", "Cervix", "Ovary", "Uterus"];
    for i=1:length(models)
        model = models{i,1};
        if contains(model.id, removeOrgans)
           sexorgan(i,1) = 1;
        end
    end
end

sexorganIndex = find(sexorgan);
models(sexorganIndex) = [];
cmodels = models;

end
