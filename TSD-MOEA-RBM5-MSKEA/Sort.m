function [SubPopulation,Rank] = Sort(ns,SubPopulation)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Rank = cell(1,ns);
    for i = 1 : ns
        if ~isempty(SubPopulation{i})
            % [~,FrontNo,CrowdDis,~,~] = EnviSelect(SubPopulation{i},length(SubPopulation{i}),MultiNonzeros);
            [~,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection(SubPopulation{i},length(SubPopulation{i}));
            % [~,FrontNo,CrowdDis] = EnvironmentalSelection(SubPopulation{i},length(SubPopulation{i}));
            [~,rank]             = sortrows([FrontNo',-CrowdDis']);
            SubPopulation{i}     = SubPopulation{i}(rank);
            Rank{i}              = 1 : length(rank);
        end
    end
end