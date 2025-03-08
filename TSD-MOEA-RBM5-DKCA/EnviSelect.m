function [SubPopulation,FrontNo,CrowdDis,ns,MultiNonzeros] = EnviSelect(Population,N,MultiNonzeros)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicated solutions
    % [~,uni] = unique(Population.objs,'rows');
    % Population = Population(uni);
    % N          = min(N,length(Population));

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;

    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;

    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);

    %% Divide Population into subPopulation
    [task_num, ~]        = unique(Population.adds([], 1));
    ns                   = size(MultiNonzeros, 1);
    SubPopulation        = cell(1, ns);
    j = 1;
    for i = 1: ns
        if ismember(i, task_num)
            try
                SubPopulation{i} = Population(Population.adds([], 1) == task_num(j));
                j = j + 1;
            catch err
                disp(err);
            end
        else
            SubPopulation{i} = [];
        end
    end
end