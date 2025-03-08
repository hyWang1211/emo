function [Population,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection(Population,N)
% The environmental selection of MSKEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lei Chen

    %% Delete duplicated solutions
    [~,uni] = unique(Population.objs,'rows');
    Population = Population(uni);
    N          = min(N,length(Population));
    % Calculate the fitness of each solution
    
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = false(1,length(FrontNo));
    Next(FrontNo<MaxFNo) = true;
    
    PopObj = Population.objs;
    fmax   = max(PopObj(FrontNo==1,:),[],1);
    fmin   = min(PopObj(FrontNo==1,:),[],1);
    PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);

    %% Environmental selection
    Last = find(FrontNo==MaxFNo);
    del  = Truncation(PopObj(Last,:),length(Last)-N+sum(Next));
    Next(Last(~del)) = true;
    % Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    %% Divide Population into subPopulation (original subpopulation update)
    
    % ns = size(MultiNonzeros, 1);
    % [element, ~] = unique(Population.adds([], 1));
    % ns = size(element, 1);
    % SubPopulation = cell(1, ns);
    % for i = 1: ns
    %     SubPopulation{i} = Population(Population.adds([], 1) == element(i));
    % end
    % MultiNonzeros = MultiNonzeros(element, :);
    % [task_num, ~]        = unique(Population.adds([], 1));
    % ns                   = size(MultiNonzeros, 1);
    % SubPopulation        = cell(1, ns);
    % j = 1;
    % for i = 1: ns
    %     if ismember(i, task_num)
    %         try
    %             SubPopulation{i} = Population(Population.adds([], 1) == task_num(j));
    %             j = j + 1;
    %         catch err
    %             disp(err);
    %         end
    %     else
    %         SubPopulation{i} = [];
    %     end
    % end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end