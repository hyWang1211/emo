function [Mask, GroupIdx] = HelpSubPopulation(MultiNonzeros, Fitness, Mask, GroupIdx, GroupSize, ns, Problem)
    MaxNonzeros = any(MultiNonzeros, 1);
    MinNonzeros = all(MultiNonzeros, 1);

    %% Maximum and Minimum Initialization
    Mask((ns - 3 - 1)*GroupSize + 1: (ns - 3)*GroupSize, :) = repmat(MaxNonzeros,GroupSize,1);
    GroupIdx{ns - 3} = (ns - 3 - 1)*GroupSize + 1: (ns - 3)* GroupSize;

    Mask((ns - 2 - 1)*GroupSize + 1: (ns - 2)*GroupSize, :) = repmat(MinNonzeros,GroupSize,1);
    GroupIdx{ns - 2} = (ns - 2 - 1)*GroupSize + 1: (ns - 2)* GroupSize;
    %% Uniform Initialization
    Mask((ns - 1 - 1)*GroupSize + 1: (ns - 1)*GroupSize, :) = UniformPoint(GroupSize,Problem.D,'Latin') > 0.5;
    GroupIdx{ns - 1} = (ns - 1 - 1)*GroupSize + 1: (ns - 1)* GroupSize;
    %% Use SparseEA Initialization
    for i = (ns - 1)*GroupSize + 1 : Problem.N
        Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
    end
    GroupIdx{ns} = (ns - 1)*GroupSize + 1: Problem.N;

    %% Only Maximum
    % Mask((ns - 1)*GroupSize + 1: ns*GroupSize, :) = repmat(MaxNonzeros,GroupSize,1);
    % GroupIdx{ns} = (ns - 1)*GroupSize + 1: ns* GroupSize;

    %% Only Minimum
    % Mask((ns - 1)*GroupSize + 1: ns*GroupSize, :) = repmat(MinNonzeros,GroupSize,1);
    % GroupIdx{ns} = (ns - 1)*GroupSize + 1: ns* GroupSize;

    %% Only Uniform
    % Mask((ns - 1)*GroupSize + 1: ns*GroupSize, :) = UniformPoint(GroupSize,Problem.D,'Latin') > 0.5;
    % GroupIdx{ns} = (ns - 1)*GroupSize + 1: ns* GroupSize;

    %% Only SparseEA
    % for i = (ns - 1)*GroupSize + 1 : Problem.N
    %     Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
    % end
    % GroupIdx{ns} = (ns - 1)*GroupSize + 1: Problem.N;

    %% Use Uniform and SparseEA
    % Mask((ns - 1 - 1)*GroupSize + 1: (ns - 1)*GroupSize, :) = UniformPoint(GroupSize,Problem.D,'Latin') > 0.5;
    % GroupIdx{ns} = (ns - 1 - 1)*GroupSize + 1: (ns - 1)* GroupSize;
    % 
    % for i = (ns - 1)*GroupSize + 1 : Problem.N
    %     Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
    % end
    % GroupIdx{ns} = (ns - 1)*GroupSize + 1: Problem.N;
end

%
% 1. 算子方向向量融合
% 2. 迁移的算子和判断
% 3. 种群的删减
%
