function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,Nonzeros)
% The operator of SparseEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [N,D]       = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    
    OffMask = [];
    for i = 1 : N/2
        OffMask1                    = Parent1Mask(i, :);
        OffMask2                    = Parent2Mask(i, :);
        innerposition               = (Nonzeros == OffMask1) & (Nonzeros == OffMask2);

        % 
        % Crossover
        k     = rand(1,D) < 0.5;
        OffMask1(k) = Parent2Mask(i, k);
        OffMask2(k) = Parent1Mask(i, k);
        % if rand < 0.5
        %     % execute Uniform crossover
        %     k                           = innerposition;
        %     % k(repmat(rand(1,1)>1,1,D))  = false; default crossover
        %     % probability is 1
        % 
        %     OffMask1(k)                 = Parent2Mask(i, k);
        %     OffMask2(k)                 = Parent1Mask(i, k);
        % else
        %     % execute Uniform crossover in other position
        %     k                           = ~innerposition;
        %     OffMask1(k)                 = Parent2Mask(i, k);
        %     OffMask2(k)                 = Parent1Mask(i, k);
        % end

        %% Mutation
        mutate_pos = randi(Problem.D, 1, 2);
        OffMask1(mutate_pos(1)) = ~OffMask1(mutate_pos(1));
        OffMask2(mutate_pos(2)) = ~OffMask2(mutate_pos(2));

        OffMask = [OffMask; OffMask1;OffMask2];
        
    end
    
    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        OffDec = OperatorGA(Problem,ParentDec);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(N/2,D);
    end
end