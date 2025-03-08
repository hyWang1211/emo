function SubOffspring = ImprovedCreateOff1(Problem,Population,rmp,MultiNonzeros)
% function [SubPopulation, MultiNonzeros] = ImprovedCreateOff1(Problem,Population,rmp,MultiNonzeros)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parent selection
    Parent11 = [];
    Parent21 = [];
    Parent12 = [];
    Parent22 = [];
    for i = 1 : floor(length(Population)/2)
        P1 = Population(i);
        P2 = Population(i+floor(length(Population)/2));
        if (P1.add{1} == P2.add{1}) || (rand<rmp)
            Parent11 = [Parent11,P1];
            Parent21 = [Parent21,P2];
        else
            Parent12 = [Parent12,P1];
            Parent22 = [Parent22,P2];
        end
    end

    %% Offspring generation
    if ~isempty(Parent11)
        Parent11Dec  = Parent11.adds([], 3, 0);
        Parent21Dec  = Parent21.adds([], 3, 0);

        % Parent11Mask = Parent11.adds([], 2);
        % Parent21Mask = Parent21.adds([], 2);

        % For dec
        try
            OffDec1      = GAreal(Parent11Dec,Parent21Dec,Problem.lower,Problem.upper,1,20,1,20);
        catch err
            disp(err);
        end
        
        % For mask
        OffMask1     = GAbinary(Parent11,Parent21,MultiNonzeros);
        cutSzie = floor(size(Parent11, 2) / 2);
        TaksNumberForOff1 = [Parent11(1:cutSzie).adds([], 1, 0); Parent21(cutSzie+1:end).adds([], 1, 0)];
        % TaksNumberForOff1 = [Parent11.adds([], 1, 0); Parent21.adds([], 1, 0)];
    else
        OffDec1 = [];
        OffMask1 = [];
        TaksNumberForOff1 = [];
    end

    if ~isempty(Parent12)
        Parent12Dec     = Parent12.adds([], 3, 0);
        Parent22Dec     = Parent22.adds([], 3, 0);

        % Parent12Mask     = Parent12.adds([], 2);
        % Parent22Mask     = Parent22.adds([], 2);

        % For dec
        OffDec2     = GAreal(Parent12Dec,Parent22Dec,Problem.lower,Problem.upper,0,20,1,20);
        % For mask
        OffMask2    = GAbinary(Parent12,Parent22,MultiNonzeros);

        cutSzie = floor(size(Parent12, 2) / 2);
        TaksNumberForOff2 = [Parent12(1:cutSzie).adds([], 1, 0); Parent22(cutSzie+1:end).adds([], 1, 0)];
        % TaksNumberForOff2 = [Parent12.adds([], 1, 0); Parent22.adds([], 1, 0)];
    else
        OffDec2 = [];
        OffMask2 = [];
        TaksNumberForOff2 = [];
    end

    SubOffspring = Problem.Evaluation([OffDec1.*OffMask1;OffDec2.*OffMask2]);
    SubOffspring.adds([TaksNumberForOff1;TaksNumberForOff2], 1, 0);
    SubOffspring.adds([OffMask1;OffMask2],2, 0);
    SubOffspring.adds([OffDec1;OffDec2],3, 0);
    
    % original subpopulation update
    % [task_num, ~]          = unique(SubOffspring.adds([], 1));
    % ns                     = size(MultiNonzeros, 1);
    % SubPopulation          = cell(1, ns);
    % j = 1;
    % for i = 1 : ns
    %     if ismember(i, task_num)
    %         try
    %             SubPopulation{i} = SubOffspring((SubOffspring.adds([], 1) == task_num(j)));
    %             j = j + 1;
    %         catch err
    %             disp(err);
    %         end
    %     else
    %         SubPopulation{i} = [];
    %     end
    %     % SubPopulation{i} = SubOffspring((SubOffspring.adds([], 1) == task_num(i)));
    %     % TempMultiNonzeros(i) = MultiNonzeros(task_num(i));
    % end
    % MultiNonzeros = TempMultiNonzeros;
end

% original
% function Offspring = GAbinary(Parent1,Parent2,MultiNonzeros)
% 
%     Parent1Mask = Parent1.adds([], 2, 0);
%     Parent2Mask = Parent2.adds([], 2, 0);
%     Parent1GuideVector = Parent1.adds([], 1, 0);
%     Parent2GuideVector = Parent2.adds([], 1, 0);
%     [N,D]       = size(Parent1Mask);
% 
%     %% Crossover for mask
%     OffMask = Parent1Mask;
%     for i = 1 : N
%         P1 = Parent1Mask(i, :);
%         P2 = Parent2Mask(i, :);
%         try
%             G1 = MultiNonzeros(Parent1GuideVector(i, :)', :);
%             G2 = MultiNonzeros(Parent2GuideVector(i, :)', :);
%         catch err
%             disp(err);
%         end
% 
%         if rand < 0.5
%             % Use maximum crossover
%             index1 = xor(P1, P2);
%             index2 = G1 | G2;
% 
%             index3 = find_two_nonzero_positions(index1);
% 
%             if length(index3) == 1
%                 if index2(index3) == 1
%                     OffMask(i, index3) = 1;
%                 else
%                     OffMask(i, index3) = 0;
%                 end
%             elseif length(index3) == 2
%                 if all(index2(index3) == [1 1])
%                     OffMask(i, index3) = 1;
%                 elseif all(index2(index3) == [0 0])
%                     OffMask(i, index3) = 0;
%                 else
%                     OffMask(i, index3(1)) = 1;
%                     OffMask(i, index3(2)) = 0;
%                 end
%             elseif isempty(index3)
%                 %% Uniform crossover
%                 [N,D] = size(P1);
%                 k     = rand(N,D) < 0.5;
%                 k(repmat(rand(N,1)>1,1,D)) = false;
%                 OffMask(i, k) = P2(k);
%             end
%         else
%             % Use minimum crossover
%             index1 = xor(P1, P2);
%             index2 = G1 & G2;
% 
%             index3 = find_two_nonzero_positions(index1);
% 
%             if length(index3) == 1
%                 if index2(index3) == 1
%                     OffMask(i, index3) = 1;
%                 else
%                     OffMask(i, index3) = 0;
%                 end
%             elseif length(index3) == 2
%                 if all(index2(index3) == [1 1])
%                     OffMask(i, index3) = 1;
%                 elseif all(index2(index3) == [0 0])
%                     OffMask(i, index3) = 0;
%                 else
%                     OffMask(i, index3(1)) = 0;
%                     OffMask(i, index3(2)) = 1;
%                 end
%             elseif isempty(index3)
%                 %% Uniform crossover
%                 [N,D] = size(P1);
%                 k     = rand(N,D) < 0.5;
%                 k(repmat(rand(N,1)>1,1,D)) = false;
%                 OffMask(i, k) = P2(k);
%             end
%         end
%     end
% 
%     %% Mutation for mask
%     Site = rand(N,D) < 1/D;
%     Offspring = OffMask;
%     try
%         Offspring(Site) = ~Offspring(Site);
%     catch err
%         disp(err)
%     end
% end

% Use maximum only
function Offspring = GAbinary(Parent1,Parent2,MultiNonzeros)

    Parent1Mask = Parent1.adds([], 2, 0);
    Parent2Mask = Parent2.adds([], 2, 0);
    Parent1GuideVector = Parent1.adds([], 1, 0);
    Parent2GuideVector = Parent2.adds([], 1, 0);
    [N,D]       = size(Parent1Mask);

    %% Crossover for mask
    OffMask = Parent1Mask;
    for i = 1 : N
        P1 = Parent1Mask(i, :);
        P2 = Parent2Mask(i, :);
        % try
        %     G1 = MultiNonzeros(Parent1GuideVector(i, :)', :);
        %     G2 = MultiNonzeros(Parent2GuideVector(i, :)', :);
        % catch err
        %     disp(err);
        % end

        % Use maximum crossover
        index1 = xor(P1, P2);
        % index2 = G1 | G2;

        index3 = find_one_nonzero_positions(index1);

        % if length(index3) == 1
        %     if G2(index3) == G1(index3)
        %         OffMask(i, index3) = G1(index3);
        %     else
        %         if rand < 0.5
        %             OffMask(i, index3) = G1(index3);
        %         else
        %             OffMask(i, index3) = G2(index3);
        %         end
        %     end
        % 
        % elseif length(index3) == 2
        %     if all(index2(index3) == [1 1])
        %         OffMask(i, index3) = 1;
        %     elseif all(index2(index3) == [0 0])
        %         OffMask(i, index3) = 0;
        %     else
        %         OffMask(i, index3(1)) = 1;
        %         OffMask(i, index3(2)) = 0;
        %     end
        % elseif isempty(index3)
            %% Uniform crossover
            [N,D] = size(P1);
            k     = rand(N,D) < 0.5;
            k(repmat(rand(N,1)>1,1,D)) = false;
            OffMask(i, k) = P2(k);
        % end
    end

    %% Mutation for mask
    Site = rand(N,D) < 1/D;
    Offspring = OffMask;
    try
        Offspring(Site) = ~Offspring(Site);
    catch err
        disp(err)
    end
end

% function Offspring = GAbinary(Parent1,Parent2,MultiNonzeros)
% 
%     Parent1Mask = Parent1.adds([], 2, 0);
%     Parent2Mask = Parent2.adds([], 2, 0);
%     Parent1GuideVector = Parent1.adds([], 1, 0);
%     Parent2GuideVector = Parent2.adds([], 1, 0);
%     [N,D]       = size(Parent1Mask);
% 
%     %% Crossover for mask
%     OffMask = Parent1Mask;
%     for i = 1 : N
%         P1 = Parent1Mask(i, :);
%         P2 = Parent2Mask(i, :);
%         try
%             G1 = MultiNonzeros(Parent1GuideVector(i, :)', :);
%             G2 = MultiNonzeros(Parent2GuideVector(i, :)', :);
%         catch err
%             disp(err);
%         end
% 
%         % Use maximum crossover
%         % index1 = xor(P1, P2);
%         if rand < 0.5
%             index1 = find(P1 & ~P2);
%             index2 = find_one_nonzero_positions(index1);
% 
%             if ~isempty(index2)
%                 if G1(index2) == 0
%                     % flip to 0
%                     OffMask(i, index2) = 0;
%                 else
%                     OffMask(i, index2) = 1;
%                 end
%             else
%                 [N,D] = size(P1);
%                 k     = rand(N,D) < 0.5;
%                 k(repmat(rand(N,1)>1,1,D)) = false;
%                 OffMask(i, k) = P2(k);
%             end
%         else
%             index1 = find(~P1 & P2);
%             index2 = find_one_nonzero_positions(index1);
% 
%             if ~isempty(index2)
%                 if G2(index2) == 0
%                     %flip to 0
%                     OffMask(i, index2) = 0;
%                 else
%                     OffMask(i, index2) = 1;
%                 end
%             else
%                 [N,D] = size(P1);
%                 k     = rand(N,D) < 0.5;
%                 k(repmat(rand(N,1)>1,1,D)) = false;
%                 OffMask(i, k) = P2(k);
%             end
%         end
%     end
% 
%     %% Mutation for mask
%     Site = rand(N,D) < 1/D;
%     Offspring = OffMask;
%     try
%         Offspring(Site) = ~Offspring(Site);
%     catch err
%         disp(err)
%     end
% end

function Offspring = GAreal(Parent1,Parent2,lower,upper,proC,disC,proM,disM)
% Genetic operators for real and integer variables

    %% Simulated binary crossover
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;

    %% Polynomial mutation
    Lower = repmat(lower,N,1);
    Upper = repmat(upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end

function result = find_one_nonzero_positions(index1)
    nonzero_positions = find(index1);
    if isempty(nonzero_positions)
        result = [];
        return;
    end
    
    if length(nonzero_positions) == 1
        result = nonzero_positions;
        return;
    end
    result = nonzero_positions(randperm(length(nonzero_positions), 1));
end

function result = find_two_nonzero_positions(index1)
    nonzero_positions = find(index1);
    if isempty(nonzero_positions)
        result = [];
        return;
    end
    
    if length(nonzero_positions) == 1
        result = [nonzero_positions, nonzero_positions];
        return;
    end
    result = nonzero_positions(randperm(length(nonzero_positions), 2));
end