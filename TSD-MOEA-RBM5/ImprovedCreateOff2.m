function SubOffspring = ImprovedCreateOff2(Problem,Population,rmp,MultiNonzeros, ntp, Fitness)
    [N,D]       = size(Population.adds([], 2, 0));
    Parent1 = Population(1:N/2);
    Parent2 = Population(N/2+1:end);

    NormalPop = [];
    MFEAPop1 = [];
    MFEAPop2 = [];
    SparseEAPop = [];
    % cnt1 = 0;
    % cnt2 = 0;
    % cnt3 = 0;
    % assign different variation strategy
    for i = 1: floor(length(Population)/2)
        P1 = Parent1(i);
        P2 = Parent2(i);

        skill_p1 = P1.add{1};
        skill_p2 = P2.add{1};
        
        if (skill_p1 <= ntp) && (skill_p2 <= ntp)
            %% Use MFEA framework and custom variation operator.
            if (skill_p1 == skill_p2) || (rand < rmp)
                MFEAPop1 = [MFEAPop1, P1, P2];
            else
                MFEAPop2 = [MFEAPop2, P1, P2];
            end
        % elseif (skill_p1 > 6) && (skill_p2 > 6)
        %     %% Use the operator in SparseEA and randomly select a parent's skill to inherit.
        %     SparseEAPop = [SparseEAPop, P1, P2];
        else
            %% Use the normal operator for sparse problems.
            % NormalPop = [NormalPop, P1, P2];
            SparseEAPop = [SparseEAPop, P1, P2];
        end
    end

    FinalMask = [];
    FinalDec = [];
    TaskNumber = [];

    %% MFEA
    if ~isempty(MFEAPop1)
        %% generate MFEA population, case 1.
        Parent1 = MFEAPop1(1:2:end);
        Parent2 = MFEAPop1(2:2:end);

        %% Crossover and mutation for mask, case 2.
        SubOffMask = GAbinaryCustom(Parent1, Parent2, MultiNonzeros);

        %% Crossover and mutation for dec, case 3.
        if any(Problem.encoding~=4)
            SubOffDec = GAreal(Parent1.adds([], 3, 0),Parent2.adds([], 3, 0),Problem.lower,Problem.upper,1,20,1,20);
        else
            SubOffDec = ones(length(Parent1),D);
        end

        
        FinalMask = [FinalMask; SubOffMask];
        FinalDec = [FinalDec; SubOffDec];

        cutSzie = floor(size(Parent1, 2) / 2);
        TaskNumber = [TaskNumber;Parent1(1:cutSzie).adds([], 1, 0); Parent2(cutSzie+1:end).adds([], 1, 0)];
    end

    if ~isempty(MFEAPop2)
        
        %% generate MFEA population, case 2. (There is no crossover)
        Parent1 = MFEAPop2(1:2:end);
        Parent2 = MFEAPop2(2:2:end);
        
        %% mutation for mask
        SubOffMask = GAbinary(Parent1.adds([], 2, 0), Parent2.adds([], 2, 0), 0, 1);
        
        %% mutation for dec
        if any(Problem.encoding~=4)
            SubOffDec = GAreal(Parent1.adds([], 3, 0),Parent2.adds([], 3, 0),Problem.lower,Problem.upper,1,20,1,20);
        else
            SubOffDec = ones(length(Parent1),D);
        end

        
        FinalMask = [FinalMask; SubOffMask];
        FinalDec = [FinalDec; SubOffDec];

        % TaskNumber = [TaskNumber;Parent1(1:2:end).adds([], 1, 0); Parent2(2:2:end).adds([], 1, 0)];
        cutSzie = floor(size(Parent1, 2) / 2);
        TaskNumber = [TaskNumber;Parent1(1:cutSzie).adds([], 1, 0); Parent2(cutSzie+1:end).adds([], 1, 0)];
    end

    %% SparseEA2
    if ~isempty(SparseEAPop)
        %% Crossover for mask
        Parent1 = SparseEAPop(1:2:end);
        Parent2 = SparseEAPop(2:2:end);
        Parent1Mask = Parent1.adds([], 2, 0);
        Parent2Mask = Parent2.adds([], 2, 0);
        Parent1Dec = Parent1.adds([], 3, 0);
        Parent2Dec = Parent2.adds([], 3, 0);
        
        %% Crossover and mutation for dec
         if any(Problem.encoding~=4)
             [SubOffDec,groupIndex,chosengroups] = GLP_OperatorGAhalf(Problem,Parent1Dec,Parent2Dec,4);	% 4 -- numberofgroups
             SubOffDec(:,Problem.encoding==4) = 1;
         else
             SubOffDec = ones(size(Parent1Dec));
         end

        %% Crossover for mask
        SubOffMask = Parent1Mask;
        for i = 1 : size(Parent1Mask, 1)
            if rand < 0.5
                index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
                index = index(TS(-Fitness(index)));
                SubOffMask(i,index) = 0;
            else
                index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
                index = index(TS(Fitness(index)));
                SubOffMask(i,index) = Parent2Mask(i,index);
            end
        end
        
        %% Mutation for mask
        if any(Problem.encoding~=4)
            chosenindex = groupIndex == chosengroups;
            for i = 1 : size(Parent1Mask, 1)
                if rand < 0.5
                    index = find(SubOffMask(i,:)&chosenindex(i,:));
                    index = index(TS(-Fitness(index)));
                    SubOffMask(i,index) = 0;
                else
                    index = find(~SubOffMask(i,:)&chosenindex(i,:));
                    index = index(TS(Fitness(index)));
                    SubOffMask(i,index) = 1;
                end
            end                    
        end  

        FinalMask = [FinalMask;SubOffMask];
        FinalDec = [FinalDec;SubOffDec];
        cutSzie = floor(size(Parent1, 2) / 2);
        % TaskNumber = [TaskNumber;Parent1(1:2:end).adds([], 1, 0); Parent2(2:2:end).adds([], 1, 0)];
        TaskNumber = [TaskNumber;Parent1(1:cutSzie).adds([], 1, 0); Parent2(cutSzie+1:end).adds([], 1, 0)];
    end

    %% Normal
    if ~isempty(NormalPop)
        SubOffMask = [];
        SubOffDec = [];

        Parent1 = NormalPop(1:2:end);
        Parent2 = NormalPop(2:2:end);

        %% Variation for mask
        OffMask = GAbinary(Parent1.adds([], 2, 0),Parent2.adds([], 2, 0),1,1);
        SubOffMask = [SubOffMask OffMask];

        %% Variation for dec
        OffDec = GAreal(Parent1.adds([], 3, 0),Parent2.adds([], 3, 0),Problem.lower,Problem.upper,1,20,1,20);
        SubOffDec = [SubOffDec OffDec];

        FinalMask = [FinalMask;SubOffMask];
        FinalDec = [FinalDec;SubOffDec];

        cutSzie = floor(size(Parent1, 2) / 2);
        % TaskNumber = [TaskNumber;Parent1(1:2:end).adds([], 1, 0); Parent2(2:2:end).adds([], 1, 0)];
        TaskNumber = [TaskNumber;Parent1(1:cutSzie).adds([], 1, 0); Parent2(cutSzie+1:end).adds([], 1, 0)];
    end

    SubOffspring = Problem.Evaluation(FinalDec.*FinalMask);
    SubOffspring.adds(TaskNumber, 1, 0);
    SubOffspring.adds(FinalMask,2, 0);
    SubOffspring.adds(FinalDec,3, 0);
end

function Offspring = GAbinaryCustom(Parent1,Parent2,MultiNonzeros)

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
        try
            G1 = MultiNonzeros(Parent1GuideVector(i, :)', :);
            G2 = MultiNonzeros(Parent2GuideVector(i, :)', :);
        catch err
            disp(err);
        end

        if rand < 0.5
            % Use maximum crossover
            index1 = xor(P1, P2);
            index2 = G1 | G2;

            index3 = find_two_nonzero_positions(index1);

            if length(index3) == 1
                if index2(index3) == 1
                    OffMask(i, index3) = 1;
                else
                    OffMask(i, index3) = 0;
                end
            elseif length(index3) == 2
                if all(index2(index3) == [1 1])
                    OffMask(i, index3) = 1;
                elseif all(index2(index3) == [0 0])
                    OffMask(i, index3) = 0;
                else
                    OffMask(i, index3(1)) = 1;
                    OffMask(i, index3(2)) = 0;
                end
            elseif isempty(index3)
                %% Uniform crossover
                [N,D] = size(P1);
                k     = rand(N,D) < 0.5;
                k(repmat(rand(N,1)>1,1,D)) = false;
                OffMask(i, k) = P2(k);
            end
        else
            % Use minimum crossover
            index1 = xor(P1, P2);
            index2 = G1 & G2;

            index3 = find_two_nonzero_positions(index1);

            if length(index3) == 1
                if index2(index3) == 1
                    OffMask(i, index3) = 1;
                else
                    OffMask(i, index3) = 0;
                end
            elseif length(index3) == 2
                if all(index2(index3) == [1 1])
                    OffMask(i, index3) = 1;
                elseif all(index2(index3) == [0 0])
                    OffMask(i, index3) = 0;
                else
                    OffMask(i, index3(1)) = 0;
                    OffMask(i, index3(2)) = 1;
                end
            elseif isempty(index3)
                %% Uniform crossover
                [N,D] = size(P1);
                k     = rand(N,D) < 0.5;
                k(repmat(rand(N,1)>1,1,D)) = false;
                OffMask(i, k) = P2(k);
            end
        end
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

function Offspring = GAbinary(Parent1,Parent2,proC,proM)
% Genetic operators for binary variables

    %% Uniform crossover
    [N,D] = size(Parent1);
    k     = rand(N,D) < 0.5;
    k(repmat(rand(N,1)>proC,1,D)) = false;
    Offspring    = Parent1;
    Offspring(k) = Parent2(k);
    
    %% Bit-flip mutation
    Site = rand(N,D) < proM/D;
    Offspring(Site) = ~Offspring(Site);
end

%% Utility function
function index = TS(Fitness)
% Binary tournament selection
    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
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
