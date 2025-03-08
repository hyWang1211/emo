function SubPopulation = CreateOff(Problem,Population, MultiNonzeros, rmp, ns)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parent selection

    % store parent skill factor
    Parent11 = [];
    Parent21 = [];
    Parent12 = [];
    Parent22 = [];

    % for case1 to store mask and dec
    MaskParent11 = [];
    MaskParent21 = [];
    DecParent11  = [];
    DecParent21  = [];

    % for case2 to store mask and dec
    MaskParent12 = [];
    MaskParent22 = [];
    DecParent12  = [];
    DecParent22  = [];
    
    for i = 1 : floor(length(Population)/2)
        P1 = Population(i);
        P2 = Population(i+floor(length(Population)/2));
        if (P1.add{1} == P2.add{1}) || (rand<rmp)
            % for mask
            MaskParent11 = [MaskParent11; P1.add{2}];
            MaskParent21 = [MaskParent21; P2.add{2}];
            % for dec
            DecParent11  = [DecParent11; P1.add{3}];
            DecParent21  = [DecParent21; P2.add{3}];
            % for parent
            Parent11     = [Parent11; P1.add{1}];
            Parent21     = [Parent21; P2.add{1}];
        else
            %% only mutation
            % for mask
            MaskParent12 = [MaskParent12; P1.add{2}];
            MaskParent22 = [MaskParent22; P2.add{2}];
            % for dec
            DecParent12  = [DecParent12; P1.add{3}];
            DecParent22  = [DecParent22; P2.add{3}];
            % for parent
            Parent12     = [Parent12; P1.add{1}];
            Parent22     = [Parent22; P2.add{1}];
        end
    end

    %% Offspring generation
    if ~isempty(Parent11)
        % for mask
        OffMask1    = GAbinary(MaskParent11, MaskParent21, Parent11, Parent21, MultiNonzeros, 1, 1);
        % for dec
        OffDec1     = GAreal(DecParent11, DecParent21, Problem.lower, Problem.upper, 1, 20, 1, 20);
        OffSpring1  = OffMask1.*OffDec1;
        TaksNumberForOff1 = [Parent11; Parent21];
    else
        OffMask1    = [];
        OffDec1     = [];
        OffSpring1  = [];
        TaksNumberForOff1 = [];
    end
    if ~isempty(Parent12)
        %% only mutation
        % for mask
        % OffMask2    = GAbinary(MaskParent12, MaskParent22, 0, 0.5);
        OffMask2    = GAbinary(MaskParent12, MaskParent22, Parent12, Parent22, MultiNonzeros, 1, 1);
        % for dec
        OffDec2     = GAreal(DecParent12, DecParent22, Problem.lower, Problem.upper, 0, 20, 1, 20);
        OffSpring2  = OffMask2.*OffDec2;
        TaksNumberForOff2 = [Parent12; Parent22];
    else
        OffMask2    = [];
        OffDec2     = [];
        OffSpring2  = [];
        TaksNumberForOff2 = [];
    end
    SubOffspring = Problem.Evaluation([OffSpring1;OffSpring2]);
    SubOffspring.adds([TaksNumberForOff1;TaksNumberForOff2], 1);
    SubOffspring.adds([OffMask1;OffMask2],2);
    SubOffspring.adds([OffDec1;OffDec2],3);
    for i = 1 : ns
        SubPopulation{i} = SubOffspring((SubOffspring.adds([], 1) == i));
    end
end

function Offspring = GAbinary(Parent1,Parent2, Parent11, Parent21, MultiNonzeros, ~, proM)
% Genetic operators for binary variables

    %% replace crossover
    [N,D] = size(Parent1);
    P1_nonzeroPos = [];
    P2_nonzeroPos = [];
    P1_zeroPos    = [];
    P2_zeroPos    = [];
    Parent11Base = MultiNonzeros(Parent11, :);
    Parent21Base = MultiNonzeros(Parent21, :);
    for i = 1 : size(Parent1, 1)
        % Nonzero
        P1_npos = find(Parent11Base(i, :) == 1);
        r1 = randi(size(P1_npos));
        P1_nonzeroPos = [P1_nonzeroPos; P1_npos(r1)];

        P2_npos = find(Parent21Base(i, :) == 1);
        r2 = randi(size(P2_npos));
        P2_nonzeroPos = [P2_nonzeroPos; P2_npos(r2)];

        % Zero
        P1_zpos = find(Parent11Base(i, :) == 0);
        r1 = randi(size(P1_zpos));
        P1_zeroPos = [P1_zeroPos; P1_zpos(r1)];

        P2_zpos = find(Parent21Base(i, :) == 0);
        r2 = randi(size(P2_zpos));
        P2_zeroPos = [P2_zeroPos; P2_zpos(r2)];
    end
    
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    
    Offspring1(P1_nonzeroPos) = Parent2(P2_nonzeroPos);
    Offspring1(P1_zeroPos)    = Parent2(P2_zeroPos);

    Offspring2(P2_nonzeroPos) = Parent1(P1_nonzeroPos);
    Offspring2(P2_zeroPos)    = Parent1(P1_zeroPos);

    Offspring     = [Offspring1;Offspring2];
    
    %% Bit-flip mutation
    Site = rand(2*N,D) < proM/D;
    Offspring(Site) = ~Offspring(Site);
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
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
             
    %% Polynomial mutation
    Lower = repmat(lower,2*N,1);
    Upper = repmat(upper,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end