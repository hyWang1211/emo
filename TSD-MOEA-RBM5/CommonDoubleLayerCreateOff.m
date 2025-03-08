function SubPopulation = CommonDoubleLayerCreateOff(Problem,Population,rmp,ns)

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
        % Parent1Dec     = Parent11.decs;
        % Parent2Dec     = Parent21.decs;

        Parent11Dec  = Parent11.adds([], 3);
        Parent21Dec  = Parent21.adds([], 3);

        Parent11Mask = Parent11.adds([], 2);
        Parent21Mask = Parent21.adds([], 2);

        % For dec
        OffDec1      = GAreal(Parent11Dec,Parent21Dec,Problem.lower,Problem.upper,1,20,1,20);
        % For mask
        OffMask1     = GAbinary(Parent11Mask,Parent21Mask,1,1);

        TaksNumberForOff1 = [Parent11.adds([], 1); Parent21.adds([], 1)];
    else
        OffDec1 = [];
        OffMask1 = [];
        TaksNumberForOff1 = [];
    end

    if ~isempty(Parent12)
        Parent12Dec     = Parent12.adds([], 3);
        Parent22Dec     = Parent22.adds([], 3);

        Parent12Mask     = Parent12.adds([], 2);
        Parent22Mask     = Parent22.adds([], 2);        

        % For dec
        OffDec2     = GAreal(Parent12Dec,Parent22Dec,Problem.lower,Problem.upper,0,20,1,20);
        % For mask
        OffMask2    = GAbinary(Parent12Mask,Parent22Mask,1,1);

        TaksNumberForOff2 = [Parent12.adds([], 1); Parent22.adds([], 1)];
    else
        OffDec2 = [];
        OffMask2 = [];
        TaksNumberForOff2 = [];
    end

    SubOffspring = Problem.Evaluation([OffDec1.*OffMask1;OffDec2.*OffMask2]);
    SubOffspring.adds([TaksNumberForOff1;TaksNumberForOff2], 1);
    SubOffspring.adds([OffMask1;OffMask2],2);
    SubOffspring.adds([OffDec1;OffDec2],3);
    for i = 1 : ns
        SubPopulation{i} = SubOffspring((SubOffspring.adds([], 1) == i));
    end
end

function Offspring = GAbinary(Parent1,Parent2,proC,proM)
% Genetic operators for binary variables

    %% Uniform crossover
    [N,D] = size(Parent1);
    k     = rand(N,D) < 0.5;
    k(repmat(rand(N,1)>proC,1,D)) = false;
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    Offspring1(k) = Parent2(k);
    Offspring2(k) = Parent1(k);
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