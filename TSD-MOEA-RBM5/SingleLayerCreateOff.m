function SubPopulation = SingleLayerCreateOff(Problem,Population,rmp,ns)

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
        Parent1Dec     = Parent11.decs;
        Parent2Dec     = Parent21.decs;

        OffDec1        = OperatorGA(Problem,[Parent1Dec;Parent2Dec]);
        % OffDec1(:,end) = [Parent1Dec(:,end);Parent2Dec(:,end)];
        TaksNumberForOff1 = [Parent11.adds([], 1); Parent21.adds([], 1)];
    else
        OffDec1 = [];
        TaksNumberForOff1 = [];
    end

    if ~isempty(Parent12)
        Parent1Dec     = Parent12.decs;
        Parent2Dec     = Parent22.decs;

        OffDec2        = OperatorGA(Problem,[Parent1Dec;Parent2Dec],{0,20,1,20});
        % OffDec2(:,end) = [Parent1Dec(:,end);Parent2Dec(:,end)];
        TaksNumberForOff2 = [Parent12.adds([], 1); Parent22.adds([], 1)];
    else
        OffDec2 = [];
        TaksNumberForOff2 = [];
    end

    SubOffspring = Problem.Evaluation([OffDec1;OffDec2]);
    SubOffspring.adds([TaksNumberForOff1;TaksNumberForOff2], 1);
    % SubOffspring.adds([OffMask1;OffMask2],2);
    % SubOffspring.adds([OffDec1;OffDec2],3);
    for i = 1 : ns
        SubPopulation{i} = SubOffspring((SubOffspring.adds([], 1) == i));
    end


    % SubOffspring = Divide(Problem.Evaluation([OffDec1;OffDec2]),length(Problem.SubD));
end