classdef TSDMOEARBM5MSKEA < ALGORITHM
% <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% A Two-stage Optimization Framework for Sparse Large-scale Multiobjective Optimization
% ntp         --- 6    --- Number of template SubPopulation. Default = 6
% nhp         --- 4    --- Number of help SubPopulation. Default = 4
% rmp        --- 1    --- tranfer frequency

    methods
        function main(Algorithm,Problem)

            [ntp, nhp, rmp] = Algorithm.ParameterSet(6, 4, 1);
            % Select a template; here, we keep the same with SparseEA
            Template = zeros(1, Problem.D);
            % randomNumbers = randperm(Problem.D, ceil(Problem.D *1));
            % Template(randomNumbers) = 1;
            % Template Perturbation, and update ns
            [MultiNonzeros, ~, Fitness, prob_matrix] = TemplatePerturbation(Template, ntp, Problem);
            % ns = size(MultiNonzeros, 1);
            % Population initialization
            Dec = unifrnd(repmat(Problem.lower, Problem.N, 1),repmat(Problem.upper, Problem.N, 1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = zeros(Problem.N, Problem.D);
            GroupSize = floor(Problem.N / (ntp + nhp));
            GroupIdx = cell(1, ntp + nhp);

            
            %% TemplatePerturbationSubPopulationGeneration2
            length_cols = Problem.D;
            for i = 1: ntp
                rows = (i - 1)*GroupSize + 1: i* GroupSize;
                length_rows = length(rows);
                
                Mask(rows, :) = rand(length_rows, length_cols) < prob_matrix(i, :);
                GroupIdx{i} = (i - 1)*GroupSize + 1: i* GroupSize;
            end

            %% HelpSubPopulationInitialization
            [Mask, GroupIdx] = HelpSubPopulation(MultiNonzeros, Fitness, Mask, GroupIdx, GroupSize, ntp + nhp, Problem);
            
            Population =  Problem.Evaluation(Mask.*Dec);
            SubPopulation = cell(1, ntp + nhp);
            for i = 1: ntp + nhp
                % The property 'add' and function 'adds' in file SOLUTION.m need to be changed
                Population(GroupIdx{i}).adds(ones(size(GroupIdx{i},2), 1)* i, 1, 0);
                Population(GroupIdx{i}).adds(Mask(GroupIdx{i}, :), 2, 0);
                Population(GroupIdx{i}).adds(Dec(GroupIdx{i}, :), 3, 0);
                SubPopulation{i} = Population(GroupIdx{i});
            end

            % Get Front Number
            [~,~,Mask,FrontNo,~] = SPEA2_EnvironmentalSelection1(Population,Population.adds([], 3, 0),Population.adds([], 2, 0),Problem.N);
            Last_temp_num=0;
            sv=zeros(1,Problem.D);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                delta= Problem.FE/Problem.maxFE;
                %-------------update sv-----------%
                First_Mask=Mask(FrontNo==1,:);
                [temp_num,~]=size(First_Mask);
                temp_vote=sum(First_Mask,1);
                sv(1,:)=(Last_temp_num/(Last_temp_num+temp_num))*sv(1,:)+(temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                Last_temp_num=temp_num;
                %-------------update pv by sv-----------%
                if delta<0.618
                    Fitness=Fitness.*(1-sv)*sqrt(delta)+Fitness;
                end

                [SubPopulation, Rank] = Sort(ntp + nhp, SubPopulation);
                Population            = [SubPopulation{:}];
                % MatingPool            = TournamentSelection(2,length(Population)*2,[Rank{:}]);
                MatingPool            = TournamentSelection(2,Problem.N*2,[Rank{:}]);
                ParentPool            = Population(MatingPool);

                SubOffspring = ImprovedCreateOff2(Problem, ParentPool, rmp, MultiNonzeros, ntp, Fitness);
                % SubOffspring = ImprovedCreateOff1(Problem, ParentPool, rmp, MultiNonzeros);
                
                MixPopulation = [SubPopulation{:} SubOffspring];
                % [Population, ~, ~] = SPEA2_EnvironmentalSelection(MixPopulation,length([SubPopulation{:}]));
                [Population, ~, ~] = SPEA2_EnvironmentalSelection(MixPopulation,Problem.N);
                % [Population,~,~] = EnvironmentalSelection(MixPopulation,Problem.N);
     
                TempSubPopulation = cell(1, ntp + nhp);
                % emptySubPopulationIdx = [];

                % Divide
                for i = 1: ntp + nhp
                    TempSubPopulation{i} = Population(Population.adds([], 1, 0) == i);
                    TempSubPopulation{i}.adds(ones(length(TempSubPopulation{i}), 1) *i, 1, 1);
                end

                %% Fill the empty SubOffspring(use maping method) and evaluate them(Method 2-direct copy), (not complete).
                SubPopulation = TempSubPopulation;
                % Use rbm to update MultiNonzeros
                for i = 1 : ntp
                    if ~isempty(SubPopulation{i})
                        %% use RBM to learn a subspace
                        [other, rbm, ~, ~] = ModelTraining(SubPopulation{i}.adds([], 2, 0), SubPopulation{i}.adds([], 3, 0));

                        MultiNonzeros(i, other) = rbm.recover(rbm.reduce(MultiNonzeros(i, other)));
                        % MultiNonzeros(i, allZero) = 0;
                        % MultiNonzeros(i, allOne) = 1;
                    else
                        % disp('err');
                    end
                end

                [~,~,Mask,FrontNo,~] = SPEA2_EnvironmentalSelection1(Population,Population.adds([], 3, 0),Population.adds([], 2, 0),Problem.N);
            end
        end
    end
end