classdef TSDMOEARBM5 < ALGORITHM
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

            %% Optimization
            % while Algorithm.NotTerminated([SubPopulation{:}])
            count = 0;
            GV = [];
            Alive = [];
            Problem.FE = 0;
            while Algorithm.NotTerminated(Population)
                count = count + 1;
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

                    % if isempty(TempSubPopulation{i})
                    %     emptySubPopulationIdx = [emptySubPopulationIdx i];
                    % end
                end
                % notEmptySubPopulationIdx = setdiff(1:ns, emptySubPopulationIdx);

                %% Fill the empty SubOffspring(use maping method) and evaluate them(Method 1-linear mapping)
                % for i = 1: length(emptySubPopulationIdx)
                %     subIdx = emptySubPopulationIdx(i);
                %     originalSubPopulation = SubPopulation{subIdx};
                % 
                %     % random select a subP from notEmptySubPopulationIdx to
                %     % help originalSubPopulation
                %     helpSubPopulation = TempSubPopulation{notEmptySubPopulationIdx(randi(length(notEmptySubPopulationIdx)))};
                % 
                %     minIndividual = min(min(length(originalSubPopulation), length(helpSubPopulation)), 5);
                %     % disp(minIndividual);
                % 
                %     % sampled solutions for transfer
                %     selectedIdx = randperm(length(helpSubPopulation), minIndividual);
                %     newHelpSubPopulation = helpSubPopulation(selectedIdx);
                % 
                %     selectedIdx2 = randperm(length(originalSubPopulation), minIndividual);
                %     newOriginalSubPopulation = originalSubPopulation(selectedIdx2);
                % 
                %     targetDomain = newOriginalSubPopulation.decs;
                %     sourceDomain = newHelpSubPopulation.decs;
                % 
                %     % matrix learning(A*M = B)
                %     M = pinv(sourceDomain' *sourceDomain) *sourceDomain' *targetDomain;
                % 
                %     % get transfer solutions
                %     [FrontNo, ~] = NDSort(helpSubPopulation.objs, length(helpSubPopulation));
                %     transferMatrix = helpSubPopulation(FrontNo == 1).decs;
                % 
                %     % transfer
                %     targetDecs = transferMatrix *M;
                % 
                %     % evaluation
                %     TempSubPopulation{subIdx} = Problem.Evaluation(targetDecs);
                %     TempSubPopulation{subIdx}.adds(ones(length(TempSubPopulation{subIdx}), 1) *i, 1, 1);
                %     TempSubPopulation{subIdx}.adds(targetDecs, 2, 1);
                %     TempSubPopulation{subIdx}.adds(double(logical(targetDecs)), 3, 1);
                % end

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

                
                % if (Problem.FE >= Problem.maxFE)
                %     save('TSD_count', 'count');
                % end
                % if count == 400 || count == 1400 || count == 2400
                %     GV = [GV;MultiNonzeros];
                %     % save('HSDEA_GV', 'MultiNonzeros');
                %     alive = zeros(1, 6);
                %     for i = 1: ntp
                %         if ~isempty(SubPopulation{i})
                %             alive(i) = 1;
                %         end
                %     end
                %     Alive = [Alive; alive];
                % 
                %     if count == 2400
                %         save('HSDEA_GV', 'GV');
                %         save('HSDEA_Alive', 'Alive');
                %         break;
                %     end
                % end
            end
        end
    end
end