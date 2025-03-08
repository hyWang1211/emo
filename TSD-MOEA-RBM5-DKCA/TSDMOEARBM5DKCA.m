classdef TSDMOEARBM5DKCA < ALGORITHM
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
            k = 4;
            N2 = 4;

            
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
            t           = 0.3;
            T           = t*Problem.D;
            Mask        = eye(Problem.D);
            Dec         = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
            Dec2        = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
            dim         = DimSelection(Problem,Mask,Dec,Dec2,Problem.D,T);
            dim(dim==0) = [];
            d1          = size(dim,2);  
            init_time   = 4;
            TempPop_d1  = [];
            TMask_d1   = [];
            for i = 1 : init_time
                Mask_dim = eye(d1);   %d1*d1
                Mask_d1  = zeros(d1,Problem.D);   %d1*D
                for j = 1 : d1
                     Mask_d1(j,dim(find(Mask_dim(j,:)))) = 1;
                end
                % if REAL
                    Dec_d1 = unifrnd(repmat(Problem.lower,d1,1),repmat(Problem.upper,d1,1));
                % else
                    % Dec_d1 = ones(d1,Problem.D);
                % end
                Population_d1 = Problem.Evaluation(Dec_d1.*Mask_d1);  
                TMask_d1      = [TMask_d1;Mask_dim];
                TempPop_d1    = [TempPop_d1,Population_d1]; 
            end
            if init_time*d1 < Problem.N
                num_e = init_time*d1;
            else
                num_e = Problem.N;
            end
            [Population_d1,Mask_dim,FrontNo_d1,CrowdDis_d1] = EnvironmentalSelection_withoutDec([Population_d1,TempPop_d1],[Mask_dim;TMask_d1],num_e);


            generation = 1;
            dim         = DimSelection(Problem,Mask,Dec,Dec2,Problem.D,T);
            dim(dim==0) = [];
            Dec = Population.adds([], 2, 0);
            while Algorithm.NotTerminated(Population)
                [SubPopulation, Rank] = Sort(ntp + nhp, SubPopulation);
                Population            = [SubPopulation{:}];
                MatingPool            = TournamentSelection(2,Problem.N*2,[Rank{:}]);
                ParentPool            = Population(MatingPool);

                SubOffspring = ImprovedCreateOff2(Problem, ParentPool, rmp, MultiNonzeros, ntp, Fitness);
                
                MixPopulation = [SubPopulation{:} SubOffspring];
                [Population, ~, ~] = SPEA2_EnvironmentalSelection(MixPopulation,Problem.N);
                
                TempSubPopulation = cell(1, ntp + nhp);

                MatingPool_d1 = TournamentSelection(2,2*Problem.N,FrontNo_d1,-CrowdDis_d1);
                [~,Mask_dim_off] = Operator1(Problem,Dec(MatingPool_d1,:),Mask_dim(MatingPool_d1,:),Fitness, 1);
                [Population_d1,Mask_dim,FrontNo_d1,CrowdDis_d1] = EnvironmentalSelection_withoutDec([Population_d1,SubOffspring],[Mask_dim;Mask_dim_off],Problem.N);
                Dec = Population.adds([], 2, 0);

                for i = 1 : size(Mask_dim,1)
                    non_zero_temp = find(Mask_dim(i,:));  
                    non_zero(1,i) = size(non_zero_temp,2);  
                end
                non_zero_num = mode(non_zero(1,:),2);
                if generation > 1   
                    if num_temp == non_zero_num   
                        same_sparsity = same_sparsity+1;
                    else
                        same_sparsity = 1;
                    end
                else
                    same_sparsity = 1;
                end
                num_temp   = non_zero_num;
                generation = generation+1;
                if same_sparsity > k
                    dim_base = DimBase(Mask_dim,dim,num_temp);
                    for loc = dim_base
                        Fitness(loc) = Fitness(loc) - ceil(1/N2*Fitness(loc));
                    end
                end




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
                    else
                    end
                end


            end
        end
    end
end