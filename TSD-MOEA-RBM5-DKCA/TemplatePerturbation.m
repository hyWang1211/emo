function [MultiNonzeros, true_ns, Fitness, prob_matrix] = TemplatePerturbation(Template, ns, Problem)
    MultiNonzeros   = zeros(ns, Problem.D);
    Fitness         = zeros(1,Problem.D);
    % Create Perturbation matrix
    PerturbationMatrix = zeros(Problem.D, Problem.D);
    for i = 1: Problem.D
        PerturbationMatrix(i, i) = ~Template(i);
    end

    PerturbationMatrix = [Template; PerturbationMatrix];

    for i = 1: ns
        % Sampling the real Dec matrix
        Dec = unifrnd(repmat(Problem.lower, Problem.D + 1, 1),repmat(Problem.upper, Problem.D + 1, 1));
        % Dec = repmat(Problem.lower, Problem.D + 1, 1) + rand(Problem.D + 1, length(Problem.lower)) .* (repmat(Problem.upper, Problem.D + 1, 1) - repmat(Problem.lower, Problem.D + 1, 1));
        % Combine Dec with PerturbationMatrix to get Population
        Population = Problem.Evaluation(PerturbationMatrix.*Dec);

        %% Nonzero variable detection
        % Sort
        [FrontNo, ~] = NDSort(Population.objs, Population.cons, Problem.D + 1);
        Fitness = Fitness + FrontNo(2:end);
        % Detection
        MultiNonzeros(i, :) = (FrontNo(2:end) < FrontNo(1) & diag(PerturbationMatrix(2:end, :))' == 1) ...
            | (FrontNo(2:end) > FrontNo(1) & diag(PerturbationMatrix(2:end, :))' == 0);
    end
    % [MultiNonzeros, ~] = unique(MultiNonzeros, 'rows');
    true_ns = size(MultiNonzeros, 1);

    % golbal_prob probability for nonzero variables
    golbal_prob = 0.5 *mean(MultiNonzeros, 1);

    % local_prob probability for nonzero variables
    local_prob = 0.5 *(MultiNonzeros *0.8 + ~MultiNonzeros *0.2);

    % final probability for nonzero vairables
    prob_matrix = local_prob + golbal_prob;
end

% function [MultiNonzeros, true_ns, Fitness, prob_matrix] = TemplatePerturbation(Template, ns, Problem)
%     MultiNonzeros   = zeros(ns, Problem.D);
%     Fitness         = zeros(1,Problem.D);
%     % Create Perturbation matrix
%     PerturbationMatrix = zeros(Problem.D, Problem.D);
%     for i = 1: Problem.D
%         PerturbationMatrix(i, i) = ~Template(i);
%     end
% 
%     % PerturbationMatrix = [Template; PerturbationMatrix];
% 
%     for i = 1: ns
%         % Sampling the real Dec matrix
%         % Dec = unifrnd(repmat(Problem.lower, Problem.D + 1, 1),repmat(Problem.upper, Problem.D + 1, 1));
%         Dec = repmat(Problem.lower, Problem.D, 1) + rand(Problem.D, length(Problem.lower)) .* (repmat(Problem.upper, Problem.D, 1) - repmat(Problem.lower, Problem.D, 1));
%         % Combine Dec with PerturbationMatrix to get Population
%         Population = Problem.Evaluation(PerturbationMatrix.*Dec);
% 
%         %% Nonzero variable detection
%         % Sort
%         [FrontNo, ~] = NDSort(Population.objs, Population.cons, Problem.D);
%         [~, I] = sort(FrontNo);
%         FrontNoIdx = I(ceil(length(I) / 10));
%         % FrontNoIdx = I(20);
% 
%         Fitness = Fitness + FrontNo;
%         % Detection
%         MultiNonzeros(i, :) = (FrontNo(1:end) < FrontNo(FrontNoIdx) & diag(PerturbationMatrix(1:end, :))' == 1) ...
%             | (FrontNo(1:end) > FrontNo(FrontNoIdx) & diag(PerturbationMatrix(1:end, :))' == 0);
%     end
%     % [MultiNonzeros, ~] = unique(MultiNonzeros, 'rows');
%     true_ns = size(MultiNonzeros, 1);
% 
%     % golbal_prob probability for nonzero variables
%     golbal_prob = 0.5 *mean(MultiNonzeros, 1);
% 
%     % local_prob probability for nonzero variables
%     local_prob = 0.5 *(MultiNonzeros *0.8 + ~MultiNonzeros *0.2);
%     prob_matrix = local_prob + golbal_prob;
%     % MultiNonzeros = MultiNonzeros(1: ceil(4*true_ns/5), :);
%     % true_ns = size(MultiNonzeros, 1);
% end
