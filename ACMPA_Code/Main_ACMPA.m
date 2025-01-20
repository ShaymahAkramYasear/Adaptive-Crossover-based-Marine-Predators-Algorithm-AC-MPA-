%% #############################################################################################
% Adaptive-Crossover-based-Marine-Predators-Algorithm-AC-MPA-
% This is the original code of Adaptive Crossover-based Marine Predators Algorithm (AC-MPA). 
% The code was developed by Assist. Prof. Dr. Shaymah Akram Yasear.
% 
% The mathematical details of AC-MPA can be found in: Shaymah Akram Yasear, 
% Adaptive crossover-based marine predators algorithm for global optimization problems, Journal of 
% Computational Design and Engineering, Volume 11, Issue 4, August 2024, Pages 124–150, 
% https://doi.org/10.1093/jcde/qwae060
% 
% Contact details of author: shayma.akram.yasear@gmail.com
% https://www.researchgate.net/profile/Shaymah-Yasear
% https://www.scopus.com/authid/detail.uri?authorId=57194554063
% https://www.webofscience.com/wos/author/record/AAN-5105-2021
%% #############################################################################################

function []= Main_ACMPA()
    % Configuration parameters
    max_runs = 30;
    popSize = 30;
    problem_sizes = [30 50 100];  % Array of different dimensions to test
    num_prbs = 30;               % Number of test problems
    max_FES = 500;               % Maximum function evaluations

    % Main loop for different problem dimensions
    for dimNo = 1:length(problem_sizes)
        current_dimension = problem_sizes(dimNo);
        
        % Set up options structure
        options = struct();
        options.MaxFES = max_FES;
        options.popSize = popSize;
        options.nVar = current_dimension;
        options.lb = -100 * ones(1, current_dimension);
        options.ub = 100 * ones(1, current_dimension);
        options.objf = str2func('cec17_func');

        % Initialize results arrays for current dimension
        objective_values = zeros(num_prbs, max_runs);
        decision_variables = cell(num_prbs, max_runs);
        convergence_curves = cell(num_prbs, max_runs);

        % Loop through each problem from CEC2017
        for fn = 1:num_prbs
            % Skip specific problems for 2D case
            if (current_dimension == 2) && ismember(fn, [4 6 8 16])
                fprintf('Skipping problem %d for dimension %d\n', fn, current_dimension);
                continue;
            end

            options.CEC_fun_no = fn;

            % Run algorithm multiple times
            for runNo = 1:max_runs
                try
                    % Run ACMPA algorithm
                    [ObjectiveValue, DecisionVariable, Convergence_curve] = ACMPA(options);

                    % Store results
                    objective_values(fn, runNo) = ObjectiveValue;
                    decision_variables{fn, runNo} = DecisionVariable;
                    convergence_curves{fn, runNo} = Convergence_curve;

                    % Display progress
                    fprintf('Dimension: %d, Problem: %d, Run: %d, Objective: %.4e\n', ...
                        current_dimension, fn, runNo, ObjectiveValue);

                catch ME
                    % Error handling
                    fprintf('Error occurred - Dimension: %d, Problem: %d, Run: %d\n', ...
                        current_dimension, fn, runNo);
                    fprintf('Error message: %s\n', ME.message);
                end
            end

            % Save results after each problem
            results = struct('objective_values', objective_values, ...
                           'decision_variables', decision_variables, ...
                           'convergence_curves', convergence_curves, ...
                           'dimension', current_dimension);
            
        end

        % Save complete results for current dimension
        save(sprintf('ACMPA_Results_Dim%d_Complete.mat', current_dimension), 'results');
    end
end
