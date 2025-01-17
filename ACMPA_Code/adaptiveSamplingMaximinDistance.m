function X = adaptiveSamplingMaximinDistance(lb, ub, dim, N, Iteration_Limit)
    % Inputs:
    % - lb, ub: Search space boundaries
    % - N: Number of points to sample
    % - Iteration_Limit: Maximum number of iterations for finding the best candidate
    
    % Output:
    % - X: Set of sampled points
    
    % Initialization
    X = []; % Initialize an empty set X
   
    % Generate a random initial point x0 within the bounds of the search space S
    x0 = lb + (ub - lb) .* rand(1,dim);
    X = [X; x0]; % Add x0 to X
    
    % Iterative Sampling
    for i = 2:N
        best_min_distance = 0;
        best_candidate = [];
        for j = 1:Iteration_Limit
            % Generate a random candidate point x'
            x_candidate = lb + (ub - lb) .* rand(1,dim);
            
            % Calculate D(x') according to Equation (9)
            distance = calculateDistance(x_candidate, X);
            
            if distance > best_min_distance
                best_min_distance = distance;
                best_candidate = x_candidate;
            end
        end
        X = [X; best_candidate]; % Add best_candidate to X
    end
end

function distance = calculateDistance(x, X)
    % Calculate the minimum distance of point x to the existing set of points X
    distances = sqrt(sum((X - x).^2, 1));
    distance = min(distances);
end
