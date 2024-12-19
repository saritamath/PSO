% portfolioCost.m - Function to calculate portfolio risk and return
function z = portfolioCost(x, model)
    % Calculate portfolio risk using covariance matrix
    risk = sqrt(x' * model.sigma * x);
    
    % Calculate portfolio return (negative for minimization)
    ret = -x' * model.r;
    
    % Return both objectives
    z = [risk; ret];
end

% generateRandomSolution.m - Function to create valid random portfolio
function x = generateRandomSolution(model)
    nAsset = model.nAsset;
    K = model.K;
    
    % Select K random assets
    selectedAssets = randperm(nAsset, K);
    
    % Initialize weights
    x = zeros(nAsset, 1);
    
    % Generate random weights for selected assets
    weights = rand(K, 1);
    weights = weights / sum(weights);  % Normalize to sum to 1
    
    % Assign weights to selected assets
    x(selectedAssets) = weights;
    
    % Apply bounds and renormalize
    x = max(min(x, model.delta), model.epsilon);
    x = x / sum(x);
end

% applyConstraints.m - Function to enforce portfolio constraints
function x = applyConstraints(x, model)
    % Keep only top K assets (cardinality constraint)
    [~, indices] = sort(abs(x), 'descend');
    x(indices(model.K+1:end)) = 0;
    
    % Apply investment bounds
    x = max(min(x, model.delta), model.epsilon);
    
    % Ensure sum of weights equals 1
    x = x / sum(x);
end

% dominatesSolution.m - Function to check Pareto dominance
function dominates = dominatesSolution(a, b)
    % Check if solution a dominates solution b
    % Returns true if a is at least as good as b in all objectives
    % and better in at least one objective
    dominates = all(a <= b) && any(a < b);
end