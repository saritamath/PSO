% MainPSO.m - Complete Portfolio Optimization using PSO with local functions
function MainPso
clc;
clear;
close all;

%% Problem Parameters
%QP=Hang_Seng_37;
 %QP=Nikkei_225;
 %QP=DAX_100_85;
 %QP=S_P_100_98;
 %QP=FTSE_100_89;
 %QP=SP500_478;
 QP=EuroStoxx50;
r=QP.Rr;
sigma=QP.Dd;
nAsset = length(r); % Number of assets
K = 10; % Cardinality (maximum number of selected assets)

% Create random covariance matrix for example
sigma = rand(nAsset); 
sigma = (sigma + sigma') / 2 + nAsset * eye(nAsset); % Make symmetric and positive definite

% Initialize asset returns and constraints
r = rand(nAsset, 1); % Mean return of assets
epsilon = 0.01 * ones(nAsset, 1); % Lower bound on weights
delta = 0.9 * ones(nAsset, 1); % Upper bound on weights
%S = zeros(nAsset, 1); % Pre-allocation
S=2;
%% Model Structure
model.sigma = sigma;
model.r = r;
model.K = K;
model.epsilon = epsilon;
model.delta = delta;
model.S = S;
model.nAsset = nAsset;

%% PSO Parameters
nPop = 100;     % Population size (number of particles)
MaxIt = 1000;   % Maximum number of iterations
w = 0.7298;    % Inertia weight
c1 = 1.4962;   % Personal learning coefficient
c2 = 1.4962;   % Global learning coefficient

%% Initialize PSO Population
% Create empty particle structure
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create population array
particle = repmat(empty_particle, nPop, 1);

% Initialize global best
GlobalBest.Cost = inf(1, 2);  % [Risk, -Return]
GlobalBest.Position = [];

% Initialize particles
for i = 1:nPop
    % Generate random valid position
    particle(i).Position = generateRandomSolution(model);
    
    % Initialize velocity
    particle(i).Velocity = zeros(size(particle(i).Position));
    
    % Evaluate cost
    particle(i).Cost = portfolioCost(particle(i).Position, model);
    
    % Update personal best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update global best
    if dominatesSolution(particle(i).Cost, GlobalBest.Cost)
        GlobalBest.Position = particle(i).Position;
        GlobalBest.Cost = particle(i).Cost;
    end
end

%% PSO Main Loop
bestCosts = zeros(MaxIt, 2);

for it = 1:MaxIt
    for i = 1:nPop
        % Update velocity
        particle(i).Velocity = w * particle(i).Velocity ...
            + c1 * rand(size(particle(i).Position)) .* (particle(i).Best.Position - particle(i).Position) ...
            + c2 * rand(size(particle(i).Position)) .* (GlobalBest.Position - particle(i).Position);
        
        % Update position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply constraints
        particle(i).Position = applyConstraints(particle(i).Position, model);
        
        % Evaluate cost
        particle(i).Cost = portfolioCost(particle(i).Position, model);
        
        % Update personal best
        if dominatesSolution(particle(i).Cost, particle(i).Best.Cost)
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update global best
            if dominatesSolution(particle(i).Cost, GlobalBest.Cost)
                GlobalBest.Position = particle(i).Position;
                GlobalBest.Cost = particle(i).Cost;
            end
        end
    end
    
    % Store best cost
    bestCosts(it,:) = GlobalBest.Cost;
    
    % Display iteration info
    fprintf('Iteration %d: Risk = %.4f, Return = %.4f\n', ...
        it, GlobalBest.Cost(1), -GlobalBest.Cost(2));
end

%% Results Visualization
figure;
plot(bestCosts(:,1), -bestCosts(:,2), 'o-');
xlabel('Portfolio Risk');
ylabel('Portfolio Return');
title('Pareto Front Approximation');
grid on;

% Display final portfolio weights
fprintf('\nOptimal Portfolio Weights:\n');
for i = 1:nAsset
    if GlobalBest.Position(i) > 1e-6
        fprintf('Asset %d: %.4f\n', i, GlobalBest.Position(i));
    end
end
headers = {'Risk', 'Return'};
% Ajouter les données (avec en-têtes si nécessaire)
headers = {'Risk', 'Return'};
data = [bestCosts(:,1), -bestCosts(:,2)];

% Ajouter les en-têtes
headers = {'Risk', 'Return'};
%xlswrite('EUROSTOX.xlsx', headers, 'sheet2', 'A1'); % Écrire les en-têtes

% Ajouter les données (Risk et Return)
% Ajouter les en-têtes
headers = {'Risk', 'Return'};
xlswrite('SP500PjjSO.xlsx', headers, 'sheet2', 'A1'); % Écrire les en-têtes
% Ajouter les données (Risk et Return)
data = [bestCosts(:,1), -bestCosts(:,2)];
data = [bestCosts(:,1), -bestCosts(:,2)];

%xlswrite('EuJJrostoxPSO.xlsx', data, 'sheet2', 'A2'); % Écrire les données


end
%% Local Functions
% Note: All function definitions must come after the main script code

function z = portfolioCost(x, model)
    % Calculate portfolio risk using covariance matrix
    risk = sqrt(x' * model.sigma * x);
    
    % Calculate portfolio return (negative for minimization)
    ret = -x' * model.r;
    
    % Return both objectives
    z = [risk; ret];
end

function x = generateRandomSolution(model)
    % Select K random assets
    selectedAssets = randperm(model.nAsset, model.K);
    
    % Initialize weights
    x = zeros(model.nAsset, 1);
    
    % Generate random weights for selected assets
    weights = rand(model.K, 1);
    weights = weights / sum(weights);  % Normalize to sum to 1
    
    % Assign weights to selected assets
    x(selectedAssets) = weights;
    
    % Apply bounds and renormalize
    x = max(min(x, model.delta), model.epsilon);
    x = x / sum(x);
end

function x = applyConstraints(x, model)
    % Keep only top K assets (cardinality constraint)
    [~, indices] = sort(abs(x), 'descend');
    x(indices(model.K+1:end)) = 0;
    
    % Apply investment bounds
    x = max(min(x, model.delta), model.epsilon);
    
    % Ensure sum of weights equals 1
    x = x / sum(x);
end

function dominates = dominatesSolution(a, b)
    % Determines if solution a dominates solution b in the context of
    % multi-objective minimization (both risk and negative return are minimized)
    % 
    % Parameters:
    %   a: First solution vector [risk; negative_return]
    %   b: Second solution vector [risk; negative_return]
    %
    % Returns:
    %   dominates: logical scalar indicating if a dominates b
    
    % Get the number of objectives
    n_obj = numel(a);
    
    % Initialize counters for our comparison
    is_better_in_any = false;
    is_worse_in_any = false;
    
    % Compare each objective
    for i = 1:n_obj
        if a(i) < b(i)
            is_better_in_any = true;
        elseif a(i) > b(i)
            is_worse_in_any = true;
            break;  % Can exit early if we find any worse objective
        end
    end
    
    % Solution a dominates b if it's better in at least one objective
    % and not worse in any objective
    dominates = is_better_in_any && ~is_worse_in_any;
end