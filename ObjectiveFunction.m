function [risk, ret] = ObjectiveFunction(x, sigma, r)
    % Calcul du risque
    risk = x' * sigma * x;

    % Calcul du rendement
    ret = r' * x;
end

