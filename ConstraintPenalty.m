function penalty = ConstraintPenalty(x, Z, S, K, epsilon, delta)
    % Pénalité pour la contrainte de budget
    budget_penalty = abs(sum(x) - 1);

    % Pénalité pour la contrainte de cardinalité
    cardinality_penalty = abs(sum(Z) - K);

    % Pénalité pour les bornes de quantité
    quantity_penalty = sum(max(0, epsilon .* Z - x)) + sum(max(0, x - delta .* Z));

    % Pénalité pour la préaffectation
    preassign_penalty = sum(abs(S - Z));

    % Total des pénalités
    penalty = budget_penalty + cardinality_penalty + quantity_penalty + preassign_penalty;
end
