function penalty = ConstraintPenalty(x, Z, S, K, epsilon, delta)
    % P�nalit� pour la contrainte de budget
    budget_penalty = abs(sum(x) - 1);

    % P�nalit� pour la contrainte de cardinalit�
    cardinality_penalty = abs(sum(Z) - K);

    % P�nalit� pour les bornes de quantit�
    quantity_penalty = sum(max(0, epsilon .* Z - x)) + sum(max(0, x - delta .* Z));

    % P�nalit� pour la pr�affectation
    preassign_penalty = sum(abs(S - Z));

    % Total des p�nalit�s
    penalty = budget_penalty + cardinality_penalty + quantity_penalty + preassign_penalty;
end
