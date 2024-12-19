function particles = InitializeParticles(nPop, nAsset)
    for i = 1:nPop
        particles(i).Position = rand(nAsset, 1); % Positions initiales
        particles(i).Velocity = zeros(nAsset, 1); % Vitesses initiales
        particles(i).Cost = Inf;
        particles(i).Best.Position = [];
        particles(i).Best.Cost = Inf;
        particles(i).Z = randi([0, 1], nAsset, 1); % Variables binaires
    end
end
