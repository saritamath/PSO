[out, risks, returns] = RunPSO(model, nPop, MaxIt);

disp('Meilleure solution :');
disp(out.BestSol.Position);

disp('Risque associé :');
disp(out.BestSol.Risk);

disp('Rendement associé :');
disp(out.BestSol.Return);

% Tracer la frontière de Pareto
figure;
plot(risks, returns, 'bo');
hold on;
plot(out.ParetoFront(:,1), out.ParetoFront(:,2), 'r-', 'LineWidth', 2);
xlabel('Risque');
ylabel('Rendement');
title('Frontière de Pareto');
grid on;
