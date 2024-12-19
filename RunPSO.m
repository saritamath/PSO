[out, risks, returns] = RunPSO(model, nPop, MaxIt);

disp('Meilleure solution :');
disp(out.BestSol.Position);

disp('Risque associ� :');
disp(out.BestSol.Risk);

disp('Rendement associ� :');
disp(out.BestSol.Return);

% Tracer la fronti�re de Pareto
figure;
plot(risks, returns, 'bo');
hold on;
plot(out.ParetoFront(:,1), out.ParetoFront(:,2), 'r-', 'LineWidth', 2);
xlabel('Risque');
ylabel('Rendement');
title('Fronti�re de Pareto');
grid on;
