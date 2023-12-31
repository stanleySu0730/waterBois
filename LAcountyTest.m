% Given LA County COVID-19 data
global ConfCase;
ConfCase = [13823 15165 16449 17567 18545 19159 19567 20460 21017 22522 23233 24262 24936 25699 26238 27866 28665 29526 30334 31241 31703 32269 33247 34552 35447 36324 37374];
global TotalPopulation;
TotalPopulation = 9651332;

% Parameter estimation using fmincon
LowerBounds = [1.00e-9, 1.00e-9, 1.00e-9];
UpperBounds = [1.00e-7, 1.00e-3, 1.00e-3];
xstart = 0.5 * (LowerBounds + UpperBounds);

problem = createOptimProblem('fmincon', 'objective', @SIV_RUN_ODE45, 'x0', xstart, 'lb', LowerBounds, 'ub', UpperBounds);
problem.options = optimoptions('fmincon', 'MaxFunEvals', 9999, 'MaxIter', 9999);

numstartpoints = 10;
ms = MultiStart('UseParallel', true, 'Display', 'iter');
[b, fval, exitflag, output, manymins] = run(ms, problem, numstartpoints);

% Extract solutions
SIVParameters = zeros(length(manymins), 3);
for i = 1:length(manymins)
    SIVParameters(i, :) = manymins(i).X;
end

% Plotting the "best" solution
best_params = manymins(1).X;
[t, y] = ode45(@SIV, 1:length(ConfCase), [TotalPopulation - ten_dayCumulative - 0.95 * 20 * CumulativeTo10DayStart; ten_dayCumulative; 0.95 * 20 * CumulativeTo10DayStart; ConfCase(1)], best_params(1), best_params(2), best_params(3));
figure;
subplot(2, 2, 1);
plot(1:length(ConfCase), y(:, 4), 'b', 'LineWidth', 2);
hold on;
scatter(1:length(ConfCase), ConfCase, 'r', 'filled');
title('Cumulative Cases');
xlabel('Days');
ylabel('Cases');
legend('Model', 'Data');

subplot(2, 2, 2);
plot(1:length(ConfCase), y(:, 1), 'g', 'LineWidth', 2);
title('Susceptible');
xlabel('Days');
ylabel('Population');

subplot(2, 2, 3);
plot(1:length(ConfCase), y(:, 2), 'm', 'LineWidth', 2);
title('Infectious');
xlabel('Days');
ylabel('Population');

subplot(2, 2, 4);
plot(1:length(ConfCase), y(:, 3), 'c', 'LineWidth', 2);
title('Recovered');
xlabel('Days');
ylabel('Population');

% Define the SIV model
function dx = SIV(t, x, beta, delta, alpha)
    dx = zeros(4, 1);

    dx(1) = -beta * x(1) * x(2);
    dx(2) = beta * x(1) * x(2) - delta * x(2);
    dx(3) = alpha * x(2);
    dx(4) = beta * x(1) * x(2);
end

function value = SIV_RUN_ODE45(z)
    % LA County COVID-19 data
    global ConfCase TotalPopulation
    global TotalPopulation
    
    ten_dayCumulative = 5370;
    CumulativeTo10DayStart = 8453;

    % if length(ConfCase) < 1
    %     error('ConfCase array is empty or contains insufficient data.');
    % end

    I0 = ten_dayCumulative;
    R0 = 0.95 * 20 * CumulativeTo10DayStart;
    CI0 = ConfCase(1);
    S0 = TotalPopulation - I0 - R0;

    initialvalues = [S0; I0; R0; CI0];
    
    tspan = 1:length(ConfCase);

    [T, y] = ode45(@SIV, tspan, initialvalues, [], z(1), z(2),z(3));

    CI = y(:, 4);
    diff = CI - ConfCase.';
    value = norm(diff, 2);
end