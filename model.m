global distances;

% Parameter estimation using fmincon
LowerBounds = [1.00e-9, 1.00e-9, 1.00e-9];
UpperBounds = [1.00e-7, 1.00e-3, 1.00e-3];
xstart = 0.5 * (LowerBounds + UpperBounds);

problem = createOptimProblem('fmincon', 'objective', @SIV_RUN_ODE45, 'x0', xstart, 'lb', LowerBounds, 'ub', UpperBounds);
problem.options = optimoptions('fmincon', 'MaxFunEvals', 9999, 'MaxIter', 9999);

numstartpoints = 10;
ms = MultiStart('UseParallel', true, 'Display', 'iter');
[b, fval, exitflag, output, manymins] = run(ms, problem, numstartpoints);
for i = 1:length(incomeArray)
    % Get relevant values from arrays 
    Income = incomeArray(i);
    Age = ageArray(i);

    % Input parameters
    delta = 0.2;
    epsilon = 0.1;
    k_i = 0.05;
    k_A = 0.02;
    P = 100;
    C = 50;

    m = delta * (1 - (epsilon * exp(-k_i * Income) + m) * (P / (C + exp(-k_A * Age))));

    squared_errors = calculate_squared_errors(z, ConfCase, TotalPopulation, distances);
    betas = calculate_betas(z, distances);
end

function betas = calculate_betas(z, distances)
    num_compartments = (length(z) - 2) / 8;
    
    beta_star = 1.67391974668301e-08;
    n_values = z(1:num_compartments);
    constant = z(num_compartments + 1);
    eta_vax = z(end - num_compartments + 1:end);
    betas = zeros(num_compartments, num_compartments);
    
    % Loop to calculate betas matrix
    for i = 1:num_compartments
        for j = 1:num_compartments
            distance_between_i_and_j = distances(i+2, j+1);
            % Calculate betas based on the given formula
            betas(i,j) = beta_star / n_values(j) * 1 / (1 + (distance_between_i_and_j) / constant) * eta_vax(j);
        end
    end
end

function SIV_RUN_ODE45(betas, ConfCase, TotalPopulation)
    % ten_dayCumulative = ...
    % CumulativeTo10DayStart = ...
    I0 = ten_dayCumulative;
    E0 = 20 * I0;
    R0 = 0.95 * 20 * CumulativeTo10DayStart;
    CI0 = ConfCase(1);
    S0 = TotalPopulation - I0 - R0;

    % Time span
    t_span = 1:length(ConfCase);

    initialvalues = [S0; I0; zeros(length(z)-2, 1)]; 

    % Solve the ODE system
    [t, y] = ode45(@(t, x) SIV(t, x, betas), t_span, initialvalues);

    % Calculate squared errors
    % squared_errors = sum((model_cases - ConfCase').^2);
end
