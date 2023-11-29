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

function dx = SIV(t, x, betas, sigmas, k, alphas, Incomes, gammas, lambdas)
    num_compartments = length(x) / 5;
    dx = zeros(5 * num_compartments, 1);
    
    for j = 1:num_compartments
        Sj = x((j-1)*5 + 1);
        Ej = x((j-1)*5 + 2);
        Ij = x((j-1)*5 + 3);
        Rj = x((j-1)*5 + 4);
        Vj = x((j-1)*5 + 5);
        
        sum_beta_I = 0;
        for i = 1:num_compartments
            sum_beta_I = sum_beta_I + betas(i, j) * x((i-1)*5 + 3);
        end
        
        dx((j-1)*5 + 1) = -sum_beta_I * Sj + sigmas(j) * Rj;
        dx((j-1)*5 + 2) = sum_beta_I * Sj - k * Ej;
        dx((j-1)*5 + 3) = k * Ej - alphas(j) * Incomes(j) * Ij;
        dx((j-1)*5 + 4) = alphas(j) * Incomes(j) * Ij - sigmas(j) * Rj;
        dx((j-1)*5 + 5) = gammas * (lambdas(1) * Ej + lambdas(2) * Ij + lambdas(3) * Rj);
    end
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

function value = SIV_RUN_ODE45(z, culmulativeCases, initialvalues, tspan, betas, sigmas, k, alphas, Incomes, gammas, lambdas)
    [~, y] = ode45(@(t, x) SIV(t, x, betas, sigmas, k, alphas, Incomes, gammas, lambdas), tspan, initialvalues);
    CI = y(:, 3);
    diff = CI - culmulativeCases;
    value = norm(diff, 2);
end
