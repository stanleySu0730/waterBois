% Given SF County COVID-19 data
global ConfCase;
ConfCase = [13823 15165 16449 17567 18545 19159 19567 20460 21017 22522 23233 24262 24936 25699 26238 27866 28665 29526 30334 31241 31703 32269 33247 34552 35447 36324 37374];
global TotalPopulation;
TotalPopulation = 9651332;
global distances;
folder_path = 'C:\Users\stanl\MatLab\waterBois';
file_name1 = 'distance.csv';
file_path1 = fullfile (folder_path,file_name1);
distances = readtable(file_path1);
file_name2 = 'cases_wastewater_vaccine.csv';
file_path2 = fullfile (folder_path,file_name2);
data = readtable(file_path2);

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

% Define SIV model
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

function value = SIV_RUN_ODE45(z)
    % SF County COVID-19 data
    global ConfCase TotalPopulation
    
    ten_dayCumulative = 5370;
    CumulativeTo10DayStart = 8453;

    I0 = ten_dayCumulative;
    E0 = 20 * I0;
    R0 = 0.95 * 20 * CumulativeTo10DayStart;
    CI0 = ConfCase(1);
    S0 = TotalPopulation - I0 - R0;

    initialvalues = [S0; I0; zeros(length(z)-2, 1)]; 
    num_compartments = (length(z) - 2) / 8; % Calculate the number of compartments
    
    beta_star = 1.67391974668301e-08; 
    n_values = z(1:num_compartments); 
    constant = z(num_compartments + 1); 
    eta_vax = z(end - num_compartments + 1:end);

    %retrieve parameters

    %Input distance matrix
    
    betas = zeros(num_compartments, num_compartments);
    
    % Calculate betas
    for i = 1:num_compartments
        for j = 1:num_compartments
            % Access distance information from the distances matrix
            distance_between_i_and_j = distances(i+2, j+1);
            betas(i,j) = beta_star / n_values(j) * 1 / (1 + (distance_between_i_and_j) / constant) * eta_vax(j);
        end
    end

    % Calculate m function based on Income and Age
    Income = 10000; % Example value for Income
    Age = 30; % Example value for Age
    
    % Input parameters 
    delta = 0.2; 
    epsilon = 0.1;
    P = 100; 
    C = 50;
    
    m = delta * (1 - (epsilon * exp(-k_i * Income) + m) * (P / (C + exp(-k_A * Age))));

end