% Define parameters
n_values = [10, 20, 30];  % population densities for each zip code
constant = 0.5;  % divide distances by this when calculating betas
eta_vax = [0.8, 0.9, 0.7];    % rate of vaccination
distances = [0.1, 0.2, 0.3; 0.4, 0.5, 0.6; 0.7, 0.8, 0.9];   % distance matrix - should be symmetric
beta_star = 1.67391974668301e-08;  % beta value for the general population

sigmas = [0.1, 0.15, 0.12];
k = 0.05;
alphas = [0.2, 0.25, 0.18];
Incomes = [100, 150, 120];
gammas = 0.1;
lambdas = [0.3, 0.4, 0.2];

% Calculate betas
betas = calculate_betas(n_values, constant, eta_vax, distances,beta_star);

%%
% Define initial conditions for three compartments (S, E, I, R, V)
initial_S1 = 800;
initial_E1 = 10;
initial_I1 = 5;
initial_R1 = 100;
initial_V1 = 50;

initial_S2 = 700;
initial_E2 = 20;
initial_I2 = 3;
initial_R2 = 50;
initial_V2 = 30;

initial_S3 = 900;
initial_E3 = 5;
initial_I3 = 8;
initial_R3 = 80;
initial_V3 = 20;

% Combine initial conditions into a single vector
initial_conditions = [
    initial_S1; initial_E1; initial_I1; initial_R1; initial_V1;
    initial_S2; initial_E2; initial_I2; initial_R2; initial_V2;
    initial_S3; initial_E3; initial_I3; initial_R3; initial_V3
];

%%
% Define time span
tspan = [0 10];

% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions,[], betas, sigmas, k, alphas, Incomes, gammas, lambdas);
%%
% Extract compartments from the solution
num_compartments = length(initial_conditions) / 5;
S = y(:, 1:5:end);
E = y(:, 2:5:end);
I = y(:, 3:5:end);
R = y(:, 4:5:end);
V = y(:, 5:5:end);

% Plotting
figure;

% Plot for compartment S
subplot(3, 2, 1);
plot(t, S);
title('Compartment S over time');
xlabel('Time');
ylabel('S');

% Plot for compartment E
subplot(3, 2, 2);
plot(t, E);
title('Compartment E over time');
xlabel('Time');
ylabel('E');

% Plot for compartment I
subplot(3, 2, 3);
plot(t, I);
title('Compartment I over time');
xlabel('Time');
ylabel('I');

% Plot for compartment R
subplot(3, 2, 4);
plot(t, R);
title('Compartment R over time');
xlabel('Time');
ylabel('R');

% Plot for compartment V
subplot(3, 2, 5);
plot(t, V);
title('Compartment V over time');
xlabel('Time');
ylabel('V');

function dx = SIV(t, x, betas, sigmas, k, alphas, Incomes, gammas, lambdas)
    num_compartments = length(x) / 5;
    dx = zeros(5 * num_compartments, 1);
    
    for j = 1:num_compartments
        Sj = x((j-1)*5 + 1);
        Ej = x((j-1)*5 + 2);
        Ij = x((j-1)*5 + 3);
        Rj = x((j-1)*5 + 4);
        Vj = x((j-1)*5 + 5);
        
        Beta = betas * reshape(x(3:5:end), num_compartments, 1);

        dx((j-1)*5 + 1) = -Beta(j) * Sj + sigmas(j) * Rj;
        dx((j-1)*5 + 2) = Beta(j) * Sj - k * Ej;
        dx((j-1)*5 + 3) = k * Ej - alphas(j) * Incomes(j) * Ij;
        dx((j-1)*5 + 4) = alphas(j) * Incomes(j) * Ij - sigmas(j) * Rj;
        dx((j-1)*5 + 5) = gammas * (lambdas(1) * Ej + lambdas(2) * Ij + lambdas(3) * Rj);
    end
end

function betas = calculate_betas(n_values, constant, eta_vax, distances, beta_star)
    num_compartments = length(n_values);
    betas = zeros(num_compartments, num_compartments);
    
    % Loop to calculate betas matrix
    for i = 1:num_compartments
        for j = 1:num_compartments  
            distance_between_i_and_j = distances(i, j); % Assuming the distance matrix is properly formatted
            
            % Calculate betas based on the given formula for transmission matrix
            betas(i, j) = beta_star / n_values(j) * 1 / (1 + (distance_between_i_and_j) / constant) * eta_vax(j);
        end
    end
end