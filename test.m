% Define parameters
Ages = [28,29,30];
n_values = [10, 20, 30];  % population densities for each zip code
eta_vax = [0.8, 0.9, 0.7];    % rate of vaccination
distances = [0, 0.2, 0.3; 0.3, 0, 0.6; 0.2, 0, 0.6];   % distance matrix - should be symmetric
Incomes = [100, 150, 120];

constant = 0.5;  % divide distances by this when calculating betas
alphas = 0.1;
beta_star = 1.67391974668301e-08;  % beta value for the general population
sigmas = [0.1, 0.15, 0.12];
k = 0.05;
gammas = 0.1;
lambdas = [0.3, 0.4, 0.2];
h = 0.5;
deltas = 0.5;
epsilons = 0.5;
k_A = 0.00001;
k_I = 0.00001;
mu = 0.01;

% Calculate betas
betas = calculate_betas(n_values, constant, eta_vax, distances,beta_star);

% Define initial conditions for three compartments (S, E, I, R, V, CI, D)
    initial_conditions = [
        800; 10; 5; 100; 50; 0; 0;  % Initial conditions for compartment 1
        700; 20; 3; 50; 30; 0; 0;    % Initial conditions for compartment 2
        900; 5; 8; 80; 20; 0; 0;     % Initial conditions for compartment 3
    ];
%%
% Define time span
tspan = [0 10];

% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigmas, k, alphas, h, deltas, epsilons, k_A, k_I, Incomes, Ages, mu);
%%
% Extract compartments from the solution
num_compartments = length(initial_conditions) / 7;
S = y(:, 1:7:end);
E = y(:, 2:7:end);
I = y(:, 3:7:end);
R = y(:, 4:7:end);
V = y(:, 5:7:end);
CI = y(:, 6:7:end);
D = y(:, 7:7:end);

% Plotting
figure;
% Plotting compartments
subplot(3, 3, 1);
plot(t, S);
title('Compartment S over time');
xlabel('Time');
ylabel('S');

subplot(3, 3, 2);
plot(t, E);
title('Compartment E over time');
xlabel('Time');
ylabel('E');

subplot(3, 3, 3);
plot(t, I);
title('Compartment I over time');
xlabel('Time');
ylabel('I');

subplot(3, 3, 4);
plot(t, R);
title('Compartment R over time');
xlabel('Time');
ylabel('R');

subplot(3, 3, 5);
plot(t, V);
title('Compartment V over time');
xlabel('Time');
ylabel('V');

subplot(3, 3, 6);
plot(t, CI);
title('Compartment CI over time');
xlabel('Time');
ylabel('CI');

subplot(3, 3, 7);
plot(t, D);
title('Compartment D (deaths) over time');
xlabel('Time');
ylabel('D');

function dx = SIV(t, x, betas, sigmas, k, alphas, h, deltas, epsilons, k_A, k_I, Incomes, Ages,mu)
    num_compartments = length(x) / 7; % Update the number of compartments
    dx = zeros(7 * num_compartments, 1); % Update size of dx
    
    for j = 1:num_compartments
        Sj = x((j-1)*7 + 1);
        Ej = x((j-1)*7 + 2);
        Ij = x((j-1)*7 + 3);
        Rj = x((j-1)*7 + 4);
        Vj = x((j-1)*7 + 5);
        CIj = x((j-1)*7 + 6);
        Dj = x((j-1)*7 + 7); % New compartment for deaths
        
        % Extract the row of betas corresponding to the j-th compartment
        Beta_row = betas(j, :);

        % Calculate Beta using the entire row and infectious compartments
        Infectious = x(3:7:end); % Consider cumulative infections and deaths compartments too
        Beta = Beta_row * Infectious;

        Income = Incomes(j);
        Age = Ages(j);
        m = deltas * (1 - epsilons * 1 / (1 + exp(k_A * Age)) * exp(k_I * Income));
        
        dx((j-1)*7 + 1) = -sum(Beta .* x(3:7:end)) * Sj;
        dx((j-1)*7 + 2) = sum(Beta .* x(3:7:end)) * Sj - k * Ej - h * Ej;
        dx((j-1)*7 + 3) = k * Ej - m * Ij;
        dx((j-1)*7 + 4) = m * Ij + h * Ej - sigmas(j) * Rj;
        dx((j-1)*7 + 5) = alphas * Ij;
        dx((j-1)*7 + 6) = sum(Beta .* x(3:7:end)) * Sj + sigmas(j) * Rj; % Equation for cumulative infections
        dx((j-1)*7 + 7) = mu * Ij; % Equation for deaths
    end
end

function betas = calculate_betas(n_values, totalpop, distances, beta_star, A, B, k)
    % Calculate the transmission rates (betas) between compartments
    
    % Determine the number of compartments based on the length of n_values
    num_compartments = length(n_values);
    
    % Initialize the betas matrix
    betas = zeros(num_compartments, num_compartments);
    
    % Loop to calculate betas matrix
    for i = 1:num_compartments
        for j = 1:num_compartments  
            distance_between_i_and_j = distances(i, j);
            
            % Calculate betas
            pop_density_term = n_values(j);
            distance_term = 1 / (1 + distance_between_i_and_j);
            exponential_term = exp(-k * totalpop);
            betas(i, j) = beta_star * pop_density_term * (A + B * exponential_term) * distance_term;
        end
    end
end
