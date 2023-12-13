%%
% read the data
data=readtable('cases_wastewater_vaccine.csv');
load('distances.mat');
load('info.mat');%column2 - zipcodes, column3 - income, column4 - population, column5 - age
load('sanfrancisco.mat');
info = info(2:7,:);
% info([2,5],:);
Incomes = info(:,3);
n_values = info(:,4);
Ages = info(:,5);
%%
% Define parameters       
tspan=[0, 20];              

beta_star = 1.04249253978453E-07;  % beta value for the general population

sigmas = 0.252323055009390;
k = 1.00280055177554E-5;
alphas = 22.106129320812855;
h = 0.25;
deltas = 0.200518657975840;
epsilons = 0.0822347534479802;
k_A = 9.99983721742511E-7;
k_I = 1.0085048810571665E-12;
mu = 0.01; %unknown

% Calculate betas
betas = calculate_betas(n_values, constant, eta_vax, distances,beta_star);

%%
% Define initial conditions for three compartments (S, E, I, R, V)
I1=38476-38323;            
E1=.40*(38323-38153);          
R1=.95*(38323-38153);
V1=20465995.29; 
S1=79997-I1-R1-E1;

I2=4866-4830;            
E2=.40*(4830-4783);          
R2=.95*(4830-4783);
V2=182255.4581; 
S2=25015-I2-R2-E2;

% Combine initial conditions into a single vector
initial_conditions = [
    S1;E1;I1;R1;V1;0;0;
    S2;E2;I2;R2;V2;0;0
];

%%
% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigmas, k, alphas, h, deltas, epsilons, k_A, k_I, Incomes, Ages,mu);
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
title('Compartment D over time');
xlabel('Time');
ylabel('D');

function dx = SIV(t, x, betas, sigmas, k, alphas, h, deltas, epsilons, k_A, k_I, Incomes, Ages, mu) 
    % Determine the number of compartments based on the length of the state vector x
    num_compartments = length(x) / 7; 
    
    % Initialize the output vector for the derivatives
    dx = zeros(7 * num_compartments, 1);
    
    % Iterate through each compartment to calculate the derivatives
    for j = 1:num_compartments
        % Extract variables for the j-th compartment
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

        % Calculate parameters related to the j-th compartment
        Income = Incomes(j);
        Age = Ages(j);
        m = deltas * (1 - epsilons * 1 / (1 + exp(k_A * Age)) * exp(k_I * Income));
        
        % Define the differential equn_valuations for each compartment
        dx((j-1)*7 + 1) = -sum(Beta .* x(3:7:end)) * Sj + sigmas * Rj;
        dx((j-1)*7 + 2) = sum(Beta .* x(3:7:end)) * Sj - k * Ej - h * Ej;
        dx((j-1)*7 + 3) = k * Ej - m * Ij - mu * Ij;
        dx((j-1)*7 + 4) = m * Ij + h * Ej - sigmas * Rj;
        dx((j-1)*7 + 5) = alphas * Ij;
        dx((j-1)*7 + 6) = sum(Beta .* x(3:7:end)) * Sj; % Equation for cumulative infections
        dx((j-1)*7 + 7) = mu * Ij; % Equation for deaths
    end
end

function betas = calculate_betas(n_values, totalpop, distances, beta_star, A, B, k_p)
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
            exponential_term = exp(-k_p * totalpop);
            betas(i, j) = beta_star * pop_density_term * (A + B * exponential_term) * distance_term;
        end
    end
end