%%
% read the data
data=readtable('cases_wastewater_vaccine.csv');
load('distance_small.mat');
load('info.mat');%column2 - zipcodes, column3 - income, column4 - population, column5 - age
load('pop.mat');
load('sanfrancisco_new.mat');
sanfrancisco = flipud(sanfrancisco);
load('marina.mat');
pop = caseswastewatervaccine(1:2,:);
marina = flipud(marina);
% info = info(2:7,:);
% info([2,5],:);
info = info(2:3,:);
Incomes = info(:,3);
n_values = info(:,4);
Ages = info(:,5); 

% Extract relevant data from sanfrancisco and marina
cases_sf = sanfrancisco(:, 1);
viralRNA_sf = sanfrancisco(:, 4);
tspan_sf = 10:size(sanfrancisco, 1);

cases_marina = marina(:, 1);
viralRNA_marina = marina(:, 4);
tspan_marina = 10:size(marina, 1);
%%
% Define parameters       
tspan = 10:size(sanfrancisco, 1);             

beta_star = 1.942E-6;  % beta value for the general population - estimated earlier
totalpop = 56503; % population of sf + marina

sigma = 1.0028E-5;  % rate at which recovered people lose immunity. 
k = 0.25;  % rate at which exposed become infected.
alpha = 15000; % shedding rate of virus per infected person
delta = 0.38; % recovery rate of infected people
epsilon = 0.05; % death rate
h = 0.1; % recovery rate of exposed people
k_p = 0.000439453125; % used in population function
A = 60.5756327138564;% parameter for population density dependence function
B = 3.561063241554381;% parameter for population density dependence function

% Calculate betas
betas = calculate_betas(n_values, totalpop,distance1,beta_star,A,B,k_p,pop);
%%
% Define initial conditions for three compartments (S, E, I, R, V)
I1 = marina(10,1) - marina(1,1);  %new cases in first 10 days of data = initial infected
E1 = .40 * (marina(10,1) - marina(1,1));  % initial exposed
R1 = .95 * (marina(1,1));  % initial recovered
Vnan1 = find(isnan(marina(:,4)));  % find which virus data is missing to check
V1 = marina(10,4);  
S1 = 25015 - I1 - R1 - E1;
CI1 = marina(10,1);

I2=sanfrancisco(10,1) - sanfrancisco(1,1);            
E2=.40*(sanfrancisco(10,1) - sanfrancisco(1,1));          
R2=.95*(sanfrancisco(1,1));
Vnan2 = find(isnan(sanfrancisco(:,4))); 
V2= sanfrancisco(10,4); 
S2=31488-I2-R2-E2;
CI2 = sanfrancisco(10,1);

% Combine initial conditions into a single vector
initial_conditions = [
    S1;E1;I1;R1;V1;CI1;
    S2;E2;I2;R2;V2;CI2
];

%%
% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigma, k, alpha, h, delta, epsilon);
%%
% Extract compartments from the solution for Marina
num_compartments = length(initial_conditions) / 6;
S_marina = y(:, 1);
E_marina = y(:, 2);
I_marina = y(:, 3);
R_marina = y(:, 4);
V_marina = y(:, 5);
CI_marina = y(:, 6);

% Extract compartments from the solution for SF
S_sf = y(:, 6 + 1);
E_sf = y(:, 6 + 2);
I_sf = y(:, 6 + 3);
R_sf = y(:, 6 + 4);
V_sf = y(:,  6 + 5);
CI_sf = y(:,  6 + 6);

% Plotting San Francisco data
figure;
subplot(2, 1, 1);
plot(tspan, CI_sf);
hold on;
scatter(tspan_sf, cases_sf(10:end), 'r.'); % Scatter plot for cases in CI (San Francisco)
hold off;
title('San Francisco - Compartment CI over time with Cases');
xlabel('Time');
ylabel('CI');

subplot(2, 1, 2);
daily_V_sf = [V_sf(1); diff(y(:, 5))];
semilogy(tspan, daily_V_sf, 'b.'); % Scatter plot for daily change in V (San Francisco)
hold on;
semilogy(tspan_sf, viralRNA_sf(10:end), 'r.'); % Scatter plot for viralRNA in V (San Francisco)
hold off;
title('SF - Daily change in Compartment V over time');
xlabel('Time');
ylabel('Daily Change in V');

% Plotting Marina data
figure;
subplot(2, 1, 1);
plot(tspan, CI_marina);
hold on;
scatter(tspan_marina, cases_marina(10:end), 'r.'); % Scatter plot for cases in CI (Marina)
hold off;
title('Marina - Compartment CI over time with Cases');
xlabel('Time');
ylabel('CI');

subplot(2, 1, 2);
daily_V_marina = [V_marina(1); diff(y(:, 5))];
semilogy(tspan, daily_V_marina, 'b.'); % Scatter plot for daily change in V (Marina)
hold on;
semilogy(tspan_marina, viralRNA_marina(10:end), 'r.'); % Scatter plot for viralRNA in V (Marina)
hold off;
title('Marina - Daily change in Compartment V over time');
xlabel('Time');
ylabel('Daily Change in V');

function dx = SIV(t, x, betas, sigma, k, alpha, h, delta, epsilon) 
    % Determine the number of compartments based on the length of the state vector x
    num_compartments = length(x) / 6; 
    
    % Initialize the output vector for the derivatives
    dx = zeros(6 * num_compartments, 1);
    
    % Iterate through each compartment to calculate the derivatives
    for j = 1:num_compartments
        % Extract variables for the j-th compartment
        Sj = x((j-1)*6 + 1);
        Ej = x((j-1)*6 + 2);
        Ij = x((j-1)*6 + 3);
        Rj = x((j-1)*6 + 4);
        Vj = x((j-1)*6 + 5);
        CIj = x((j-1)*6 + 6);

        % Calculate Beta using the entire row and infectious compartments
        Infectious = x(3:6:end); % Consider cumulative infections and deaths compartments too
        Beta = betas .* Infectious';
        
        % Define the differential equations for each compartment
        dx((j-1)*6 + 1) = -Beta(j) * Sj + sigma * Rj;
        dx((j-1)*6 + 2) = Beta(j) * Sj - k * Ej - h * Ej;
        dx((j-1)*6 + 3) = k * Ej - delta * Ij - epsilon * Ij;
        dx((j-1)*6 + 4) = delta * Ij + h * Ej - sigma * Rj;
        dx((j-1)*6 + 5) = alpha * Ij;
        % dx((j-1)*6 + 6) = k * Ej;
        dx((j-1)*6 + 6) = Beta(j) * Sj; % Equation for cumulative infections
    end
end

function betas = calculate_betas(n_values, totalpop, distance, beta_star, A, B, k_p,pop)
    % Calculate the transmission rates (betas) between compartments
    
    % Determine the number of compartments based on the length of n_values
    num_compartments = length(n_values);
    
    % Initialize the betas matrix
    betas = zeros(num_compartments, num_compartments);
    
    % Loop to calculate betas matrix
    for i = 1:num_compartments
        for j = 1:num_compartments  
            distance_between_i_and_j = distance(i, j);
            
            % Calculate betas
            distance_term = 1 / (1 + distance_between_i_and_j);
            exponential_term = exp(-k_p * n_values(j));
            betas(i, j) = (beta_star/totalpop) * pop(j) * (A + B * exponential_term) * distance_term;
        end
    end
end