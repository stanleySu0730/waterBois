%%
% read the data
data=readtable('cases_wastewater_vaccine.csv');
load('distance_small.mat');
load('pop.mat');
load('info.mat');%column2 - zipcodes, column3 - income, column4 - population, column5 - age
load('fremont.mat');
fremont = flipud(fremont);
load('newark.mat');
newark = flipud(newark);
% info = info(2:7,:);
% info([2,5],:);
info = info(4:5,:);
pop = caseswastewatervaccine(3:4,:);
Incomes = info(:,3);
n_values = info(:,4);
Ages = info(:,5); 

% Extract relevant data from sanfrancisco and marina
cases_fremont = fremont(:, 1);
viralRNA_fremont = fremont(:, 4);
tspan_fremont = 10:size(fremont, 1);

cases_newark = newark(:, 1);
viralRNA_newark = newark(:, 4);
tspan_newark = 10:size(newark, 1);
%%
% Define parameters       
tspan = 10:size(fremont, 1);             

beta_star = 4.042E-7;  % beta value for the general population - estimated earlier
for j = 1:size(pop)
    totalpop = totalpop + pop(j);
end

sigma = 1.0028E-5;  % rate at which recovered people lose immunity. 
k = 0.25;  % rate at which exposed become infected.
alpha = 100; % shedding rate of virus per infected person
delta = 0.38; % recovery rate of infected people
epsilon = 0.34; % death rate
h = 0.1; % recovery rate of exposed people
k_p = 0.000439453125; % used in population function
A = 60.5756327138564;% parameter for population density dependence function
B = 3.561063241554381;% parameter for population density dependence function

% Calculate betas
betas = calculate_betas(n_values, totalpop,distance1,beta_star,A,B,k_p,pop);

%%
% Define initial conditions for three compartments (S, E, I, R, V)
I1 = fremont(10,1) - fremont(1,1);
E1 = .40 * (fremont(10,1) - fremont(1,1));
R1 = .95 * (fremont(1,1));
V1 = fremont(10,4);
S1 = 63911 - I1 - R1 - E1;
CI1 = fremont(10,1);

I2=newark(10,1) - newark(10,1);            
E2=.40*(newark(10,1) - newark(1,1));          
R2=.95*( newark(1,1));
V2= newark(10,4); 
S2=44245-I2-R2-E2;
CI2 = newark(10,1);

% Combine initial conditions into a single vector
initial_conditions = [
    S1;E1;I1;R1;V1;CI1;
    S2;E2;I2;R2;V2;CI2
];

%%
% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigma, k, alpha, h, delta, epsilon);
%%
% Extract compartments from the solution for San Francisco
num_compartments = length(initial_conditions) / 6;
S_marina = y(:, 1:6);
E_marina = y(:, 2:6);
I_marina = y(:, 3:6);
R_marina = y(:, 4:6);
V_marina = y(:, 5:6);
CI_marina = y(:, 6:6);

% Extract compartments from the solution for Marina
S_sf = y(:, 6 + 1);
E_sf = y(:, 6 + 2);
I_sf = y(:, 6 + 3);
R_sf = y(:, 6 + 4);
V_sf = y(:, 6 + 5);
CI_sf = y(:, 6 + 6);

% Plotting San Francisco data
figure;
subplot(2, 1, 1);
plot(tspan, CI_sf);
hold on;
scatter(tspan_newark, cases_newark(10:end), 'r.'); % Scatter plot for cases in CI (San Francisco)
hold off;
title('Newark - Compartment CI over time with Cases');
xlabel('Time');
ylabel('CI');

subplot(2, 1, 2);
plot(tspan, V_sf);
hold on;
semilogy(tspan_newark, viralRNA_newark(10:end), 'b.'); % Scatter plot for viralRNA in V (San Francisco)
hold off;
title('Newark - Compartment V over time with Viral RNA');
xlabel('Time');
ylabel('V');

% Plotting Marina data
figure;
subplot(2, 1, 1);
plot(tspan, CI_marina);
hold on;
scatter(tspan_fremont, cases_fremont(10:end), 'r.'); % Scatter plot for cases in CI (Marina)
hold off;
title('Fremont - Compartment CI over time with Cases');
xlabel('Time');
ylabel('CI');

subplot(2, 1, 2);
plot(tspan, V_marina);
hold on;
semilogy(tspan_fremont, viralRNA_fremont(10:end), 'b.'); % Scatter plot for viralRNA in V (Marina)
hold off;
title('Fremont - Compartment V over time with Viral RNA');
xlabel('Time');
ylabel('V');

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
            exponential_term = exp(-k_p * totalpop);
            betas(i, j) = (beta_star/pop(i)) * pop(j) * (A + B * exponential_term) * distance_term;
        end
    end
end