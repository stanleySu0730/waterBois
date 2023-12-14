% read the data
data=readtable('cases_wastewater_vaccine.csv');
load('info.mat');%column2 - zipcodes, column3 - income, column4 - population, column5 - age
load('pop.mat');
load('sanfrancisco_new.mat');
sanfrancisco = flipud(sanfrancisco);
pop = caseswastewatervaccine(1,:);
% info = info(2:7,:);
% info([2,5],:);
info = info(2,:);
Incomes = info(:,3);
n_values = info(:,4);
Ages = info(:,5); 

% Extract relevant data from sanfrancisco and marina
cases_sf = sanfrancisco(:, 1);
viralRNA_sf = sanfrancisco(:, 4);
tspan_sf = 10:size(sanfrancisco, 1);
%%
% Define parameters       
tspan = 10:size(sanfrancisco, 1);             

beta_star = 1.942E-6;  % beta value for the general population - estimated earlier
totalpop = 56503; % population of sf + marina
distance = 0;
sigma = 1.0028E-5;  % rate at which recovered people lose immunity. 
k = 0.25;  % rate at which exposed become infected.
alpha = 100000; % shedding rate of virus per infected person
delta = 0.38; % recovery rate of infected people
epsilon = 0.05; % death rate
h = 0.1; % recovery rate of exposed people
k_p = 0.000439; % used in population function
A = 2500;% parameter for population density dependence function
B = 3;% parameter for population density dependence function

% Calculate betas
betas = calculate_betas(n_values, totalpop,distance,beta_star,A,B,k_p,pop);
%%
% Define initial conditions for three compartments (S, E, I, R, V)
V2= sanfrancisco(1,4); 
I2=V2/alpha;            
E2=.40*I2;          
R2=.95*(sanfrancisco(1,1));
Vnan2 = find(isnan(sanfrancisco(:,4))); 
S2=31488-I2-R2-E2;
CI2 = sanfrancisco(1,1);

% Combine initial conditions into a single vector
initial_conditions = [
    S2;E2;I2;R2;V2;CI2
];

%%
% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigma, k, alpha, h, delta, epsilon);
%%
% Extract compartments from the solution for Marina
num_compartments = length(initial_conditions) / 6;
S_sf = y(:, 1);
E_sf = y(:, 2);
I_sf = y(:, 3);
R_sf = y(:, 4);
V_sf = y(:, 5);
CI_sf = y(:, 6);

% Plotting San Francisco data
figure;
subplot(2, 1, 1);
plot(tspan, CI_sf);
hold on;
scatter(tspan, cases_sf(10:end), 'r.'); % Scatter plot for cases in CI (San Francisco)
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

function betas = calculate_betas(n_values, totalpop, distance, beta_star, A, B, k_p, pop)
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
            exponential_term = exp(-k_p * n_values(i));
            betas(i, j) = (beta_star / totalpop) * pop(j) * (A(i) + B(i) * exponential_term) * distance_term;
        end
    end
end