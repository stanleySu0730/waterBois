%%
% read the data
data=readtable('cases_wastewater_vaccine.csv');
load('distance.mat');
load('info.mat');%column2 - zipcodes, column3 - income, column4 - population, column5 - age
load('sanfrancisco.mat');
sanfrancisco = flipud(sanfrancisco);
load('fremont.mat');
fremont = flipud(fremont);
load('marina.mat');
marina = flipud(marina);
load('newark.mat');
newark = flipud(newark);
load('sanleandro.mat');
sanleandro = flipud(sanleandro);
load('unioncity.mat');
unioncity = flipud(unioncity);
load('pop.mat');
info = info(2:7,:);
pop = caseswastewatervaccine(1:end,:);
% info([2,5],:);
Incomes = info(:,3);
n_values = info(:,4);
Ages = info(:,5);
%%
% Define parameters       
tspan = 10:size(fremont, 1);             

beta_star = 3.842E-7;  % beta value for the general population - estimated earlier
totalpop = 0;
for j = 1:size(pop)
    totalpop = totalpop + pop(j);
end

sigma = 1.0028E-5;  % rate at which recovered people lose immunity. 
k = 0.25;  % rate at which exposed become infected.
alpha = 25000; % shedding rate of virus per infected person
delta = 0.38; % recovery rate of infected people
epsilon = 0.05; % death rate
h = 0.1; % recovery rate of exposed people
k_p = 0.000439453125; % used in population function
A = 60.5756327138564;% parameter for population density dependence function
B = 3.561063241554381;% parameter for population density dependence function

% Calculate betas
betas = calculate_betas(n_values, totalpop,distance,beta_star,A,B,k_p,pop);

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

I3 = fremont(10,1) - fremont(1,1);
E3 = .40 * (fremont(10,1) - fremont(1,1));
R3 = .95 * (fremont(1,1));
Vnan3 = find(isnan(fremont(:,4)));
V3 = fremont(10,4);
S3 = 63911 - I3 - R3 - E3;
CI3 = fremont(10,1);

I4=newark(10,1) - newark(10,1);            
E4=.40*(newark(10,1) - newark(1,1));          
R4=.95*( newark(1,1));
Vnan4 = find(isnan(newark(:,4)));
V4= newark(10,4); 
S4=44245-I4-R4-E4;
CI4 = newark(10,1);

I5 = sanleandro(10,1) - sanleandro(1,1);
E5 = .40 * (sanleandro(10,1) - sanleandro(1,1));
R5 = .95 * (fremont(1,1));
Vnan5 = find(isnan(sanleandro(:,4)));
V5 = sanleandro(10,4);
S5 = 48666 - I5 - R5 - E5;
CI5 = sanleandro(10,1);

I6=unioncity(10,1) - unioncity(10,1);            
E6=.40*(unioncity(10,1) - unioncity(1,1));          
R6=.95*( unioncity(1,1));
Vnan6 = find(isnan(unioncity(:,4)));
V6= unioncity(10,4); 
S6=72446-I6-R6-E6;
CI6 = unioncity(10,1);

% Combine initial conditions into a single vector
initial_conditions = [
    S1;E1;I1;R1;V1;CI1;
    S2;E2;I2;R2;V2;CI2;
    S3;E3;I3;R3;V3;CI3;
    S4;E4;I4;R4;V4;CI4;
    S5;E5;I5;R5;V5;CI5;
    S6;E6;I6;R6;V6;CI6;
];

%%
% Solve the ODE
[t, y] = ode45(@SIV, tspan, initial_conditions, [], betas, sigma, k, alpha, h, delta, epsilon);
%%
% Extract compartments from the solution
num_compartments = length(initial_conditions) / 6;
S = y(:, 1:6:end);
E = y(:, 2:6:end);
I = y(:, 3:6:end);
R = y(:, 4:6:end);
V = y(:, 5:6:end);
CI = y(:, 6:6:end);

% Calculate Daily_V for each region
daily_V = zeros(size(V));
for i = 1:size(V, 2)
    daily_V(1, i) = V(1, i); % Set the first value of daily_V as V(1)
    for j = 2:length(V)
        daily_V(j, i) = V(j, i) - V(j-1, i); % Calculate the difference for subsequent values
    end
end
% Plotting SEIRV and CI compartments for each region
regions = {'Marina', 'San Francisco', 'Fremont', 'Newark', 'San Leandro', 'Union City'};

% Create a new figure for SEIRV compartments
figure;

% Plotting S, E, I, R, V, CI compartments for all regions on a single graph
subplot(3, 2, 1);
hold on;
for i = 1:6
    plot(t, S(:, i) ./ pop(i), 'DisplayName', regions{i});
end
hold off;
title('S Compartment for All Regions');
xlabel('Time');
ylabel('Susceptible / Population');
legend('Location', 'best');

subplot(3, 2, 2);
hold on;
for i = 1:6
    plot(t, E(:, i) ./ pop(i), 'DisplayName', regions{i});
end
hold off;
title('E Compartment for All Regions');
xlabel('Time');
ylabel('Exposed / Population');
legend('Location', 'best');

subplot(3, 2, 3);
hold on;
for i = 1:6
    plot(t, I(:, i) ./ pop(i), 'DisplayName', regions{i});
end
hold off;
title('I Compartment for All Regions');
xlabel('Time');
ylabel('Infected / Population');
legend('Location', 'best');

subplot(3, 2, 4);
hold on;
for i = 1:6
    plot(t, R(:, i) ./ pop(i), 'DisplayName', regions{i});
end
hold off;
title('R Compartment for All Regions');
xlabel('Time');
ylabel('Recovered / Population');
legend('Location', 'best');

subplot(3, 2, 5);
hold on;
for i = 1:6
    semilogy(t, abs(daily_V(:, i)), '.', 'DisplayName', regions{i});
end
hold off;
title('V Compartment for All Regions');
xlabel('Time');
ylabel('Population');
legend('Location', 'best');

subplot(3, 2, 6);
hold on;
for i = 1:6
    plot(t, CI(:, i) ./ pop(i), 'DisplayName', regions{i});
end
hold off;
title('CI Compartment for All Regions');
xlabel('Time');
ylabel('Cumulative Infections / Population');
legend('Location', 'best');

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
            exponential_term = exp(-k_p * n_values(i));
            betas(i, j) = (beta_star/totalpop) * pop(j) * (A + B * exponential_term) * distance_term;
        end
    end
end