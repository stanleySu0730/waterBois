%%
load('sanfrancisco.mat');   
sanfrancisco = flipud(sanfrancisco);
%%
% Extract the first and fourth columns
cases = sanfrancisco(:, 1);
viralRNA = sanfrancisco(:, 4);
indices = 1:length(cases); % Creating indices for x-axis

% Plotting the first graph for cases
figure;
subplot(2, 1, 1); % Create subplot 1
scatter(indices, cases, 'b', 'filled'); % Plot cases
ylabel('Cases'); % Label for y-axis
xlabel('Index'); % Label for x-axis
title('Graph of Cases');

% Plotting the second graph for viralRNA as dots
subplot(2, 1, 2); % Create subplot 2 
scatter(indices, viralRNA, 'r', 'filled'); % Plot viralRNA
ylabel('Viral RNA'); % Label for y-axis
xlabel('Index'); % Label for x-axis
title('Graph of Viral RNA');