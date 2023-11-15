close all;
clear all;

rng("default") %fix rng seed for consistency.

%% Solving model and generate synthetic data for viral RNA measurement in WW
xinit = [2e6, 1, 0 , 1]; %initial value
h = 0.01; %time step
T = 100; %last time point
time = 0:h:T; %time vector

[~,truth] = ode45(@(t,x) SIV(t,x),time,xinit); %solve the model and store variables in "truth"

% The "truth" matrix stores state variables over time based on the prescribed time vector "time".
% However, because of the step size "h", the time vector contains state measurement 'almost continuously'.
% So we need to figure out what the daily measurements based on this continuous dynamic variable.

% First, set up a variable holder for the daily "truth".
daily_truth_ww = zeros(1,T); % 

for ii = 1:T
    % Because of the time step h, every day contains exactly 100 elements. 
    % This means the daily measurement on day 1 is just the recorded number on day 1 - day 0.
    % In other words, the daily measurement on day 1 would be truth(101) - truth(1).
    % We only want the third column truth(:,3), which only contains information on viral RNA in wastewater.
    daily_truth_ww(ii) = truth(100*(ii)+1,3) - truth(100*(ii-1)+1,3);
end

% Add some errors to the truth daily observation "daily_truth" using randn (normal distribution).
daily_ww_obs = daily_truth_ww + randn(size(daily_truth_ww)).*daily_truth_ww/2; %the "daily_truth_ww/2 is just a scaling. You can try different scaling.
tspan = 1:1:length(daily_ww_obs); %just getting a different (daily) time vector for plotting.

figure(1); hold on; box on;
% Plot every 3 points. Change the number in the middle (1: number : end) to increase the frequency for plotting.
scatter(tspan(1:3:end),log10(daily_ww_obs(1:3:end)),30,'r'); %scatter plot of the daily observation
plot(tspan(1:3:end),log10(daily_ww_obs(1:3:end)),':r'); %connecting the scatter plot

% Plot the daily truth measurements (no error) in the same frequency.
plot(tspan(1:3:end),log10(daily_truth_ww(1:3:end)),'r');

% axis and legend stuffs
ylim([5 12])
xlim([0 tspan(end)])
xlabel('Days','FontSize',14)
ylabel('Log10 daily viral RNA copies in Wastewater','FontSize',14)

legend('Synthetic WW Data','','True WW Dynamics','Location','NorthEast','FontSize',16)

%% Generate synthetic data for the number of infectious individuals
% This is similar to above, but is done for the number of infectious individuals.
% Note that in the basic SI-V model, you will need an additional auxiliary compartment C, to track the cumulative number of infected individuals over one entire outbreak. 
% The daily number of cases is then the difference of C over each day. 
% For example, the number of cases on day 1, is C(day 1) - C(day 0).
daily_truth_case = zeros(1,T);

for ii = 1:T
    daily_truth_case(ii) = truth(100*(ii)+1,4) - truth(100*(ii-1)+1,4);
end

daily_case_obs = daily_truth_case + randn(size(daily_truth_case)).*daily_truth_case/2;

figure(2); hold on; box on;
scatter(tspan(1:3:end),log10(daily_case_obs(1:3:end)),30,'r');
plot(tspan(1:3:end),log10(daily_case_obs(1:3:end)),':r');

plot(tspan,log10(daily_truth_case),'r');

% ylim([5 12])
xlim([0 tspan(end)])
xlabel('Days','FontSize',14)
ylabel('Log10 daily case','FontSize',14)

legend('Synthetic Case Data','','True Case Dynamics','Location','NorthEast','FontSize',16)
hold off;

%% Corelation - with true dynamics (no error)
% first with the entire true dynamics - just to get an idea of the shape.
% The code below is one standard way to get the R-score
    y = daily_truth_case'; 
    x = daily_truth_ww';
    X = [ones(length(x),1) x];
    b = X\y;

    yCalc2 = X*b;

    figure(3); box on; hold on;
    plot(x,yCalc2,'r')
    scatter(x,y,30,'r')
    axis tight

    ylabel('Daily true case','FontSize',14);
    xlabel('Daily true ww','FontSize',14)

    %calculate R2
    Rsq2 = 1 - sum((y-yCalc2).^2)/sum((y-mean(y)).^2);
    R = Rsq2^(1/2);

    legend(strcat('Correlation R =',{' '},num2str(R)),'FontSize',16,'Location', 'NorthEast')

% Correlation on increasing data set

    segment_perc = 0.02:0.01:1;
    R_vec = nan(length(segment_perc),1);
    data_size = length(daily_truth_case);

    for ii = 1:length(segment_perc)

        y = daily_truth_case(1:floor(segment_perc(ii)*data_size))'; 
        x = daily_truth_ww(1:floor(segment_perc(ii)*data_size))';
        X = [ones(length(x),1) x];
        b = X\y;

        yCalc2 = X*b;

        %calculate R2
        Rsq2 = 1 - sum((y-yCalc2).^2)/sum((y-mean(y)).^2);
        R_vec(ii) = Rsq2^(1/2);

    end

    figure(4); box on; hold on;
    scatter(segment_perc, R_vec,30,'r');
    plot(segment_perc, R_vec,'r');
    xlabel('Fraction of total data usage','FontSize',14)
    ylabel('Correlation R - true dynamics','FontSize',14)

%% Corelation - with synthetic data (with error)
% Similar to above, but the correlation is done with synthetic data.

synthetic_case_data = daily_case_obs; %(1:1:end);
synthetic_ww_data = daily_ww_obs;%(1:1:end);

% Correlation on increasing data set

    segment_perc = 0.02:0.01:1;
    R_vec = nan(length(segment_perc),1);
    data_size = length(synthetic_case_data);

    for ii = 1:length(segment_perc)

        y = synthetic_case_data(1:floor(segment_perc(ii)*data_size))'; 
        x = synthetic_ww_data(1:floor(segment_perc(ii)*data_size))';
        X = [ones(length(x),1) x];
        b = X\y;

        yCalc2 = X*b;

        %calculate R2
        Rsq2 = 1 - sum((y-yCalc2).^2)/sum((y-mean(y)).^2);
        R_vec(ii) = Rsq2^(1/2);

    end

    figure(5); box on; hold on;
    scatter(segment_perc, R_vec,30,'r');
    plot(segment_perc, R_vec,'r');
    xlabel('Fraction of total data usage','FontSize',14)
    ylabel('Correlation R - synthetic data','FontSize',14)

% what it looks like with the entire curve of synthetic data

    y = synthetic_case_data'; 
    x = synthetic_ww_data';
    X = [ones(length(x),1) x];
    b = X\y;

    yCalc2 = X*b;

    figure(6); box on; hold on;
    plot(x,yCalc2,'r')
    scatter(x,y,30,'r')
    axis tight

    ylabel('Daily synthetic case','FontSize',14);
    xlabel('Daily synthetic ww','FontSize',14)

    %calculate R2
    Rsq2 = 1 - sum((y-yCalc2).^2)/sum((y-mean(y)).^2);
    R = Rsq2^(1/2);

    legend(strcat('Correlation R =',{' '},num2str(R)),'FontSize',16,'Location', 'NorthEast')

    xlim([0 inf])
    ylim([0 inf])

%%
function dx = SIV(t,x)
% A simple SI-V model for wastewater
% S' = -beta*S*I
% I' = beta*S*I - delta*I
% V' = alpha*I
% C' = beta*S*I (this auxiliary variable is used to track cumulative cases).

beta = 5e-7; delta = 1/8; alpha = 1e5; %just some reasonable model parameters
dx  = zeros(4,size(x,2));

dx(1,:) = -beta.*x(1,:).*x(2,:);
dx(2,:) = beta.*x(1,:).*x(2,:) - delta*x(2,:);
dx(3,:) = alpha*x(2,:);
dx(4,:) = beta.*x(1,:).*x(2,:);

end
