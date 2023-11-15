% Objective function for parameter estimation
function value = SIV_RUN_ODE45(z)
    % LA County COVID-19 data

    day = length(ConfCase);
    tspan = 1:1:day;

    I0 = ten_dayCumulative;
    R0 = 0.95 * 20 * CumulativeTo10DayStart;
    CI0 = ConfCase(1);
    S0 = TotalPopulation - I0 - R0;

    initialvalues = [S0; I0; R0; CI0; z(1); z(2); z(3)];

    [t, y] = ode45(@(t, y) SIV(t, y), tspan, initialvalues);

    CI = y(:, 4);
    diff = CI - reshape(ConfCase, size(CI));
    value = norm(diff, 2);
end

% Parameter estimation using MultiStart and fmincon
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

% Plot the "best" solution

% Perform simulations using the "best" parameters

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