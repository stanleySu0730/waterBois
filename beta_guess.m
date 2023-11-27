%%
close;
clear all;

rng('default') 

%% How many runs for parameter estimation do you want:
NoStartPoints=10;

%% http://publichealth.lacounty.gov/media/coronavirus/
% LA County Data - April 20, 2020 (START DATE)  through May 16, 2020
ConfCase = [13823 15165	16449 17567	18545 19159	19567 20460	21017 22522	23233 24262	24936 25699	26238 27866	28665 29526	30334 31241	31703 32269	33247 34552	35447 36324	37374];
Deaths = [619 666 732 798 850 896 916 948 1004 1065	1119 1174 1212 1231	1260 1315 1369 1420	1470 1515 1531 1570 1617 1663 1711 1752	1793];
ten_dayCumulative = 5370;
CumulativeTo10DayStart = 8453;

% Total Population considered for outbreak
TotalPopulation = 9651332; 

day = length(ConfCase);        
tspan = 1:1:day;               

I0 = ten_dayCumulative;                %% Add together all the cases 10 days prior to the start date
E0 = .40*CumulativeTo10DayStart;       %% Fit to data
R0 = .95*20*CumulativeTo10DayStart;    %% Changed to a percentage of cumulative infected up to 10 days before start date
CI0 = ConfCase(1);                     %% total infections

S0 = TotalPopulation-I0-R0-E0;

initialvalues = [S0,I0,E0,R0,CI0];

%% Upper and Lower bounds for parameters you are fitting.
%Parameter vector we are approximating:
%z = [beta] 

LowerBounds=[1.00e-9];      %Lowerbounds for the parameters you are estimating
UpperBounds=[1.00e-07];     %Upperbounds for the parameters you are estimating

xstart=.5*(LowerBounds+UpperBounds);                            

%% MultiStart and fmincon 
problem = createOptimProblem('fmincon','objective', ...
                              @SIR_RUN_ODE45,'x0',xstart, ...
                              'lb',LowerBounds,'ub',UpperBounds);

problem.options = optimoptions(problem.options, ...
                               'MaxFunEvals',9999, ...
                               'MaxIter',9999);
numstartpoints = NoStartPoints;                              
ms = MultiStart('UseParallel',true,'Display','iter');         
[b,fval,exitflag,output,manymins] = run(ms,problem,numstartpoints);  

for i=1:length(manymins)
    SIRParameters(i,:)=manymins(i).X;       
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;           
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;     
end

%% Plot
[t,y] = ode45(@(t,y) SEIR_Model(t,y,SIRParameters(1,:)),tspan,initialvalues);

S = y(:,1);
I = y(:,2);
R = y(:,3);
E = y(:,4);
CI = y(:,5);

figure(1)
tiledlayout(1,5)
nexttile
hold all
plot(tspan,CI,'LineWidth',2,'LineStyle',':')
scatter(tspan,ConfCase,'filled')
title('Cumulative Cases')
xlabel('Days')
nexttile
plot(tspan,E,'LineWidth',2,'LineStyle',':')
title('Exposed')
xlabel('Days')
nexttile
plot(tspan,S,'LineWidth',2,'LineStyle',':')
title('Susceptible')
xlabel('Days')
nexttile
plot(tspan,I,'LineWidth',2,'LineStyle',':')
title('Infectious')
xlabel('Days')
nexttile
plot(tspan,R,'LineWidth',2,'LineStyle',':')
title('Recovered')
xlabel('Days')

%% Run the ode45 solvers  
function value = SIR_RUN_ODE45(z) 
ConfCase = [13823 15165	16449 17567	18545 19159	19567 20460	21017 22522	23233 24262	24936 25699	26238 27866	28665 29526	30334 31241	31703 32269	33247 34552	35447 36324	37374];
Deaths = [619 666 732 798 850 896 916 948 1004 1065	1119 1174 1212 1231	1260 1315 1369 1420	1470 1515 1531 1570 1617 1663 1711 1752	1793];
ten_dayCumulative = 5370;
CumulativeTo10DayStart = 8453;

% Total Population considered for outbreak
TotalPopulation = 9721138; 

day = length(ConfCase);        
tspan = 1:1:day;               

I0 = ten_dayCumulative;                   %% Add together all the cases 10 days prior to the start date
E0 = .40*20*ten_dayCumulative;                  %% Think about exposed 
R0 = .95*20*CumulativeTo10DayStart;    %% Changed to a percentage of cumulative infected up to 10 days before start date
CI0 = ConfCase(1);                     %% Total infections

S0 = TotalPopulation-I0-R0-E0;

initialvalues = [S0,I0,E0,R0,CI0];

[t,y] = ode45(@(t,y) SEIR_Model(t,y,z),tspan,initialvalues);        
CI = y(:,4);
       
diff = CI - reshape(ConfCase,size(CI));

value = norm(diff,2);
end 

%%
function dydt = SEIR_Model(t,y,z)
beta = z(1);        
kappa = 1/2;
eta = 1/10; 
gamma = 1/6;
sigma = 1/4;
epsilon = 0.02; 

% Initiate DE variables
dydt = zeros(5,1);

S = y(1);
E = y(2);
I = y(3);
R = y(4);

dydt(1) = -beta*S*I;
dydt(2) = beta*S*I-kappa*E-eta*E;
dydt(3) = kappa*E-gamma*(1-epsilon)*I;
dydt(4) = gamma*(1-epsilon)*I+eta*E-sigma*R;
dydt(5) = beta*S*I;

end

