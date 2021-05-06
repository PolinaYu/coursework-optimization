clear;
% 
%          # social distancing (0-1)
%          # 0   = no social distancing
%          # 0.1 = masks
%          # 0.2 = masks and hybrid classes
%          # 0.3 = masks, hybrid, and online classes
% Мы будет рассматривать начало развития эпидемии в Германии, поэтому
% примем u=0
% SEIR model


%% Set Parameters
u =  0;
I0 = 1e-5; % Initial proportion of infected.
beta = 3;  % S to E coefficient( wk^-1)
alfa = 0.5; %E to I coeff( wk^-1)
gamma = 0.5; % I to R ( wk^-1)
tmax = 104; % Number of weeks
dt = .01;  %size of time steps in weeks
Imax = 1.1;% max proportion of infecred(is actually 1, but we need a buffer)


%%Initialize Vectors
t = 0:dt:tmax; % time vector
Nt = length(t); % number of time steps

S = zeros(1, Nt); % Susceotible vector
E = zeros(1, Nt); % Exposed vector
I = zeros(1, Nt); % Infection vector
R = zeros(1, Nt); % Removed vector 
I(1) = I0; % set initial infection value


%%Calculations
for it = 1:Nt-1
   S(it) = 1 - E(it) - I(it)- R(it);
   dE = (1-u)*beta*I(it)*S(it)-alfa*E(it); % rate of change of exposed per week
   E(it+1) =E(it)+ dE*dt;
   dI = alfa*E(it)- gamma*I(it); % rate of change per week
   I(it+1) =I(it)+ dI*dt;
   dR =gamma*I(it);
   R(it+1) = R(it)+dR*dt;
end
S(Nt)= 1-E(Nt)-I(Nt)-R(Nt);


%%Plots
plotchoice = 5; % 1 =S, 2 = E, 3= I, 4 = R, 5 = all
switch plotchoice
    case 1
        plot(t,S, '-b', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time (weeks)')
        ylabel('Proportion Susceptible')
        title('Proportion of Susceptible vs. Time')
    case 2
        plot(t,E, '-g', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time (weeks)')
        ylabel('Proportion Exposed')
        title('Proportion of Exposed vs. Time')
    case 3
        plot(t,I, '-r', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time (weeks)')
        ylabel('Proportion Infected')
        title('Proportion of infections vs. Time')
    case 4
        plot(t, R, '-m', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time (weeks)')
        ylabel('Proportion Removed')
        title('Proportion of Removed vs. Time')
    case 5
         plot(t, S, '-b',...
              t, E, '-g',...
              t, I, '-r',...
              t, R, '-m', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time (weeks)')
        ylabel('Proportions of S,E, I and R')
        title('Proportion of S, E, I and R vs. Time')
        legend('S', 'E', 'I', 'R',...
            'Location', 'SouthEastOutside')
        
end



