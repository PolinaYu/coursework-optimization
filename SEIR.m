function SEIR
%Math-Works website is also useful for this. Wikipedia explains the 4 odes
%nicely.
u =  0;
I0 = 1e-5; % Initial proportion of infected.
beta = 3;  % S to E coefficient( wk^-1)
alfa = 0.5; %E to I coeff( wk^-1)
gamma = 0.5; % I to R ( wk^-1)
tmax = 104; % Number of weeks
dt = .01;  %size of time steps in weeks
t = 0:dt:tmax; % time vector
Nt = length(t); % number of time steps

    function dF = system(t, x)     %set up the function to be solved in a seperate function environmen
        x(1) = 1-x(2)-x(3)-x(4);
        ds = -(1 - u)*x(1)*x(3)*beta;  %set of 4 odes, I just found them on wikipedia to then use here.
        de = (1 - u)*beta*x(3)*x(1) - alfa*x(2);
        di = alfa*x(2) - gamma *x(3);
        dr = gamma*x(3);
        dF =[ds; de; di; dr];
    end
options = odeset('Refine', 10, 'RelTol', 1e-4);     %allows us to tailor ode45's options
[t,y] = ode45(@system, t, [I0 0 0 0], options); %call ode45 to solve the model "rigid".

plot(t, y)
title('SEIR model')
legend('S','E','I', 'R')
end
