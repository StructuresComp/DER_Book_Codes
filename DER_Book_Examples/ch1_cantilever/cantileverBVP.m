%Numerical solution to the cantilevered rod problem
%Solving the boundary value problem using bvp4c
%Alyssa Novelia 9/7/2017

% Rod length (m)
RodLength = 0.1;

% Volumetric density (kg m^-3)
rho = 1000;

% Cross-sectional radius of rod (m)
r0 = 1e-3;

% Length density (kg m^-1)
rho0 = rho * pi * r0^2;

% Young's modulus (Pa)
Y = 1e6;

% Gravity (m s^-2)
g = 9.81;

% Flexural rigidity (N m^-2)
EI = Y * pi * r0^4 / 4;

% Nondimensional constant
alpha = rho0 * g * RodLength^3 / EI;

% Equation of motion
odefun = @(x,y) [y(2); alpha*(1-x)*cos(y(1))];

% Boundary condition
hangingbc = @(ya,yb) [ya(1); yb(2)];

% Guess solution
solinit = bvpinit(linspace(0,1,10),[-1 0]);

% Calculate solution
sol = bvp4c(odefun,hangingbc,solinit);

% Solution interval
x = linspace(0,1);

%Plot solution
y = deval(sol,x);

xi = x*RodLength;
dxi = diff(xi);
th = y(1,:);

X = zeros(1,numel(th));
Y = X;

for iii = 2:numel(X)
    X(iii) = X(iii-1) + dxi(iii-1)*cos(th(iii));
    Y(iii) = Y(iii-1) + dxi(iii-1)*sin(th(iii));
end

plot(X,Y);

save('hangingrodtransient.mat', 'X','Y');

