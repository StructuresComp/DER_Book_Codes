%Alyssa Novelia
%Kirchoff rod under end moment

%Rod parameters
young = 1.3e6;
radius = 0.00079;
poiss = 0.5;
lambda1 = young*pi*(radius^4)/4;
lambda2 = lambda1;
lambda3 = (young/(2*(1+poiss)))*pi*(radius^4)/2;
tfin = 0.2;

%States = [nu1; nu2; nu3; d1x; d1y; d1z; d2x; d2y; d2z; ...
%   d3x; d3y; d3z; x; y; z];
Ydot = @(t,Y)[((lambda2-lambda3)/lambda1) * Y(2) * Y(3); ...
    ((lambda3-lambda1)/lambda2) * Y(3) * Y(1); ...
    ((lambda1-lambda2)/lambda3) * Y(1) * Y(2); ...
    Y(3) * Y(7) - Y(2) * Y(10); ...
    Y(3) * Y(8) - Y(2) * Y(11); ...
    Y(3) * Y(9) - Y(2) * Y(12); ...
    Y(1) * Y(10) - Y(3) * Y(4); ...
    Y(1) * Y(11) - Y(3) * Y(5); ...
    Y(1) * Y(12) - Y(3) * Y(6); ...
    Y(2) * Y(4) - Y(1) * Y(7); ...
    Y(2) * Y(5) - Y(1) * Y(8); ...
    Y(2) * Y(6) - Y(1) * Y(9); ...
    Y(10); ...
    Y(11); ...
    Y(12)];
IC = [-20; -20; 20; 0; 1; 0; 0; 0; 1; 1; 0; 0; 0; 0; 0];

%Solve EOM
[T Y] = ode45(Ydot,[0 tfin], IC);

%Plot result
nu1 = Y(:,1);
nu2 = Y(:,2);
nu3 = Y(:,3);
d1x = Y(:,4);
d1y = Y(:,5);
d1z = Y(:,6);
d2x = Y(:,7);
d2y = Y(:,8);
d2z = Y(:,9);
d3x = Y(:,10);
d3y = Y(:,11);
d3z = Y(:,12);
x = Y(:,13);
y = Y(:,14);
z = Y(:,15);
%%
%Plot centerline of rod
figure; hold on;
refconfig = plot3([0 0],[0 0],[0 tfin],'linewidth',2);
curconfig = plot3(y,z,x,'linewidth',2);
axis equal;
view([-35 10]);
xlabel('y');
ylabel('z');
zlabel('x');

save('kirchoffrodexact.mat', 'x','y','z');
