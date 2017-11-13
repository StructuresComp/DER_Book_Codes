% Khalid Jawed, khalidjm@seas.ucla.edu
% May 2016
% Jacobian verification in DER using a 4-noded rod
%%
clear all;
close all;
clc;

fprintf('Discrete elastic rods: Jacobian Verification for a rod with 4-node\n');

%%
global EI EA GJ 
global nv ne
global refLen voronoiRefLen 
global tangent
global undeformedTwist kappaBar

%% Set up the intial configuration of the four-noded rod
% This section is arbitrarily chosen
node0 = [0; 0; 0];
node1 = [1.0; 0; 0];
node2 = [2.0; 0.2; 0];
node3 = [3.0; 0.2; 0.2];
theta1 = 0;
theta2 = pi/4;
theta3 = pi/4;

%% Arrange the degrees of freedom vector
nv = 4; % number of nodes
ne = nv - 1; % number of edges

% Store into x (size =15 from 4 nodes and 3 edges)
x = zeros(3*nv + ne, 1);
x(1:3) = node0;
x(5:7) = node1;
x(9:11) = node2;
x(13:15) = node3;
x(4) = theta1;
x(8) = theta2;
x(12) = theta3;

% Compute undeformed twist
theta_f = x(4:4:end);
theta_e = [0; theta_f(1:end-1)];
undeformedTwist = theta_f - theta_e;

% Elastic stiffness
% We set all to unity for simplicity
EI = 1;
EA = 1;
GJ = 1;

%% Reference length and Voronoi length
refLen = zeros(ne, 1);
for c=1:ne
    dx = x(4*c+1:4*c+3) - x(4*(c-1)+1:4*(c-1)+3);
    refLen(c) = norm(dx);
end
voronoiRefLen = zeros(nv, 1);
for c=1:nv
    if c==1
        voronoiRefLen(c) = 0.5 * refLen(c);
    elseif c==nv
        voronoiRefLen(c) = 0.5 * refLen(c-1);
    else
        voronoiRefLen(c) = 0.5 * (refLen(c-1) + refLen(c));
    end
end

%% Reference director & material director
d1 = zeros(ne, 3); % reference director, u (or d1)
d2 = zeros(ne, 3); % reference director, v (or d2)
tangent = zeros(ne, 3); % tangent
for c=1:ne
    dx = x(4*c+1:4*c+3) - x(4*(c-1)+1:4*(c-1)+3);
    tangent(c,:) = dx / norm(dx);
end

% Figure out a good choice for d1(1)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0, t1);
if (abs(d1Tmp) < 1.0e-6)
    t1 = [0 1 0];
    d1Tmp = cross(t0, t1);
end
d1(1,:) = d1Tmp;

d1_l = d1(1,:);
d2(1,:) = cross(t0, d1_l);
for c=2:ne
    t1 = tangent(c,:);
    d1_l = parallel_transport(d1_l, t0, t1);
    d1_l = (d1_l - dot(d1_l, t1) * t1);
    d1_l = d1_l / norm(d1_l);
    d1(c,:) = d1_l;
    d2_l = cross(t1, d1_l);
    d2(c,:) = d2_l;
    t0 = t1;
end

% Compute material directors
theta0 = x(4:4:end); % Extract only the theta components
[m1, m2] = computeMaterialDirectors(d1, d2, theta0); % compute m1 and m2
refTwist = zeros(ne,1);
% refTwist = getRefTwist(d1, tangent, refTwist); % Should be all zero

% Natural curvature
kappaBar = getkappaBar( x, m1, m2 );

%% Verify Jacobian
% Compute forces and 'analytical' Jacobians for the given configuration
[Fbend, JAnalyticalBend] = getFb(x, m1, m2);
[Ftwist, JAnalyticalTwist] = getFt(x, refTwist);
[Fstretch, JAnalyticalStretch] = getFs(x);

% Compute forces and numerically evaluate Jacobian
% Create containers for Jacobians first)
JNumericalBend = JAnalyticalBend*0;
JNumericalTwist = JAnalyticalTwist*0;
JNumericalStretch = JAnalyticalStretch*0;
for k1=1: length(x)
    xPerturbed = x;
    xPerturbed(k1) = xPerturbed(k1) + 1.0e-4; % add some perturbation
    
    % Compute tangents
    tangentIterate = computeTangent(xPerturbed);
    % Compute reference directors
    [d1Iterate, d2Iterate] = computeTimeParallel(d1, x, xPerturbed);
    % Compute reference twist
    refTwistIterate = getRefTwist(d1Iterate, tangentIterate, refTwist);
    thetaIterate = x(4:4:end);
    % Compute material directors
    [m1Iterate, m2Iterate] = computeMaterialDirectors(d1Iterate, d2Iterate, thetaIterate);

    % Compute elastic forces for perturbed configuration
    [Fbend_Perturbed, ~] = getFb(xPerturbed, m1Iterate, m2Iterate);
    [Ftwist_Perturbed, ~] = getFt(xPerturbed, refTwistIterate);
    [Fstretch_Perturbed, ~] = getFs(xPerturbed);
    
    dx = xPerturbed(k1) - x(k1);
    % Fill in 'numerical' Jacobian (bending)
    dF = Fbend_Perturbed - Fbend;
    JNumericalBend(k1, :) = dF / dx;
    % Fill in 'numerical' Jacobian (twisting)
    dF = Ftwist_Perturbed - Ftwist;
    JNumericalTwist(k1, :) = dF / dx;
    % Fill in 'numerical' Jacobian (stretching)
    dF = Fstretch_Perturbed - Fstretch;
    JNumericalStretch(k1, :) = dF / dx;
end

%% Multiply by -1 to get Jacobian of the energy
JAnalyticalBend = - JAnalyticalBend;
JNumericalBend = - JNumericalBend;
JAnalyticalTwist = - JAnalyticalTwist;
JNumericalTwist = - JNumericalTwist;
JAnalyticalStretch = - JAnalyticalStretch;
JNumericalStretch = - JNumericalStretch;

%% Plot
FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81; 237 28 36; 0 174 239; 247 230 130; ...
    247 208 130]/255; % colors

[n, ~] = size(JAnalyticalBend);
nn = n^2 - 1;

%h1 = figure(2);
clf();
%subplot(3,1,1);
figure(100);
hold on
for i = 1:n
    for j = 1:n
%         if (abs(i - j) <= 5)
            component = n * (i-1) + (j-1);
            plot(component, JAnalyticalBend(i, j), ...
                'sq', 'Color', colpos(2,:), 'MarkerFaceColor', colpos(2,:));
            plot(component, JNumericalBend(i, j), ...
                'o', 'Color', colpos(3,:), 'MarkerSize', 7);
%         end
    end
end
hold off
%box on
%xticks([0 100 nn]);
xlim([0 nn]);
%l = legend('Numerical', 'Analytical', 'location', 'NorthEast');
%ylabel('H_b');
%xlabel('(4N-1) * i + j');
%set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
%set(l,'Fontname', FONT,'FontSize',FONTSIZE/2);

figure(200);
%subplot(3,1,2);
hold on
for i = 1:n
    for j = 1:n
%         if (abs(i - j) <= 5)
            component = n * (i-1) + (j-1);
            plot(component, JAnalyticalTwist(i, j), ...
                'sq', 'Color', colpos(2,:), 'MarkerFaceColor', colpos(2,:));
            plot(component, JNumericalTwist(i, j), ...
                'o', 'Color', colpos(3,:), 'MarkerSize', 7);
%         end
    end
end
hold off
%box on
%xticks([0 100 nn]);
xlim([0 nn]);
%l = legend('Numerical', 'Analytical', 'location', 'NorthEast');
%ylabel('H_t');
%xlabel('(4N-1) * i + j');
%set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
%set(l,'Fontname', FONT,'FontSize',FONTSIZE/2);

%subplot(3,1,3);
figure(300);
hold on
for i = 1:n
    for j = 1:n
%         if (abs(i - j) <= 5)
            component = n * (i-1) + (j-1);
            plot(component, JAnalyticalStretch(i, j), ...
                'sq', 'Color', colpos(2,:), 'MarkerFaceColor', colpos(2,:));
            plot(component, JNumericalStretch(i, j), ...
                'o', 'Color', colpos(3,:), 'MarkerSize', 7);
%         end
    end
end
hold off
%box on
%xticks([0 100 nn]);
xlim([0 nn]);
%l = legend('Numerical', 'Analytical', 'location', 'NorthEast');
%ylabel('H_s');
%xlabel('(4N-1) * i + j');
%set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
%set(l,'Fontname', FONT,'FontSize',FONTSIZE/2);

%set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth 3*pHeight], ...
%    'PaperSize', [pWidth 3*pHeight]);
%saveas(h1, 'JacobianVerification.pdf');
