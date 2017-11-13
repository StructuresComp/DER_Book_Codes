clear all; %#ok<*CLALL>

%% Aesthetics
FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 50 50 50; 25 25 25]/255; % colors
syms = 'o^sv<>p';

%%
global currentX currentY currentZ M1 N1

%% Read the data
% Get file name
filename = 'ch1_helixTerminalMomentData.txt';

% Read the file
fid = fopen(filename);
C = textscan(fid,'%f%f%f%f%f','CommentStyle', '#');
fclose(fid);

t = C{1}; % time
x = C{2}; % x-coordinates
y = C{3}; % y-coordinates
z = C{4}; % z-coordinates
twist = C{5}; % discrete twist

ind = find(t == max(t));

normX = x(ind);
normY = y(ind);
normZ = z(ind);

% Exact solution
load kirchoffrodexact

%%
options = optimset('TolFun', 1.0e-12, 'TolX', 1.0e-12, ...
    'MaxIter', 1.0e5, 'MaxFunEvals', 1.0e6);

%% Fit using exact data
currentX = x; %#ok<*NASGU>
currentY = y;
currentZ = z;

% Assume that the axis of the helix is ax + by + c = 0
origin1 = [-0.012; 0.0289; -0.0056];
dirN = [-0.6396; -0.6396];
for k=1:2
    M1 = origin1;
    dirN = fminsearch(@helixDistanceDirection, dirN, options);
    A_exact = [M1; dirN];
    
    N1 = dirN;
    origin1 = fminsearch(@helixDistanceOrigin, origin1, options);
    A_exact = [origin1; dirN];
    
    fprintf('residual (exact)=%e\n', helixDistance(A_exact));
end

% % Uncomment the following figure for debugging
% figure(1);
% plot3(x, y, z, 'k-');
% hold on
% n = [A_exact(4) A_exact(5) sqrt(1 - A_exact(4)^2 - A_exact(5)^2)]';
% inPoint = [A_exact(1) A_exact(2) A_exact(3)]';
% pt1 = inPoint - n * 0.1;
% pt2 = inPoint + n * 0.1;
% plot3([pt1(1) pt2(1)], [pt1(2) pt2(2)], ...
%     [pt1(3) pt2(3)], 'ro-');
% hold off

%% Fit using DER data
currentX = normX;
currentY = normY;
currentZ = normZ;

% Assume that the axis of the helix is ax + by + c = 0
origin1 = [0.0506; 0.0025; 0.0096];
dirN = [-0.9070; -0.1715];
for k=1:2
    M1 = origin1;
    dirN = fminsearch(@helixDistanceDirection, dirN, options);
    A_DER = [M1; dirN];
    
    N1 = dirN;
    origin1 = fminsearch(@helixDistanceOrigin, origin1, options);
    A_DER = [origin1; dirN];
    
    fprintf('residual=%e\n', helixDistance(A_DER));
end

% % Uncomment the following figure for debugging
% figure(2);
% plot3(normX, normY, normZ, 'k-');
% hold on
% n = [A_DER(4) A_DER(5) sqrt(1 - A_DER(4)^2 - A_DER(5)^2)]';
% inPoint = [A_DER(1) A_DER(2) A_DER(3)]';
% pt1 = inPoint - n * 0.1;
% pt2 = inPoint + n * 0.1;
% plot3([pt1(1) pt2(1)], [pt1(2) pt2(2)], ...
%     [pt1(3) pt2(3)], 'ro-');
% hold off

%%
save('helixParams.mat', 'normX', 'normY', 'normZ', 'x', 'y', 'z', ...
    'A_exact', 'A_DER');
