clear all; %#ok<*CLALL>
close all;

%% Aesthetics
FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239]/255; % colors

%%
load helixParams

%% Transform coords
% Exact solution
n = [A_exact(4); A_exact(5); sqrt(1 - (A_exact(4))^2 - (A_exact(5))^2)];
O = A_exact(1:3); % origin
ptTransformedExact = getTransformed(O, n, x, y, z);
% Ensure the helix starts at z=0
ptTransformedExact(:,3) = ptTransformedExact(:,3) - ptTransformedExact(1,3);
% Get total length for normalization
s_Exact = cumsum(sqrt(sum(diff(ptTransformedExact,[],1).^2,2)));
L_exact = s_Exact(end); % length
ptTransformedExact = ptTransformedExact/ L_exact;

% DER
n = [A_DER(4); A_DER(5); sqrt(1 - (A_DER(4))^2 - (A_DER(5))^2)];
n = -1*n; % This is adhoc to match exact and DER solution. This is OK. We are just flipping the helix axis
O = A_DER(1:3); % origin
ptTransformedDER = getTransformed(O, n, normX, normY, normZ);
ptTransformedDER(:,3) = ptTransformedDER(:,3) - ptTransformedDER(1,3);
% Get total length for normalization
s_DER = cumsum(sqrt(sum(diff(ptTransformedDER,[],1).^2,2)));
LDER = s_DER(end); % length
ptTransformedDER = ptTransformedDER/ LDER;

h1 = figure(1);
   
plot3(ptTransformedDER(:,1), ptTransformedDER(:,2), ptTransformedDER(:,3), ...
    'o', 'Color', colpos(1,:));
hold on
plot3(ptTransformedExact(:,1), ptTransformedExact(:,2), ptTransformedExact(:,3), ...
    'k--', 'LineWidth', 2);
hold off

box on
l = legend('DER', 'Exact', 'Location', 'NorthEast');
set(l, 'Fontname', FONT,'FontSize',FONTSIZE);
xlabel('Norm. x-axis','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('Norm. y-axis','Fontname', FONT,'FontSize',FONTSIZE);
zlabel('Norm. z-axis','Fontname', FONT,'FontSize',FONTSIZE);

set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, 'Fig_Helix_Comparison.pdf');
