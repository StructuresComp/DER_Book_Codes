clear all; %#ok<*CLALL>
clc;

%% Aesthetics
FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 50 50 50; 25 25 25]/255; % colors
syms = 'o^sv<>p';

%%
fid = fopen('ch4_helixUncoilingData.txt', 'r');
C = textscan(fid,'%f%f%f%f%f%f','CommentStyle', '#');
fclose(fid);

t = C{1}; % time
x = C{2}; % x-coordinates
y = C{3}; % y-coordinates
z = C{4}; % z-coordinates
theta = C{5}; % discrete twist
refTwist = C{6}; % ref twist

uniqT = unique(t); % times at which snapshots of the rod and its twist are recorded

colorNo = 1;

h1 = figure(1);
clf();
hold on

h2 = figure(2);
clf();
hold on

for c2 = [1 21 51 1001]
    
    ind = find(t == uniqT(c2));
    num = length(ind);
    s = (1:num) / num;
    ne = num - 1;
    
    theta_i = theta( ind(1:ne-1) );
    theta_f = theta( ind(2:  ne) );    
    refTwistLocal = refTwist(ind(1:ne));
    twistLocal = theta_f - theta_i + refTwistLocal(2:end);
    twistLocal = [0; twistLocal]; %#ok<*AGROW> % twist=0 at edge=0
    xx = x(ind);
    yy = y(ind);
    zz = z(ind);
    
    
    figure(1);
    hold on
    subplot(3,4,[1:3, 5:7, 9:11]);
    lh = plot(1:ne, refTwistLocal, 'linewidth',2,'Color', colpos(colorNo,:));
    
    figure(1);
    hold on
    subplot(3,4,[4 8 12]);
    plot3(xx,yy,zz,'.-','Color', colpos(colorNo,:)); 
    
    figure(2);
    hold on
    plot(1:ne, twistLocal, 'linewidth',2,'Color', colpos(colorNo,:));
    colorNo = colorNo + 1;
end

figure(1);
subplot(3,4,[1:3, 5:7, 9:11]);
hold off;
box on
axis([1 ne -2 1.5]);
legend(['t = ' num2str(uniqT(1)) ' s'], ...
    ['t = ' num2str(uniqT(21)) ' s'], ...
    ['t = ' num2str(uniqT(51)) ' s'], ...
    ['t = ' num2str(uniqT(1001)) ' s'], ...
    'Location', 'SouthWest');
xlabel('edge number','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('reference twist','Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);


subplot(3,4,[4 8 12]);
hold off;
box on
axis equal;
xlabel('x','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('y','Fontname', FONT,'FontSize',FONTSIZE);
zlabel('z','Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);

set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 2*pWidth 2*pHeight], ...
    'PaperSize', [2*pWidth 2*pHeight]);
saveas(h1, 'Figure_refTwist.pdf');

figure(2);
hold off
box on
legend(['t = ' num2str(uniqT(1)) ' s'], ...
    ['t = ' num2str(uniqT(21)) ' s'], ...
    ['t = ' num2str(uniqT(51)) ' s'], ...
    ['t = ' num2str(uniqT(1001)) ' s'], ...
    'Location', 'NorthEast');
xlim([1 ne]);
xlabel('edge number','Fontname', FONT,'FontSize',FONTSIZE);
ylabel('twist','Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, 'Figure_Twist.pdf');
