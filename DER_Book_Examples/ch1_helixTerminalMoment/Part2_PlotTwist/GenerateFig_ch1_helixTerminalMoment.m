%% Aesthetics
FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; 50 50 50; 25 25 25]/255; % colors
syms = 'o^sv<>p';

%%
plotTime = [1 2 5 10 20]; % plot twist vs. arc length at these times

% Get file name
filename = 'ch1_helixTerminalMomentData.txt';

Length = 0.2; % Length of filament in meter
k1 = 20.0; % moment along m1
k2 = 20.0; % moment along m2
k3 = 20.0; % torsion at one end

% Read the file
fid = fopen(filename);
C = textscan(fid,'%f%f%f%f%f','CommentStyle', '#');
fclose(fid);

t = C{1}; % time
x = C{2}; % x-coordinates
y = C{3}; % y-coordinates
z = C{4}; % z-coordinates
twist = C{5}; % discrete twist

uniqT = unique(t); % times at which snapshots of the rod and its twist are recorded


legendsPlot = plotTime*0;
legendsString = num2str(plotTime', 't=%02d s');

colorNo = 1;
for c2 = 1:length(uniqT)
    
    currentTime = uniqT(c2); % current time at this loop
    currentTime = round(currentTime*100) / 100; % round it out since ...
    % we want to use the == operator next
    ind = find(plotTime == currentTime);
    
    if (numel(ind)==0), continue, end % if not in the list, ignore
    % if the current time is in the list "plotTime", plot!
    
    ind = find(t == uniqT(c2));
    num = length(ind);
    s = (1:num) / num;
    
    edgeLength = Length / num;
    arcLengthParameter = s(2:end-1); % normalized
    DiscreteTwist = twist(ind(2:end-1));
    normDiscreteTwist = DiscreteTwist / edgeLength * Length;
    
    LineWidth = 2.0;
    if (colorNo == length(plotTime)), LineWidth = 1.0; end % manual adjustment
    
    h1 = figure(1);
    hold on
    
    legendsPlot(colorNo) = plot(arcLengthParameter, normDiscreteTwist, ...
        'Color', colpos(colorNo,:), 'LineWidth', LineWidth);
    mn = mean(normDiscreteTwist);
    if (mn == 0), mn = 1; end
    ylim([0 1.2 * mn]);
    hold off
    
    % plot configuration
    h2 = figure(2); hold on;
%    clf;
%     normX = x(ind) / Length;
%     normY = y(ind) / Length;
%     normZ = z(ind) / Length;
    normX = x(ind);
    normY = y(ind);
    normZ = z(ind);
    
    plot3(normX, normY, normZ, '-', ...
        'Color', colpos(colorNo,:), 'LineWidth', LineWidth);
    axis equal
    grid on
    set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
    
    set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
        'PaperSize', [pWidth pHeight]);
    imgName = num2str(currentTime, 'Figure_config_t=%03d.pdf');
    saveas(h2, imgName);
    
    colorNo = colorNo + 1;
end

figure(1);
box on

figure(2);
hold on;
plot3(normX,normY,normZ,'-o');
grid on
xlabel('Norm. arc-length, s', 'Fontname', FONT,'FontSize',FONTSIZE);
ylabel('Discrete twist * total length / edge length', ...
    'Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
l = legend(legendsPlot, legendsString);
set(l,'Fontname', FONT,'FontSize',FONTSIZE);

set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);
saveas(h1, 'Figure_discreteTwist.pdf');
