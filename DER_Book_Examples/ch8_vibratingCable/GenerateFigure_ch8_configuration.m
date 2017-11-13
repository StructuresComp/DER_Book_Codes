%% Common block
clear all;
close all;
% clc;

FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81; 237 28 36; 0 174 239; 0 0 0; ...
    144 144 144]/255; % colors
ncolors = length(colpos);
syms = '^v<>oshv<>osh';

%%
dataDir = 'datafilesConfiguration/';
listing = dir([dataDir,'*txt']);
Nfiles = length(listing);

L = 0.1;
fprintf('Using length = %f to normalize\n', L);

desiredT = [0; 1; 2];

for c=1:Nfiles
    filename = listing(c).name;
    ind = strfind(filename, '_');
    F = str2double(filename(ind(2)+1:ind(3)-1));
    
    fid = fopen([dataDir, filename], 'r');
    data = textscan(fid, '%f%f%f%f%f','CommentStyle','#');
    fclose(fid);
    
    t = data{1};
    freq = data{2};
    x = data{3};
    y = data{4};
    z = data{5};
    
    for c2 = 1 : length(desiredT)
        ind = find(t == desiredT(c2));
        if (numel(ind) == 0)
            fprintf('No data point. Error!\n');
            continue;
        end
        
        h1 = figure(1);
        plot( x(ind)/L, y(ind)/L, 'k', 'LineWidth', 2);
        axis equal
        xlabel('Norm. x','Fontname', FONT,'FontSize',FONTSIZE);
        ylabel('Norm. y','Fontname', FONT,'FontSize',FONTSIZE);
        set(gca, 'Fontname', FONT, 'FontSize', FONTSIZE);
        box on
        
        set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
            'PaperSize', [pWidth pHeight]);
        imgName = [num2str(F, 'F=%d_Hz_'), num2str(desiredT(c2), 't=%3.1f_second'), '.pdf'];
        saveas(h1, imgName);
    end
end
