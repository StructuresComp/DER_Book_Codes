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
dataDir = 'datafilesEnergy/';
listing = dir([dataDir,'*txt']);
Nfiles = length(listing);

for c=1:Nfiles
    filename = listing(c).name;
    ind = strfind(filename, '_');
    F = str2double(filename(ind(2)+1:ind(3)-1));
    
    fid = fopen([dataDir, filename], 'r');
    data = textscan(fid, '%f%f%f%f%f%f','CommentStyle','#');
    fclose(fid);
    
    t = data{1};
    freq = data{2};
    Eb = data{3};
    Es = data{4};
    Eg = data{5};
    Ek = data{6};
    
%     ind = find(t==1);
%     Eb0 = Eb(ind(1));
%     Es0 = Eb(ind(1));
%     Eg0 = Eb(ind(1));
%     Ek0 = Eb(ind(1));
    
    h1 = figure(1);
    plot( t, Eb, 'Color', colpos(1,:), 'LineWidth', 0.5);
    hold on
    plot( t, Es, 'Color', colpos(2,:), 'LineWidth', 0.5);
    plot( t, Eg, 'Color', colpos(3,:), 'LineWidth', 0.5);
    plot( t, Ek, 'Color', colpos(4,:), 'LineWidth', 0.5);
    plot( t, Eb + Es + Eg + Ek, 'k', 'LineWidth', 2);
    hold off
    l = legend('E_b', 'E_s', 'E_g', 'E_k', 'Total');
    set(l, 'Fontname', FONT, 'FontSize', FONTSIZE);
    xlim([0 2]);
    xlabel('Time, t [s]','Fontname', FONT,'FontSize',FONTSIZE);
    ylabel('Energy, E [J]','Fontname', FONT,'FontSize',FONTSIZE);
    set(gca, 'Fontname', FONT, 'FontSize', FONTSIZE);
    box on
    
    set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
        'PaperSize', [pWidth pHeight]);
    saveas(h1, num2str(F, 'Fig_Energy_%03d_Hz.pdf'));

%     waitforbuttonpress
end

