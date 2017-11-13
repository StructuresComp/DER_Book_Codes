FONT = 'Arial';
FONTSIZE = 12;
pWidth = 4; % inches
pHeight = 3;
colpos = [247 148 30;0 166 81;237 28 36;0 174 239; ...
    0 50 150; 150 50 0]/255; % colors
syms = 'o^v<>p';

dirName = 'datafiles/';
Files = dir([dirName, '*.txt']);
N = length(Files);
nv = zeros(N,1);
tipX = nv*0;
tipY =  nv*0;

load hangingrodtransient
bvp4c_normX = X/0.1;
bvp4c_normY = Y/0.1;

for c=1:N
    fname = Files(c).name;
    ind = strfind(fname, '_');
    
    len = str2double(fname(ind(2)+1: ind(3)-1));
    nv(c) = str2double(fname(ind(4)+1: end-4));
    
    fid = fopen([dirName, fname], 'r');
    data = textscan(fid, '%f%f%f%f','CommentStyle','#');
    
    t = data{1};
    x = data{2} / len;
    y = data{3} / len;
    z = data{4} / len;
    maxT = max(t);
    ind = find(t == maxT);
    
    h1 = figure(1);
    clf();
    hold on
    
    xDER = x(ind);
    yDER = y(ind);
    
    p1 = plot(xDER, yDER, 'k-', 'LineWidth', 2); % Draw the rod centerline
    plot(xDER(3:end),yDER(3:end), 'o', 'Color', 'none', 'MarkerFaceColor', ...
        'k', 'LineWidth', 2);  % draw the free nodes
    plot(xDER(1:2),yDER(1:2), 'o', 'Color', 'none', 'MarkerFaceColor', ...
        'r', 'LineWidth', 2); % draw the fixed nodes

    tipX(c) = x(ind(end));
    tipY(c) = y(ind(end));
    
    p2 = plot(bvp4c_normX, bvp4c_normY, '--', 'Color', colpos(2,:), ...
        'LineWidth', 2); % bvp4c solution
    hold off
    title(num2str(nv(c), 'Number of vertices=%d'));
    ylim([-len 0]);
    axis equal
    box on
    set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
    l = legend([p1, p2], 'DER', 'bvp4c');
    set(l,'Fontname', FONT,'FontSize',FONTSIZE);
    
    set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
        'PaperSize', [pWidth pHeight]);
    imgName = num2str(nv(c), 'Convergence_C++_nv=%d.pdf');
    saveas(h1, imgName);
end
