function [] = save_plot2(mode,save_pics,x,what)
%close all
figure
if strcmp(mode,'img'), imshow(uint8(x)); end
if strcmp(mode,'mesh'), 
    mesh(x); zinf = max(max(abs(x))); zlim([-zinf, zinf]); 
end
if strcmp(mode,'surf'),
    surf(x); zinf = max(max(abs(x))); zlim([-zinf, zinf]); shading interp;
end
if strcmp(mode,'gamma'),
    pcolor(x); zinf = max(max(abs(x))); zlim([-zinf, zinf]); shading flat;
    if strcmp(mode,'gamma'), view([0 0 1]); axis off; end
end
if strcmp(mode,'grey'), imshow(x, []); end
if strcmp(mode,'mine'), csf = 50; image(csf*real(x)); axis off; end
respath = save_pics;
if ~exist(respath, 'dir')
  mkdir(respath);
end
addpath(respath);
fig = gcf;
savefig(fig,[respath what '.fig']);
set(gcf, 'PaperSize', [11.69 8.27]); % paper size (A4), landscape
% Extend the plot to fill entire paper.
set(gcf, 'PaperPosition', [0 0 11.69 8.27]);
saveas(gcf,[respath what '.pdf'],'pdf');
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '-depsc', [respath what '.eps']);
close(fig)
end