function [] = plot_recovery_results2(params,mode,save_pics,...
                                     x,y,rec_filt,rec_lasso)
if save_pics, % save separate figures
    save_plot2(mode,save_pics,x,'x');
    save_plot2(mode,save_pics,y,'y');
    save_plot2(mode,save_pics,rec_filt,'rec_filter');
    if nargin==7, 
        save_plot2(mode,save_pics,rec_lasso,'rec_lasso');
    end
end
close all
figure
if ~isempty(params)
    param_string = '';
    fields = fieldnames(params);
    for i=1:numel(fields)
        param_string = [param_string (fields{i}) '-' ...
            num2str(params.(fields{i})) '_'];
    end
    param_string = strrep(param_string, '.', 'p');
    param_string = strrep(param_string, '/', ':');
    set(gcf,'name', ['recovery_results_' param_string]);
else
    set(gcf,'name','recovery_results');
end
if strcmp(mode,'mesh') || strcmp(mode,'surf') || strcmp(mode,'gamma'),
    zinf = max([max(max(abs(x))) max(max(abs(y)))...
            max(max(abs(rec_filt))) max(max(abs(rec_lasso)))]);
end
if strcmp(mode,'mine'), csf = 50; end
subplot(2,2,1);
if strcmp(mode,'img'), imshow(uint8(x)); end
if strcmp(mode,'mesh'), mesh(x); zlim([-zinf, zinf]); end
if strcmp(mode,'surf') || strcmp(mode,'gamma'),
    surf(x); zlim([-zinf, zinf]); shading interp; 
    if strcmp(mode,'gamma'), view([0 0 1]); axis off; end
end
if strcmp(mode,'grey'), imshow(x, []); end
if strcmp(mode,'mine'), image(csf*x); axis off; end
title('Ground truth');
subplot(2,2,2);
if strcmp(mode,'img'), imshow(uint8(y)); end
if strcmp(mode,'mesh'), mesh(y); zlim([-zinf, zinf]); end
if strcmp(mode,'surf') || strcmp(mode,'gamma'),
    surf(y); zlim([-zinf, zinf]); shading interp;
    if strcmp(mode,'gamma'), view([0 0 1]); axis off; end
end
if strcmp(mode,'grey'), imshow(y, []); end
if strcmp(mode,'mine'), image(csf*y); axis off; end
title(['Observations, ' 'MSE=' num2str(norm(y(:)-x(:)))]);
subplot(2,2,3);
if strcmp(mode,'img'), imshow(uint8(rec_filt)); end
if strcmp(mode,'mesh'), mesh(rec_filt); zlim([-zinf, zinf]); end
if strcmp(mode,'surf') || strcmp(mode,'gamma'),
    surf(rec_filt); zlim([-zinf, zinf]); shading interp;
    if strcmp(mode,'gamma'), view([0 0 1]); axis off; end
end
if strcmp(mode,'grey'), imshow(rec_filt, []); end
if strcmp(mode,'mine'), image(csf*rec_filt); axis off; end
title(['Filtering recovery, ' 'MSE=' num2str(norm(rec_filt(:)-x(:)))]);
if nargin == 7,
    subplot(2,2,4);
    if strcmp(mode,'img'), imshow(uint8(rec_lasso)); end
    if strcmp(mode,'mesh'), mesh(rec_lasso); zlim([-zinf, zinf]); end
    if strcmp(mode,'surf') || strcmp(mode,'gamma'),
        surf(rec_lasso); zlim([-zinf, zinf]); shading interp;
        if strcmp(mode,'gamma'), view([0 0 1]); axis off; end
    end
    if strcmp(mode,'grey'), imshow(rec_lasso, []); end
    if strcmp(mode,'mine'), image(csf*rec_lasso); axis off; end
    title(['Lasso recovery, ' 'MSE=' num2str(norm(rec_lasso(:)-x(:)))]);
end
respath = '../test/rec_test-2d/';
if ~exist(respath, 'dir')
  mkdir(respath);
end
addpath(respath);
savefig(gcf,[respath param_string '.fig']);
set(gcf, 'PaperSize', [11.69 8.27]); % paper size (A4), landscape
% Extend the plot to fill entire paper.
set(gcf, 'PaperPosition', [0 0 11.69 8.27]);
saveas(gcf,[respath param_string '.pdf'],'pdf');
end