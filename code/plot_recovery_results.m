function [] = plot_recovery_results(params,x,y,rec_filt,rec_lasso)
%close all
figure
if ~isempty(params)
    param_string = '';
    fields = fieldnames(params);
    for i=1:numel(fields)
        param_string = [param_string '_' (fields{i}) '-'...
            num2str(params.(fields{i}))];
    end
    param_string = strrep(param_string, '.', 'p');
    param_string = strrep(param_string, '/', ':');
    set(gcf,'name', ['recovery_results' param_string]);
else
    set(gcf,'name','recovery_results');
end
subplot(1+(nargin==5),1,1);
T = length(x);
plot(0:T-1,real(x(1:T)),'k-',0:T-1,real(y(1:T)),'b.',...
    0:T-1,real(rec_filt(1:T)),'r-.','Linewidth',2)
legend('Ground Truth', 'Observations', 'Recovery', 'Location','Southwest');
title(['Filtering recovery, ' 'MSE=' num2str(norm(rec_filt(1:T)-x(1:T)))]);
if nargin == 5
    subplot(1+(nargin==5),1,2);
    plot(0:T-1,real(x(1:T)),'k-',0:T-1,real(y(1:T)),'b.',...
        0:T-1,real(rec_lasso(1:T)),'r-.','Linewidth',2)
    legend('Ground Truth', 'Observations', 'Recovery', 'Location','Southwest');
    title(['Lasso recovery, ' 'MSE=' num2str(norm(rec_lasso(1:T)-x(1:T)))]);
end
respath = '../test/rec_test-1d/';
if ~exist(respath, 'dir')
  mkdir(respath);
end
addpath(respath);
set(gcf, 'PaperSize', [11.69 8.27]); % paper size (A4), landscape
% Extend the plot to fill entire paper.
set(gcf, 'PaperPosition', [0 0 11.69 8.27]);
if ~isempty(params)
    savefig(gcf,[respath param_string '.fig']);
    saveas(gcf,[respath param_string '.pdf'],'pdf');
else
    savefig(gcf,[respath 'recovery_results.fig']);
    saveas(gcf,[respath 'recovery_results.pdf'],'pdf');
end
% soundsc(rec)
end