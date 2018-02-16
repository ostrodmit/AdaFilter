function [x,y,rec_filt,rec_lasso] = recovery_test2(params,input_control)
%
% RECOVERY_TEST2 Test of pointwise denoising for a 2-D signal.
%
% Input parameters
%   params : structure with test parameters, with the following fields:
%       *sce [{'<filename>','RanSin-<k>'}] : testing example. 
%           If 'RanSin-<k>, the artificial signal (random 2-d sines) will
%           be generated, where <k> is the number of sines. Otherwise,
%           image file '<filename>' will be used.
%       snr [double,1] : signal-to-noise ratio (*not* in dBs).
%       lep [{0,1},1] : if 1, use the Lepski-like adaptation procedure.
%       mode [{0-2},1] : filtering windows mode. 0 : pointwise mode,
%           the filter is set up separately for each point. 1 : blockwise
%           mode, signal is divided into blocks, each of the size twice the
%           filter bandwidth, and the blocks are processed separately. 2 :
%           same but the blocks overlap by half.
%       par [{0,1},0] : use parfor. Compatible with neighb = 1 (see below).
%       verb [{0,1},1] : show the progress.
%       neighb [{0,1},1] : boost filter fitting by initializing by solution
%           for the previously processed sample. Compatible with par=1.
%       w [{int > 0},50] : filter bandwidth (not relevant if lep=1).
%       n [{int > 0},132] : signal length. Must be at least 4*w if mode=0,
%           2*w if mode=1, and 3*w if mode=2.
%   input_control : structure of control parameters for FIT_FILTER. The
%       user may precise just the desired fields, the others are inialized
%       by default (see FIT_FILTER), with rho = 4 and eps = sigma*rho/2
%       where sigma is the noise level defined by <snr>.
%
% Output parameters
%   x : true signal.
%   y : observations.
%   rec_filt : recovery by filtering.
%   rec_lasso : recovery by the Lasso (baseline).
%
% Example : Brodatz pic, Bandwidth adaptation, blockwise :
%   params = struct('sce','../Brodatz/D75.gif','snr',1.5,'lep',1,'mode',1);
%   [x,y,rec_filt,rec_lasso] = recovery_test2(params);
%
% See also : recovery_test, fit_filter.
%
% Copyright : Dmitry Ostrovsky, 2016.

%% handling params
if ~isfield(params,'sce'), error('Scenario not given.'); end
sce = params.sce;
if ~isfield(params,'snr'), params.snr = 1; end
snr = params.snr;
if ~isfield(params,'lep'), params.lep = 1; end
lep = params.lep;
if ~isfield(params,'mode'), params.mode = 1; end
mode = params.mode;
if ~isfield(params,'par'), params.par = 0; end
par = params.par;
if ~isfield(params,'verb'), params.verb = 1; end
verb = params.verb;
if ~isfield(params,'neighb'), params.neighb = 1; end
neighb = params.neighb;
if ~isfield(params,'w'), params.w = 50; end
w = params.w;
if ~isfield(params,'n'), params.n = 132; end
n = params.n;
%%
isImg = ~(length(sce) >= 7 && strcmp(sce(1:7),'RanSin-')); % image file
if isImg
    x = double(imread(sce));
    x = x(1:n,1:n);
else % random sines
    k = str2num(sce(8:end));
    % generate a 2-d signal composed of k sines
    freqsX = rand(k,1);
    freqsY = rand(k,1);
%     freqsX = sqrt(3) + sqrt(5)*(1:k)';
%     freqsY = sqrt(3) + sqrt(7)*(1:k)';
    x = 0;
    for i = 1:k
        a = repmat(2*pi*1i*freqsX(i)*(1:n), n, 1);
        b = repmat(2*pi*1i*freqsY(i)*(1:n)', 1, n);
        % x = x + randn(1,1)*exp(a+b);
        x = x + exp(a+b);
    end
    x = real(x);
end
sigm = norm(x)/n/snr;
y = x + sigm * randn(n,n);
if ~isImg
    % normalize
    yy = conj(y).*y;
    power = mean(mean(yy));
    x = x/sqrt(power);
    y = y/sqrt(power);
    csf = 50;
end
%% generate missing control params by default
p = 2; rho = 4; eps = sigm*rho/2;
control = struct('rho',rho,'accuracy','abs','eps',eps);
if nargin == 2
    % initialize control params given on input
    fields = fieldnames(input_control);
    for i=1:numel(fields)
        control.(fields{i}) = input_control.(fields{i});
    end
end
control
%%
% in the absolute accuracy mode Lepski starts from w = 2
% in the relative accuracy mode Lepski starts from w = 6
if ~lep
    if mode == 0
        rec_filt = recover_pointwise(y,w,control,par,verb,neighb);
    elseif mode == 1
        rec_filt = recover_blockwise(y,w,control,par,verb,neighb);
    else 
        rec_filt = recover_overblock(y,w,control,par,verb,neighb);
    end
else
    rec_filt = recover_lepski(y,sigm,mode,control,par,verb,neighb);
end
rec_filt = real(rec_filt);

% run Recht's lasso
overf=2;
ssize=size(y);
osize=ssize*overf;
% penalty "as usual"
lambda=sigm*sqrt(2*log(ssize(1)*ssize(2)))/4;
rdata=make_rdata(y, osize);
clear control
control.ItrMax=100;
verb=0;
control.print=verb; 
control.Eucl=1;
sol=fncl1mo(rdata,lambda,control);
rec_lasso=real(raha_oracle(sol.x,rdata));

close all
figure
param_string = ['sce-' num2str(sce) '_snr-' num2str(snr)...
    '_n-' num2str(n) '_w-' num2str(w) '_p-' num2str(p)... 
    '_rho-' num2str(rho) '_lep-' num2str(lep) '_mode-' num2str(mode)];
param_string = strrep(param_string, '.', 'p');
param_string = strrep(param_string, '/', ':');
set(gcf,'name', param_string)
subplot(2,2,1);
if isImg
    imshow(uint8(x));
else
    image(csf*real(x));
end
title('Ground truth');
subplot(2,2,2);
if isImg
    imshow(uint8(real(y)));
else
    image(csf*real(y));
end
title('Observations');
subplot(2,2,3);
if isImg
    imshow(uint8(real(rec_filt)));
else
    image(csf*real(rec_filt));
end
title(['Filtering recovery, ' 'MSE=' num2str(norm(rec_filt(:)-x(:)))]);
subplot(2,2,4);
if isImg
    imshow(uint8(rec_lasso));
else
    image(csf*rec_lasso);
end
title(['Lasso recovery, ' 'MSE=' num2str(norm(rec_lasso(:)-x(:)))]);
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