%% Random spikes
% Generate some 1-D data: 100 samples of the sum of 3 random sines, SNR=2.
close all
rng(2,'twister'); % initialize random number generator
n=100;
[x,y,sigm] = generate_data('RanSin-3',n,2);
%% 
% Denoise via Recht's oversampled Lasso.
tic; recl = lasso_recovery(y,sigm); toc
%%
% Denoise via filtering, take rho = 6 (twice the number of frequencies).
clear params
params.sigm=sigm; % mandatory
params.rho=6; % mandatory, twice the number of frequencies
params.lep=0; % no bandwidth adaptation.
tic; recf = filter_recovery(y,params); toc
%%
% Plot the results (the plot is saved in '../test/rec_test-1d/').
plot_recovery_results(params,x,y,recf,recl);
%%
% Lasso performed better... Well, let's change to the pointwise mode. Note
% that the sample will be effectively 2 times smaller.
params.mode=0;
tic; recf = filter_recovery(y,params); toc
plot_recovery_results(params,x,y,recf);
%%
% Better now, we've almost catched up. But note that eps is 0.45 when MSE
% is 1.11. Let's try increasing the precision.
params.verb=0; % quiet mode
solver_control = struct('accuracy','rel','eps',0.05); % relative accuracy 0.05
tic; recf = filter_recovery(y,params,solver_control); toc
norm(recf-x)
%%
% Didn't help much. Actually that means we've already found a very precise
% solution the last time. Let's move on.
%% Slow modulation
n = 1500;
[x,y,sigm] = generate_data('6',n,1);
recl = lasso_recovery(y,sigm);
clear params
params.sigm = sigm; params.rho = 2; params.mode = 1; params.lep = 0;
recf = filter_recovery(y,params);
plot_recovery_results(params,x,y,recf,recl);
%%
% Note that Filtering shrinked the signal a bit too much.
plot_recovery_results(params,x,y,recf);
%%
% We can correct it by increasing rho.
params.rho = 8;
recf = filter_recovery(y,params);
plot_recovery_results(params,x,y,recf);
%% Non-stationary signal
n = 200; 
[x,y,sigm] = generate_data('7',n,4);
clear params
params.sigm = sigm;
params.rho = 2;
params.mode = 0;
params.lep = 1; % bandwidth adaptation
recf = filter_recovery(y,params);
plot_recovery_results(params,x,y,recf);
%%
% Not that good but we managed to grasp the behaviour of the signal.
%% Speech
n = 4000;
[x,y,sigm] = generate_data('Speech',n,0.7);
clear params
params.sigm = sigm;
params.rho = 4; % ~2 harmonics are enough for a human ear
params.mode = 1;
params.lep = 1;
params.par = 1; % let's go parallel
solver_control = struct('accuracy','rel','eps',0.2);
recf = filter_recovery(y,params);
recl = lasso_recovery(y,sigm);
plot_recovery_results(params,x,y,recf,recl);
norm(y-x)
%%
% However by ear one can hear some difference.
Fs=7418;
spath = '../test/rec_test-1d/speech/';
if ~exist(spath, 'dir'), mkdir(spath); end
% Manually normalizing to supress the 'data clipped' warning
x = x/max(abs(x));
y = y/max(abs(y));
recl = recl/max(abs(recl));
recf = recf/max(abs(recf));
% Writing the output
audiowrite([spath 'true.wav'],x,Fs,'BitsPerSample',16);
audiowrite([spath 'noised.wav'],y,Fs,'BitsPerSample',16);
audiowrite([spath 'lasso.wav'],recl,Fs,'BitsPerSample',16);
audiowrite([spath 'filter.wav'],recf,Fs,'BitsPerSample',16);