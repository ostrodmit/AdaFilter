%% Random sines
% Generate 2-D data: sum of 2 random sines, SNR=0.4.
close all
rng(3,'twister'); % initialize random number generator
n=80;
[x,y,sigm] = generate_data2('RanSin-2',n,0.4);
%% 
% Denoise via Recht's oversampled Lasso.
tic; recl = lasso_recovery(y,sigm); toc
%%
% Denoise via filtering, take rho = 6 (twice the number of frequencies).
clear params
params.sigm=sigm; % mandatory
params.rho=8; % mandatory
params.lep=0; % no bandwidth adaptation.
tic; recf = filter_recovery(y,params); toc
%%
% Plot the results (the plot is saved in '../test/rec_test-2d/').
%
% isImg=0; % not an image file
% isMesh=0; % don't use mesh
% plot_recovery_results2(params,isImg,isMesh,x,y,recf,recl);
mode = 'mine'; % for generated signals
save_pics=0; % don't save separate figures
plot_recovery_results2(params,mode,save_pics,x,y,recf,recl);
%% Brodatz
n = 128;
[x,y,sigm] = generate_data2('../Brodatz/D75.gif',n,1.5); % Brodatz file
recl = lasso_recovery(y,sigm);
clear params
params.sigm=sigm;
params.rho=4;
params.lep=1; % we need bandwidth adaptation, the picture has sharp borders
params.par=1; % let's go parallel
params.mode=1; % blockwise recovery
tic; recf = filter_recovery(y,params); toc
%%
% Plot results
%
% isImg=1;
mode = 'img'; % for image files
% isMesh=0;
save_pics=0; % don't save separate figures
plot_recovery_results2(params,mode,save_pics,x,y,recf,recl);
%%
% Now let's try overlapping blocks to improve the visual quality:
params.mode=2; % overlapping blocks
tic; recf = filter_recovery(y,params); toc
plot_recovery_results2(params,mode,save_pics,x,y,recf,recl);