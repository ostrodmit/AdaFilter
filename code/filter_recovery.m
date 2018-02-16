function rec = filter_recovery(y,params,solver_control)
%
% FILTER_RECOVERY Denoising of a 1-D or 2-D signal by adaptive filtering.
%
% Input parameters
%   *y : signal to denoise (column vector or matrix)
%   *params : structure with the following fields:
%       *sigm [double] : noise level. If given, the filters are fitted with 
%           absolute accuracy (see input_control). Will be used if lep=1.
%       *rho [double] : l_1-norm of the filter DFT (normalized to be ~1).
%           Used only if solver_control.constrained=1 (i.e. by default), to
%           control precision.
%       lep [{0,1},1] : if 1, use the Lepski-like adaptation procedure.
%       mode [{0-2},1] : filtering window mode. 0 : pointwise mode,
%           the filter is set up separately for each point. 1 : blockwise
%           mode, signal is divided into blocks, each of the size twice the
%           filter bandwidth, and the blocks are processed separately. 2 :
%           same but the blocks overlap by half.
%       par [{0,1},0] : use parfor. Compatible with neighb=1 (see below).
%       verb [{0,1},1] : show the progress.
%       neighb [{0,1},1] : boost filter fitting by initializing from the
%           previously processed point/block. Compatible with par=1.
%       w [{int > 0},*] : filter bandwidth (not relevant if lep=1). By
%           default is floor(n/4) if mode=0, floor(n/2) if mode= 1, and
%           floor(n/3) if mode==2.
%   solver_control : structure of control parameters for FIT_FILTER
%       The user may precise just the desired fields, except for *rho*
%       which is initialized by params.rho). Missing fields will be 
%       initialized by default values (see FIT_FILTER) but with 
%       accuracy='abs' and eps=params.sigm*params.rho/2.
%
% Output parameters
%   rec : recovery.
%
% Example : denoise a Brodatz image, with bandwidth adaptation, blockwise
%   [x,y,sigm] = generate_data2('../Brodatz/D75.gif',128,1.5);
%   params = struct('sigm',sigm,'rho',4,'lep',1,'mode',1);
%   recf = filter_recovery(y,params);
%
% See also : fit_filter, lasso_recovery.
%
% Copyright : Dmitry Ostrovsky, 2016.

%% handling params
if ~isfield(params,'sigm'), error('params.sigm must be initialized'); end
sigm = params.sigm;
if ~isfield(params,'rho'), error('params.rho must be initialized'); end
rho = params.rho;
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
if ~isfield(params,'w'),
    if mode == 0, params.w = floor(length(y)/4); end
    if mode == 1, params.w = floor(length(y)/2); end
    if mode == 2, params.w = floor(length(y)/3); end
end
w = params.w;
%% generate control params by default
accuracy = 'abs';
eps = sigm*rho/2;
control = struct('accuracy',accuracy,'eps',eps);
if nargin == 3
    % initialize control params given on input
    fields = fieldnames(solver_control);
    for i=1:numel(fields)
        control.(fields{i}) = solver_control.(fields{i});
    end
end
control.rho = rho;
%control
%% denoising
% in the absolute accuracy mode Lepski starts from w = 2
% in the relative accuracy mode Lepski starts from w = 6
if ~lep
    if mode == 0
        rec = recover_pointwise(y,w,control,par,verb,neighb);
    elseif mode == 1
        rec = recover_blockwise(y,w,control,par,verb,neighb);
    else 
        rec = recover_overblock(y,w,control,par,verb,neighb);
    end
else
%     if ~isfield(params,'sigm')
%         error('Noise level must be given if lep=1, initialize params.sigm');
%     end
    rec = recover_lepski(y,sigm,mode,control,par,verb,neighb);
end
rec = real(rec);
end