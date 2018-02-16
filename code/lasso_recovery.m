function rec = lasso_recovery(y,sigm,params)
%
% FILTER_RECOVERY Denoising of a 1-D or 2-D signal by the Lasso.
%
% Input parameters
%   *y : signal to denoise (vector or matrix)
%   *sigm [double] : noise level.
%   params : structure with the fields:
%       verb [{0,1},1] : show the progress.
%
% Output parameters
%   rec : recovery.
%
% Example : denoise a Brodatz image
%   [x,y,sigm] = generate_data2('../Brodatz/D75.gif',128,1.5);
%   recl = lasso_recovery(y,sigm);
%
% Copyright : Dmitry Ostrovsky, 2016.

if nargin<2
    error('Noise level must be given.');
end
if nargin<3 || ~isfield(params,'verb'), params.verb = 0; end
verb = params.verb;
overf=2;
if isvector(y),
    ssize=length(y);
else
    ssize=size(y);
end
osize=ssize*overf;
% penalty "as usual"
%lambda=sigm*sqrt(2*log(prod(ssize)))/4; % best
%lambda=sigm*sqrt(2*log(prod(ssize)))/2;
lambda=sigm*sqrt(2*log(prod(ssize)));
rdata=make_rdata(y, osize);
clear control
control.ItrMax=4000;
control.print=verb;
control.Eucl=1;
sol=fncl1mo(rdata,lambda,control);
rec=real(raha_oracle(sol.x,rdata));
end