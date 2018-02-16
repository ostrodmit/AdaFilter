function x = recover_pointwise(y,w,control,parallel,verb,neighb)
if nargin < 6; neighb = 1; end
if nargin < 5; verb = 1; end
if nargin < 4; parallel = 0; end
if parallel
    x = recover_pointwise_par(y,w,control,verb,neighb);
else
    x = recover_pointwise_seq(y,w,control,verb,neighb);
end