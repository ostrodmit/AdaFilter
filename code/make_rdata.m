function data=make_rdata(y,lsize)
% prepares data for the filter
ssize=size(y); 
if min(ssize(1),ssize(2))==1, data.ssize=ssize(1)*ssize(2); 
else data.ssize=size(y); 
end
if isscalar(lsize), lsize=[lsize,1]; end
if min(lsize(1),lsize(2))==1, data.lsize=lsize(1); 
else data.lsize=lsize; end

if (lsize(1)<ssize(1))||(lsize(2)<ssize(2))
    error('large_size<small_size (undersampling is not allowed)')
end

data.n=lsize(1)*lsize(2);
data.y=y;
data.mult=@raha_oracle;
data.Lf=1;
end % endof make_rdata
