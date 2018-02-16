function [f,df]=fdf_oracle(x,rdata)
no=nargout;
Ax=feval(rdata.mult,x,rdata,0);
temp=Ax-rdata.y;
ffttemp = fft(temp)*sqrt(length(temp));
f=norm(temp(:),2)^2/2;
if no>1,
    df=feval(rdata.mult,temp,rdata,1); end
end % endof fdf_oracle
