function [f,df]=r_oracle(x,rdata)
no=nargout;
Ax=raha_oracle(x,rdata);
temp=Ax-rdata.y; 
f=norm(temp(:),2)^2;
if no>1, df=raha_oracle(temp,rdata,1); end
end % endof r_oracle
