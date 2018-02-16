format long
SNR=1; 
m=1000;
stype=1;
[y,s]=ysim(m,SNR,stype);

osize=2000;
lambda=2;
data=make_rdata(y, osize);
%%%
control.ItrMax=100;
control.print=16; 
control.Eucl=0;
sol=fncl1mo(data,lambda, control);
s_hat=real(raha_oracle(sol.x,data));
t=0:m-1;
plot(t,y,'+',t,s,'-',t,s_hat,'m-.')
disp(['MSE=',num2str(norm(s_hat-s)),'; SNR=',num2str(norm(s)/norm(s_hat-s))]) 