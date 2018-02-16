format long
SNR=1; 
stype=1;
n = 500; 
k = 1000; 
% sigma = 1;

% %generating a 2d signal comprised of k sines
% freqsX = rand(k,1);
% freqsY = rand(k,1);
% s=0;
% for i = 1:k
%     a = repmat(2*pi*1i*freqsX(i)*(1:n), n, 1);
%     b = repmat(2*pi*1i*freqsY(i)*(1:n)', 1, n);
%     s = s+randn(1,1)*exp(a+b);
% end
% sigma=norm(s,'fro')/n/SNR;
% y = s + sigma * randn(n,n);
% % normalize

osize=[500,500];
lambda=4;
data=make_rdata(y, osize);
%%%
control.ItrMax=200;
control.print=16; 
control.Eucl='n';
sol=fncl1mo(data,lambda,control);
s_hat=real(raha_oracle(sol.x,data));
figure
subplot(2,2,1);
imagesc(real(s));
title('Signal');
subplot(2,2,2);
imagesc(real(y));
title('Observations');
subplot(2,2,3);
imagesc(real(s_hat));
title('MP recovery');
subplot(2,2,4);
ycut=randi(n); 
plot(1:n,real(s(ycut,:)),'-g',1:n,real(s_hat(ycut,:)),'--b',1:n,real(y(ycut,:)),'+r') 
title(['signal cut for y=',num2str(ycut)]);
disp(['MSE=',num2str(norm(s_hat-s,'fro')),'; SNR=',num2str(norm(s,'fro')/norm(s_hat-s,'fro'))]) 