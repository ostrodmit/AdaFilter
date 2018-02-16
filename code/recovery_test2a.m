csf = 50; % color scaling factor
p = 2;

sigm = 1;
rho = 4; lambda = sigm/rho;

T = 100; L = 99;
n = T+L+1;
t = T+1:T+L+1;
k = 2;

% % generating a 2d signal comprised of k sines
% freqsX = rand(k,1);
% freqsY = rand(k,1);
% x=0;
% for i = 1:k
%     a = repmat(2*pi*1i*freqsX(i)*(1:n), n, 1);
%     b = repmat(2*pi*1i*freqsY(i)*(1:n)', 1, n);
%     x = x+randn(1,1)*exp(a+b);
% end
% y = x + sigm * randn(n,n);
% % normalize
% yy = conj(y).*y;
% power = mean(mean(yy));
% x = x/sqrt(power);
% y = y/sqrt(power);
% y1 = y;
% y2 = y1(t, t);
control = struct('p', p, 'constrained', 1, 'rho', rho, 'lambda', lambda,...
    'solver', 'mp',...
    'accuracy', 'rel', 'eps', 0.25, 'tol', 1e-10,...
    'max_iter',1000,'max_cpu', 1000, 'l2_prox', 0,...
    'warm_start', 1, 'online', 1, 'last_part', 0, 'verbose', 1);
sol = fit_filter(y1, y2, control);
rec = conv2(y1, sol.filter, 'valid');
figure
subplot(2,2,1);
image(csf*real(x(t,t)));
title('True signal');
subplot(2,2,2);
image(csf*real(y2));
title('Observations');
subplot(2,2,3);
image(csf*real(rec));
title('Recovery');
subplot(2,2,4);
ycut=randi(L+1); ycuT=T+ycut;
plot(1:L+1,real(x(ycuT,t)),'-g',1:L+1,real(rec(ycut,:)),'--b',1:L+1,real(y2(ycut,:)),'+r') 
title(['Cut for y=',num2str(ycut)]);