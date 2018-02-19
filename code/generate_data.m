function [x,y,sigm] = generate_data(sce,n,snr)
%
% GENERATE_DATA Generate noisy 1-D observations.
%
% Input parameters
%   *sce [{'RanSin-<k>','CohSin-<k>','Sobolev-<s>','Speech','1'-'8'}] : 
%       signal scenario. Here, <k> is the number of harmonics. 
%       Random Sines, Coherent Sines, Sobolev smooth signal, Speech
%       (file 'mtlb.mat'), or an artificial signal generated by YSIM.
%   *n [{int > 0}] : signal length.
%   snr [double,1] : signal-to-noise ratio (*not* in dBs).
%
% Output parameters
%   x : true signal.
%   y : observations.
%   sigm : noise level (standard deviation).
%
% Examples
%   [x,y,sigm] = generate_data('RanSin-3',100,2); % 3 random sines, snr=2
%   [x,y,sigm] = generate_data('Speech',4000,1); % speech
%
% See also : generate_data2.
%
% Copyright : Dmitry Ostrovsky, 2016.

if nargin<3, snr = 1; end
t = (0:n-1)';
if length(sce)>=7 ...
        && (strcmp(sce(1:7),'RanSin-') || strcmp(sce(1:7),'CohSin-'))
    k = str2num(sce(8:end));
    if strcmp(sce(1:7),'RanSin-')
        freqs = rand(k,1);
        numfreqs = k;
    else
        spikes = rand(k,1);
        delta = 0.1;
        freqs = [spikes;spikes+delta/n]; % close random spikes
        numfreqs = 2*k;
    end
    amps = randn(numfreqs,1);
    x = zeros(n,1);
    for j = 1:numfreqs
        x = x+amps(j)*exp(2*pi*1i*freqs(j)*t')';
    end
    x = real(x);
    sigm = norm(x)/sqrt(n)/snr;
    y = x + sigm*randn(n,1);
elseif length(sce)>=8 && strcmp(sce(1:8),'Sobolev-')
    s = str2num(sce(9:end));
    % generate a random function from the Sobolev ball
    R = 1/pi^s;
    N = 100*n; % finer grid
    semiax = ones(N,1);  % ellipsoid in the trigonometric basis
    for k = 1:N,
        if mod(k,2),
            semiax(k) = (k+1)^s;
        else
            semiax(k) = k^s;
        end
        semiax(k) = semiax(k)/R;
    end
    z = randn(N,1); z = z/norm(z);
    trigc = z./semiax; % trigonometric Fourier basis
    c = [];
    for k = 1:2:N-1,
        c = [c;(trigc(k)-1i*trigc(k+1))/2]; % complex Fourier basis
    end
    c = [0; c; flipud(conj(c))]; % symmetrize
    f = N*ifft(c); % values of a Sobolev-smooth function over the N-grid
    Ngrid = linspace(0,1,N+1)';
    plot(Ngrid,f)
    ngrid = linspace(0,1,n)';
    % generate a random direction
    x = zeros(n,1);
    for i = 1:n,
        t = ngrid(i,1);
        % find closest point of a finer grid
        tmp = abs(Ngrid-t);
        [~,idx] = min(tmp);
        x(i) = f(idx);
    end
    x = real(x);
    sigm = norm(x)/sqrt(n)/snr;
    y = x + sigm*randn(n,1);
elseif strcmp(sce,'Speech') % speech
    speech = open('mtlb.mat');
    x = speech.mtlb(1:n);
    sigm = norm(x)/sqrt(n)/snr;
    y = x + sigm*randn(n,1);
elseif strcmp(sce,'Piano') % piano stroke, downsampled from 44 to 11 kHz.
    [X,fs] = audioread('data/c1_11kHz.wav');
    x = X(:,1); % only left channel;
    % downsample
    x = x(1:min(length(x),n));
    sigm = norm(x)/sqrt(length(x))/snr;
    y = x + sigm*randn(length(x),1);
elseif ismember(str2num(sce),1:8) % artificial signal
    [y,x]=ysim(n,snr,str2num(sce));
    sigm = norm(x)/sqrt(n)/snr;
else
    error('Incorrect scenario');
end
end