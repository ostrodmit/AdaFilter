function [x,y,sigm] = generate_data2(sce,n,snr)
%
% GENERATE_DATA Generate noisy 2-D observations.
%
% Input parameters
%   *sce [{'<filename>','RanSin-<k>','CohSin-<k>','SingleIdx-<s>'}]: 
%       scenario. If 'RanSin-<k>' or 'CohSin-<k>', the artificial signal
%       (random 2-d sines) will be generated, <k> being the number of
%       spikes/pairs of spikes. If 'SingleIdx-<s>', a function which is
%       (<s>,1)-Hoelder smooth along one direction, and constant across,
%       will be generated. Otherwise, image file '<filename>' will be used.
%   *n [{int > 0}] : sample size.
%   snr [double,1] : signal-to-noise ratio (*not* in dBs).
%
% Output parameters
%   x : true signal.
%   y : observations.
%   sigm : noise level (standard deviation).
%
% Examples
%   [x,y,sigm] = generate_data2('../Brodatz/D75.gif',128,1); % Brodatz file
%
% See also : generate_data.
%
% Copyright : Dmitry Ostrovsky, 2016.
if nargin<3, snr = 1; end
isSingleIdx = 0;
if length(sce) >= 10,
    if strcmp(sce(1:10), 'SingleIdx-'),
        isSingleIdx = 1;
    end
end
isRanSin = length(sce) >= 7 && strcmp(sce(1:7),'RanSin-');
isCohSin = length(sce) >= 7 && strcmp(sce(1:7),'CohSin-');
isImg = ~(isSingleIdx || isRanSin || isCohSin); % image file
if isImg
    x = double(imread(sce));
    x = x(1:n,1:n);
else % RanSin, CohSin, or SingleIdx
    if isRanSin || isCohSin,
        k = str2num(sce(8:end));
        if isRanSin % k spikes in the spectral domain
            freqsX = rand(k,1);
            freqsY = rand(k,1);
            %     freqsX = sqrt(3) + sqrt(5)*(1:k)';
            %     freqsY = sqrt(3) + sqrt(7)*(1:k)';
            numfreqs = k;
        else % k pairs of close spikes in the spectral domain
            spikesX = rand(k,1);
            spikesY = rand(k,1);
            delta = 0.1;
            freqsX = [spikesX;spikesX+delta/n]; % close random spikes
            freqsY = [spikesY;spikesY+delta/n];
            numfreqs = 2*k;
        end
        amps = randn(numfreqs,1);
        x = 0;
        for j = 1:numfreqs
            a = repmat(2*pi*1i*freqsX(j)*(1:n), n, 1);
            b = repmat(2*pi*1i*freqsY(j)*(1:n)', 1, n);
            % x = x + randn(1,1)*exp(a+b);
            x = x + amps(j)*exp(a+b);
        end
    else
        s = str2num(sce(11:end));
        % generate a random function from the Sobolev ball
        R = 1/pi^s;
        N = 100*n; % oversampling factor to compensate for the edge effects
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
            c = [c;(trigc(k)-1i*trigc(k+1))/2];% complex Fourier basis
        end
        c = [0; c; flipud(conj(c))]; % symmetrize
        f = N*ifft(c); % values of the function over the N-grid
        %plot(f)
        Ngrid = linspace(0,1,N+1)';
        ngrid = linspace(0,1,n)';
        % generate a random direction
        theta = randn(2,1); theta = abs(theta) / norm(theta,1);
        x = zeros(n,n);
        for i = 1:n,
            for j = 1:n,
                % define offset
                t = [ngrid(i,1); ngrid(j,1)];
                offset = dot(theta,t);
                % find closest point of a finer grid
                tmp = abs(Ngrid-offset);
                [~,idx] = min(tmp);
                x(i,j) = f(idx);
            end
        end
    end
    x = real(x);
end
sigm = norm(x(:))/n/snr;
y = x + sigm * randn(n,n);
if ~isImg
    % normalize
    yy = conj(y).*y;
    power = mean(mean(yy));
    x = x/sqrt(power);
    y = y/sqrt(power);
    sigm = sigm/sqrt(power);
end