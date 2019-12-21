function [x, pdf] = WAHE( f, g, Threshold )

% -------------------------------------------------------------------------
% An implementation of "Weighted Approximated Histogram Equalization."
%   T. Arici, S. Dikbas, and Y. Altunbasak, "A histogram modification
%   framework and its application for image contrast enhancement," IEEE
%   Trans. Image Process., vol. 18, no. 9, pp. 1921-1935, Sep. 2009.
%
% -------------------------------------------------------------------------
% Input variables
%   f        : Input gray scale image, single channel
%   g        : Level of Enhancement
%   Threshold: See the paper for details
%
% Output variables
%   x    : Output transformation function.
%   pdf  : Probability mass function of input pixel intensities, or
%   equivalently normalized histogram vector.
% 
% -------------------------------------------------------------------------
%                           written by Chulwoo Lee, chulwoo@mcl.korea.ac.kr


if nargin < 3
    Threshold = 2;
end

[R,C] = size(f);


%% WAHE
kappa = 0;                      % initialize kappa
count = 0;                      % initialize count
h = zeros(256,1);
for m=1:R
    for n=3:C
        temp = f(m,n) - f(m,n-2);
        kappa = kappa + abs(temp);
        if abs(temp) > Threshold
            h(f(m,n)+1) = h(f(m,n)+1) + 1;
            count = count + 1;
        end
    end
end

kappa_star = 1/(1+g);
u = count/256;                  % uniform

% We do not apply B&W stretch.
h_tilda = zeros(256,1);
for n = 1:256
    h_tilda(n) = (1-kappa_star)*u + kappa_star * h(n);
end


%% reconstruct transformation function
pdf = h_tilda/sum(h_tilda);

cdf = zeros(length(h_tilda),1); cdf(1) = pdf(1);
for k=2:length(h_tilda)
    cdf(k) = cdf(k-1) + pdf(k);
end

x = 255*cdf;

    
end
