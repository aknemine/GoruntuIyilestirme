function x = CVC( src_data, weights, wSize, K )

% -------------------------------------------------------------------------
% An implementation of "Contextual and Variational Contrast enhancement."
%   T. Celik and T. Tjahjadi, "Contextual and variational contrast
%   enhancement," IEEE Trans. Image Image Process., vol. 20, no. 12, pp.
%   3431-3441, Dec. 2011.
%
% -------------------------------------------------------------------------
% Input variables (see the paper for details)
%   src_data : can be either 2D histogram or gray scale image. This script
%   automatically detects based on its dimension.
%   weights: [alpha, beta, gamma]
%   wSize  : local window size during 2D histogram acquisition
%   K      : intensity resolution. e.g., 256
%
% Output variables
%   x    : Output transformation function.
% 
% -------------------------------------------------------------------------
%                           written by Chulwoo Lee, chulwoo@mcl.korea.ac.kr


if nargin < 2
    weights = [1/3 1/3 1/3];
end
if nargin < 3
    wSize = 7;
end
if nargin < 4
    K = 256;
end

D = diag(ones(256,1)) - diag(ones(255,1),1);
wHSize = floor(wSize/2);

% weights
if sum(weights) ~= 1
    weights = weights/sum(weights);
end
alpha = weights(1); beta = weights(2);  gamma = weights(3);


%% h_x: input 2D histogram
[R, C] = size(src_data);
if R==256 && C==256
    h_x = src_data;
else
    in_Y = src_data;
    
    % 2D histogram acquisition
    h_x = zeros(K,K);
    
    for j=1+wHSize:R-wHSize
        for i=1+wHSize:C-wHSize
            
            ref = in_Y(j,i);
            
            for jj=-wHSize:wHSize
                for ii=-wHSize:wHSize
                    if jj == 0 && ii ==0
                        continue
                    end
                    trg = in_Y(jj+j,ii+i);
                    h_x(ref+1,trg+1) = h_x(ref+1,trg+1) + 1;
                end
            end
        end
    end
    clear rer trg
end


%% h_p in Eq. (3) and Fig. 2(a)
[tmpx, tmpy] = meshgrid(1:K, 1:K);
h_p = ( sqrt( (tmpx-tmpy).*(tmpx-tmpy) ) + 1)/K;


%% Weighted h_x by h_p
h_x = h_x.*h_p;
H_x = h_x/sum(h_x(:));      % normalize
P_x = zeros(256,1);         % CDF
P_x(1) = H_x(1,1);
for k=2:K
    P_x(k) = P_x(k-1) + sum(H_x(k,1:k)) + sum(H_x(1:k,k)) - H_x(k,k);
end


%% H_u
H_u = ones(K,K)/(K*K);


%% H_t
H_t = ((alpha+beta)*eye(K)+ gamma*D*D')\(alpha*H_x+beta*H_u);
H_t = H_t/sum(H_t(:));


%% P_t
P_t = zeros(256,1);
P_t(1) = H_t(1,1);
for k=2:K
    P_t(k) = P_t(k-1) + sum(H_t(k,1:k)) + sum(H_t(1:k,k)) - H_t(k,k);        
end


%% mapping function
x = zeros(K,1);
for k=1:K
    min_idx = 0;
    min_val = 10;
    for l=1:K
        tmp_val = abs(P_t(l)-P_x(k));
        if tmp_val < min_val
            min_val = tmp_val;
            min_idx = l;
        end
    end
    x(k) = min_idx-1;
end


end

