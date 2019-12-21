function x = LDR( src_data, alpha, U )

% -------------------------------------------------------------------------
% An implementation of
%   C. Lee, C. Lee, and Chang-Su Kim, "Contrast enahancement based on
%   layered difference representation of 2D histograms," IEEE Trans. Image
%   Image Process., vol. 22, no. 12, pp. 5372-5384, Dec. 2013
%
% -------------------------------------------------------------------------
% Input variables (see the paper for details)
%   src_data : can be either 2D histogram or gray scale image. This script
%   automatically detects based on its dimension.
%   alpha    : controls the level of enhancement
%   U        : U matrix in Equation (31). If it is provided, we can save
%   the computation time.
%
% Output variables
%   x    : Output transformation function.
% 
% -------------------------------------------------------------------------
%                           written by Chulwoo Lee, chulwoo@mcl.korea.ac.kr


if nargin < 3
    % Pre-computing
    U = zeros(255,255);
    tmp_k = 1:255;
    for layer=1:255
        U(:,layer) = min(tmp_k,256-layer) - max(tmp_k-layer,0);
    end
end

if nargin < 2
    alpha = 2.5;
end


[R, C] = size(src_data);
if R==256 && C==256
    h2D_in = src_data;
else
    in_Y = src_data;
    
    % unordered 2D histogram acquisition
    h2D_in = zeros(256,256);
    
    for j=1:R
        for i=1:C
            ref = in_Y(j,i);
            
            if j~=R
                trg = in_Y(j+1,i);
                h2D_in(max(trg,ref)+1,min(trg,ref)+1) = h2D_in(max(trg,ref)+1,min(trg,ref)+1) + 1;
            end
            
            if i~=C
                trg = in_Y(j,i+1);
                h2D_in(max(trg,ref)+1,min(trg,ref)+1) = h2D_in(max(trg,ref)+1,min(trg,ref)+1) + 1;
            end
        end
    end
    clear ref trg
end


%% Intra-Layer Optimization
D = zeros(255,255);
s = zeros(255,1);


% iteration start
for layer = 1:255

    h_l = zeros(256-layer,1);

    tmp_idx = 1;
    for j=1+layer:256
        i=j-layer;

        h_l(tmp_idx,1) = log(h2D_in(j,i)+1);    % Equation (2)
        tmp_idx = tmp_idx+1;
    end
    clear tmp_idx

    s(layer,1) = sum(h_l);

    % if all elements in h_l is zero, then skip
    if s(layer,1) == 0
        continue
    end
    
    % Convolution
    m_l = conv(h_l, ones(layer,1));             % Equation (30)
    
    d_l = (m_l - min(m_l))./U(:,layer);         % Equation (33)

    if sum(d_l) == 0
        continue
    end
    D(:,layer) = d_l/sum(d_l);

end

%% Inter-Layer Aggregation
W = (s/max(s)).^alpha;                          % Equation (23)
d = D*W;                                        % Equation (24)

%% reconstruct transformation function
d = d/sum(d);       % normalization
tmp = zeros(256,1);
for k=1:255
    tmp(k+1) = tmp(k) + d(k);
end

x = 255*tmp;

end