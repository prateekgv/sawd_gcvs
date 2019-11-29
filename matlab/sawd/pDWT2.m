%==========================================================================
% function [toe_off,heel_stk] = sawd_aided_gait_segmentation(u, zupt, foot)
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is to compute windowed 
%                wavelet transform
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019, Released under the 3-Clause BSD license
%==========================================================================
function X = pDWT2(x,R,M,p,Nfft)

if ~isrow(x)
    x = x';
end

[h,g] = compute_wavelet_filter('Daubechies',p);
% disp(['h filter = [' num2str(h) ']']);
% disp(['g filter = [' num2str(g) ']']);

x = [zeros(1,R) x zeros(1,R)];          % to deal with first and last block
Y = buffer(x,  R, R*(M-1)/M, 'nodelay');
% Y = buffer(x,  R, R*(M-1)/M);
X = zeros(size(Y));
[~,n] = size(X);

Jmax = log2(Nfft)-1; Jmin = 0;
for i=1:n
    X(:,i) = Y(:,i);
    for j=Jmax:-1:Jmin
        a1 = X(1:2^(j+1),i);
        a = subsampling(cconvol(a1,h));
        d = subsampling(cconvol(a1,g));
        X(1:2^(j+1),i) = cat(1, a, d );
    end

end

end

