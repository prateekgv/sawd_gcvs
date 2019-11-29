%==========================================================================
% function y = ipDWT2(X,R,M,p,N)
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is to compute inverse 
%                windowed wavelet transform
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019, Released under the 3-Clause BSD license
%==========================================================================
function y = ipDWT2(X,R,M,p,N)

[h,g] = compute_wavelet_filter('Daubechies',p);
% disp(['h filter = [' num2str(h) ']']);
% disp(['g filter = [' num2str(g) ']']);

[Nfft, Nc] = size(X);                   % get size
Z = zeros(size(X));

Jmax = log2(Nfft)-1; Jmin = 0;
for i=1:Nc
    f1 = X(:,i);
    for j=Jmin:Jmax
        a = f1(1:2^j);
        d = f1(2^j+1:2^(j+1));
        a = cconvol(upsampling(a,1),reverse(h),1);
        d = cconvol(upsampling(d,1),reverse(g),1);
        f1(1:2^(j+1)) = a + d;
    end
    Z(:,i) = f1;
end


y = zeros(1,R/M*(Nc+M-1));
i = 0;
for k = 1:Nc
    y(i + (1:R)) = y(i + (1:R)) + Z(:,k).';
    i = i + R/M;
end
y = (1/M).*y(R+(1:N));

end