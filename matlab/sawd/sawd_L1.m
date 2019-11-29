%==========================================================================
% function [x, cost, k, k_app ] = sawd_L1(y, d, fc, lam, mu, k_app)
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB function is run the SAWD
%                application
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019, Released under the 3-Clause BSD license
%==========================================================================
function [x, cost, k, k_app ] = sawd_L1(y, d, fc, lam, mu, k_app)

global simdata;
fs = simdata.fs;

y = y(:);                               % Convert to column vector
N = length(y);   
[A, ~, C] = ABCfilt(d, fc, N);          % Create the low-pass filter
G = (mu*(A*A') + (C*C'));               % G: Banded system
L = @(x) A\(C*x);                       % L: low-pass filter

MAX_ITER = 50;
TOL_STOP = 1e-4;

cost = zeros(1, MAX_ITER);              % Cost function history

[A1,A1H] = make_transforms('DWT2',N, [2^nextpow2(fs) 1 8 2^nextpow2(fs)]); 

k = A1H(L(y));
d = zeros(size(k));
b = (1/mu) * C' * ((A*A')\(C*y));
Ab = A1H(b);

iter = 0;
last_iter = false;


while not(last_iter)
    
    iter = iter + 1;
    
    g = Ab + (k+d);
    u = g - A1H(C' * (G\(C*(A1(g)'))));
    k = soft(u-d,lam/(mu));
    d = d - (u-k);
    cost(1,iter) = 0.5*sum(abs(L(y - real(A1(k)'))).^2) + ...
        lam * sum(abs(k(:))) ;
    
    if (iter >= MAX_ITER) || ...
            ((iter>2) && (abs(cost(1,iter)-cost(1,iter-1)) < TOL_STOP))
        last_iter = true;
    end
    
end
cost = cost(1:iter);

x = real(L(A1(k)'));
x = x(:);
k_app = [k(:,2), k_app];


end

