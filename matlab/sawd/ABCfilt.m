function [A, B, C, B1, D, a, b, b1, H1norm, HTH1norm] = ABCfilt(deg, fc, N, K)
% [A, B, B1] = ABfilt(d, fc, N, K)
%
% Banded matrices for zero-phase high-pass recursive filter.
% The filter is H = inv(A) * B where the matrices A and B
% are created as 'sparse' matrices.
%
% INPUT
%   d  : degree of filter is 2d
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal
%   K  : oder of difference matrix D (need 1 <= K <= 2*d) (default K = 1)
%
% OUTPUT
%   A, B, B1 : banded filter matrices
%       with B = B1*D where D is the K-th order difference (up to sign)
%
% Use [A, B, B1, D, a, b, b1] = ABfilt(...) to return
% filter coefficient vectors a, b, b1.

% Ivan Selesnick,  NYU-Poly, 2012
% Version 2, 2016
% Version 3, 2017, August
%
% Revised in 2016 to reduce transient effects at the start and the end
% of the signal, in comparison with the original 2012 version.
%
% Revised in 2017 to allow 1 <= K <= 2*d

if nargin < 4
    K = 1;
end

if K > 2*deg
    error('ABfilt: K > 2*d')
end

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^deg;

% Define p such that P(z)P(1/z) = B(z), i.e., P'*P = B
p = 1;
for k = 1:deg
    p = conv(p, [-1 1]);
end
P = spdiags( p(ones(N-deg,1), :), 0:deg, N-deg, N);        % banded matrix

B = P' * P;
% pp = conv(p, p(end:-1:1));

q = sqrt(t);
for i = 1:deg
    q = conv(q, [1 1]);
end
Q = spdiags( q(ones(N-deg,1), :), 0:deg, N-deg, N);    % banded matrix

A = P'*P + Q'*Q;
C = Q'*Q;

if K <= deg
    d = 1;
    for i = 1:K
        d = conv(d, [-1 1]);
    end
    D = spdiags(d(ones(N,1), :), 0:K, N-K, N);               % D: banded matrix
    
    p1 = deconv(p, d);
    P1 = spdiags( p1(ones(N-deg,1), :), 0:deg-K, N-deg, N-K);  % banded matrix
    B1 = P' * P1;

    b1 = conv(p1, p(end:-1:1));
else
    % deg < K <= 2*deg
    K2 = 2*deg - K;
    d = 1;
    for i = 1:K2
        d = conv(d, [-1 1]);
    end
    B1 = spdiags(d(ones(N,1), :), 0:K2, N-K2, N)';               % D: banded matrix
    
    p1 = deconv(p, d);
    D1 = spdiags( p1(ones(N-deg,1), :), 0:deg-K2, N-deg, N-K2);  % banded matrix
    D = D1'*P;

    b1 = d;

end

a = conv(p, p(end:-1:1)) + conv(q, q(end:-1:1));
b = conv(p, p(end:-1:1));

% verify that B = B1*D
err = B - B1*D;
mae = max(abs(err(:)));
if mae > 1e-10
    disp('Error in ABfilt (B1*D not equal to B)')
end


% Calculate filter norms

imp = zeros(size(B1,2), 1);
imp(round(N/2)) = 1;                        % imp : impulse signal (located at center to avoid transients)

h1 = A \ (B1 * imp);
H1norm = sqrt( sum( abs( h1 ).^2 ) );       % norm of filter inv(A)*B1

hh = B' * ((A*A') \ (B1 * imp));
HTH1norm = sqrt( sum( abs( hh ).^2 ) );     % norm of filter B'*inv(A*A')*B1


