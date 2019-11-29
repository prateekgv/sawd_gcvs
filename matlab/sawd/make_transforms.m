%==========================================================================
% main.m
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is to invoke windowed
%                wavelet transform and its inverse
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019, Released under the 3-Clause BSD license
%==========================================================================
function [A, AH] = make_transforms(type, N, params)
%==========================================================================
% [A, AH] = make_transforms(type, N, params)
% with A AH = I
%
% INPUT
%    N : signal length
%    type : type of transform can be 'DFT' or 'STFT'
%    params : parameters for transform
% 
%   For DWT2, params = [R, M, p, Nfft]
%       R = length of frame
%       M = window overlapping factor
%       p = type of Daubechies wavelet
%       Nfft = number of DFT coefficients, Nfft >= R
% OUTPUT
%    A, AH : function handles for transform
%       AH is the conjugate transpose of A
%       Note: input to AH must be a row vector of length N
%
% DWT2 : Windowed discrete wavelet transform
%==========================================================================
switch type
        
    case 'DWT2'
        
        R = params(1);
        M = params(2);
        p = params(3);
        Nfft = params(4);
        AH = @(x) pDWT2(x,R,M,p,Nfft);
        A = @(x) ipDWT2(x,R,M,p,N);

end

