function [ w_l, W_A ] = wiener( input, output, Filt_Len )
% Estimates the optimal Linear (w_l) and Widely Linear (W_A) Wiener Coefficents
% for a set of input and outputs.
%
% (.)^H denotes the Hermitian Operator and ^* denotes the complex conjugate
% Given a true system that is given by y = (h^H)x +(g^H)x^* q where q is
% measurement noise that is white. The output y is a weighted sum of the
% input x and and x^* 
% 
% The Optimal Linear Wiener Coefficents are given by w_l = C^(-1)r 
% where C is the autocorrelation matrix C = E{xx^H} and 
% r is the cross correlation vector between the input and output, r = E{xy^*}
%
% The Optimal Widely Linear Wiener Coefficents are given by C_a^(-1)r_a
%  z = [x x^*] C_a = E{zz^H} and r_a E{zy^*}
%
%
% Author: Sithan Kanna,  Imperial College London, April 2013
%________________________________________________________________________%

[c, lags] = xcov(input, 'unbiased');
[p, ~] = xcov(input, conj(input), 'unbiased');
[r_cross, ~] = xcov(output, input, 'unbiased');
[u_cross, ~] = xcov (conj(output), input, 'unbiased');

dom = and(lags>= 0, lags< Filt_Len); 

C = toeplitz(c(dom));
P = toeplitz(real(p(dom))) + 1i.*toeplitz(imag(p(dom))); 
r = conj(r_cross(dom));
u = conj(u_cross(dom)); 

w_l = (C)^(-1)*(r); 
h_o = (C - P*conj(C)^(-1)*conj(P)) \ (r - P*conj(C)^(-1)*conj(u)); 
g_o = (conj(C) - conj(P)*C^(-1)*P) \ (conj(u) - conj(P)*C^(-1)*r); 

W_A = [h_o; g_o];
 
end

