function [z, circ] = noise_gen(CovMat, Sig_Len)
% Creates a complex valued noise  z = x + jy of length Sig_Len and
% calculates  the circularity coefficent.
%
% The variance and cross-correlation values of x and y are given by CovMat
% CovMat = |sig_x^2   sig_xy  |
%          |sig_xy    sig_y^2 |
% It accomplishes it by starting from 2 independent zero mean Gaussian 
% variables with unit variance a and b. 
% x = R*a, y = I*b + C*a, where R, I and C are constants used to set the
% powers and cross correlations of x and y. 
%
% sig_x2   = E{x^2} = (R^2)*E{a^2} = R^2 (eq1)
% sig_y2   = E{y}   = I^2 + C^2          (eq2)
% sigma_xy = E{xy}  = (RC)*E{a^2} = RC   (eq3)
% 
% Using (eq1),(eq2) and (eq3) the values for the constant are computed. 
% 
%
% Author: Sithan Kanna, Danilo Mandic, Imperial College London, April 2013
%________________________________________________________________________%

sig_x2  = CovMat(1,1);
sig_y2 = CovMat(2,2);
sig_xy = CovMat(1,2);

R = sqrt(sig_x2);
C = sig_xy/sqrt(sig_x2);
I = sqrt(sig_y2 - C^2);

CovMat_DC= [R, 0; C, I];

dual_channel = CovMat_DC*randn(2,Sig_Len);
z = dual_channel(1, :) + 1i*dual_channel(2,:);
z = z';
circ = (sig_x2  - sig_y2 + 2*1i*sig_xy)/(sig_x2 + sig_y2);
end

