function [ y ] = wl_filter(B, A, B_bar, A_bar, x )
% General Filter function for complex valued signals
% Creates and Widely Linear ARMA(p,q) signal where p is the order of the 
% Auto-Regressive part and q os the order of the Moving Average part. 
% ^* denotes the conjugate operator. 
% 
% a(1)*y(n) = b(1)*x(n) + ... + b(q)*x(n-q+1) - a(2)*y(n-1) - ... - a(p)*y(n-p+1) 
%             + b_bar(1)*x^*(n)+...+ b_bar(q_bar)*x^*(n-q_bar+1)  
%                 - a_bar(2)*y^*(n-1) - ... - a_bar(q_bar)*y^*(n-q_bar+1)
%
% Note that q (or p) need not be the same as q_bar (or p_bar). 
% The function takes care of it by appending the shorter vector with zeros.
%
% Hint: Set B_bar and A_bar = 0 to create a Linear ARMA model.
%
% ________________Important__________________________
% B, A, B_bar, A_bar need to be row vectors (1 x P)
% x needs to be a column vector (N x 1)
%
% Author: Sithan Kanna, Danilo Mandic, Imperial College London, April 2013
%________________________________________________________________________%




N = length(x); 

% Finds the Order of the Widely Linear ARMA model
p = max(length(A), length(A_bar)) -1; 
q = max(length(B), length(B_bar)); 



input = [zeros(q-1, 1); x]; 
output = zeros(N+p, 1); 



B =  [zeros(q-length(B),1); transpose(B(end:-1:1))]; 
B_bar =  [zeros(q-length(B_bar),1); transpose(B_bar(end:-1:1))]; 
A =  [zeros(p-length(A)+1,1); transpose(A(end:-1:1))]; 
A_bar =  [zeros(p-length(A_bar)+1,1); transpose(A_bar(end:-1:1))]; 



for i = 1:N
    output(p+i) = B'*input(i:i+q-1) -A'*output(i:i+p) ...
        + B_bar'*conj(input(i:i+q-1)) -  A_bar'*conj(output(i:i+p));
end
y = output(p+1:end); 

end

