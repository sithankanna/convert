addpath(genpath('wl_toolbox'));


%% Create White Noise With Different Covariances 
N = 1000; 
CovMat1 = [0.5  0.4; 0.4 0.5]; 
CovMat2 = [0.8  0.2; 0.2 0.2]; 
CovMat3 = [0.5  0; 0  0.5]; 
%CovMat1 =[1.5  0.9; 0.9 1.5]; 
%CovMat2 = [1.5  0.9; 0.9 1.5]; 
%CovMat3 = [1.5  0.9; 0.9 1.5]; 

[z1, circ1] = noise_gen(CovMat1, N); 
[z2, circ2] = noise_gen(CovMat2, N); 
[z3, circ3] = noise_gen(CovMat3, N); 
z = [z1; z2; z3]; 
cir_coeff = [circ1.*ones(N,1); circ2.*ones(N,1); circ3.*ones(N,1)]; 


%% Track it With The LMS
mu = 0.01;
Sig_Len = length(z); 
FiltLen = 1; 
w = (0.5 +  0.5*1i)*ones(FiltLen, Sig_Len); 
y =  zeros(Sig_Len, 1);
e =  zeros(Sig_Len, 1);

for k = 1:Sig_Len-1
     y(k) = conj(w(:, k))*z(k); 
     e(k) = conj(z(k)) - y(k); 
     w(:, k+1) =  w(:, k) + (mu)*conj(e(k))*z(k);    
end
%% Plot the figure - Synthetic Data
FontSize = 16;

clf
subplot(2, 1,1); 
plot(real((w)), 'Color', [0 0.4 0.8], 'LineWidth' , 1.5); 
hold on
plot(real(cir_coeff), '--', 'Color', [0.8 0 0], 'LineWidth' , 2.5); 
hold off
%ylabel('Re\{p/c\}'); 
ylim([-0.1   1.1])

title('Real Part of the Circularity Quotient', 'FontSize',FontSize)
legend('Estimated', 'True'); 
xtitle = xlabel('Sample, k', 'FontSize', FontSize ); % Create label
pos = get(xtitle,'pos'); % Read position [x y z]
set(xtitle,'pos',[2850 -0.26 1]) % Move label to right

set(gca,'FontSize',FontSize);

subplot(2, 1,2); 
plot(imag((w)), 'Color', [0 0.4 0.8], 'LineWidth' , 1.5); 
hold on
plot(imag(cir_coeff), '--', 'Color', [0.8 0 0],  'LineWidth' , 2.5); 
hold off
%ylabel('Im\{ p/c \}'); 
%ylim([-0.1   1.1])



xtitle = xlabel('Sample, k', 'FontSize', FontSize ); % Create label
pos = get(xtitle,'pos'); % Read position [x y z]
set(xtitle,'pos',[2850 -0.7 1]) % Move label to right



title('Imaginary Part of the Circularity Quotient', 'FontSize',FontSize)
legend1 = legend('Estimated', 'True'); 
set(gca,'FontSize',FontSize);


% 
% 











