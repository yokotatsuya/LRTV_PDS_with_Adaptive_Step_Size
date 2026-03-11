clear all

addpath('plotting_function');
addpath('func_LRTV');


%%% Original Data %%%
N = 256; % resolution of demo image

% Loading the RGB image as a 3rd order tensor
X0 = double(imread(['test_images/Giants/cgiant' num2str(N) '.png']));

% Domain of pixel values, we assume [0, 255]
dom(2) = 255;
dom(1) = 0;

%%% Making missing entry and adding noise for synthetic observed tensor %%%
II = size(X0);
N  = prod(II);
missing_rate = 0.3;

sig = 10; % controlling noise level by sig
noise = sig*randn(II); % generating noise

snr = SIR(X0,X0+noise); % evaluate SN ratio

idd = (randperm(N) > N*missing_rate);
Q   = reshape(idd,II);
Tms   = zeros(II);
Tms(Q)= X0(Q) + noise(Q);

%%% Recovery by using LRTVadaPDS %%%
% parameters for modeling
al = 0.5;           % weight for TV
tv = [0.5 0.5 0];   % TV level in each mode
be = 1 - al;        % weight for LR
lr = [0.4 0.4 0.2]; % LR level in each mode  
delta2 = 0.8 * (1 - missing_rate) * (sig/dom(2)).^2 * prod(II); % controling noise inequality

% parameters for optimization
Gam = [0.1 0.1]; % initial values of step size
rho = 0.05;      % hyper-parameter for step size adaptation
maxiter = 10000; % maximum iterations
tol = 1e-5;      % parameter for checking convergence
verb = 1;        % visualizing optimization procedure

% It skips to evaluate objective function for each step that requires several times of SVD
[X, histo, GAM1, GAM2, ppk, ddk, Cp, Cd] = LRTV_pds_ada_skipobj(Tms/dom(2),Q,al,tv,be,lr,delta2,dom/dom(2),Gam,rho,maxiter,tol,verb);

% It does not skip to evaluate objective function for each step that requires several times of SVD
%[X, histo, GAM1, GAM2, ppk, ddk, Cp, Cd, histo_obj] = LRTV_pds_ada(Tms/dom(2),Q,al,tv,be,lr,delta2,dom/dom(2),Gam,rho,maxiter,tol,verb);

figure(1); clf;
subplot(1,3,1)
imagesc(uint8(X0));
title('Original')

subplot(1,3,2)
imagesc(uint8(Tms));
title('Imcomplete&noisy')

subplot(1,3,3)
imagesc(uint8(X*dom(2)));
title(['LRTV (' num2str(SIR(X0,X*dom(2))) ' dB)'])