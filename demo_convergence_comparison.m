clear all;
close all;

addpath('plotting_function');
addpath('func_LRTV');


%%% Original Data %%%
N = 128; % resolution of demo image

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


% initial values of step size
gam1 = 0.01;
Gam  = [gam1 1/(8*gam1)];
maxiter = 10000; % maximum iterations
tol = 1e-8;      % parameter for checking convergence
verb = 1;        % visualizing optimization procedure


  %fix step size
  tic;
  [X, histo, GAM1, GAM2, ppk, ddk, histo_obj] = LRTV_pds_fix(Tms/dom(2),Q,al,tv,be,lr,delta2,dom/dom(2),Gam,maxiter,tol,verb);
  FIX.X    = X;
  FIX.time = toc;
  FIX.histo= histo;
  FIX.GAM1 = GAM1;
  FIX.GAM2 = GAM2;
  FIX.ppk  = ppk;
  FIX.ddk  = ddk;
  FIX.obj  = histo_obj;

  %Goldstein's adaptation
  tic;
  [X, histo, GAM1, GAM2, ppk, ddk, Cbk, histo_obj] = LRTV_pds_nips(Tms/dom(2),Q,al,tv,be,lr,delta2,dom/dom(2),Gam,maxiter,tol,verb);
  NIPS.X    = X;
  NIPS.time = toc;
  NIPS.histo= histo;
  NIPS.GAM1 = GAM1;
  NIPS.GAM2 = GAM2;
  NIPS.ppk  = ppk;
  NIPS.ddk  = ddk;
  NIPS.Cbk  = Cbk;
  NIPS.obj  = histo_obj;

  %Proposed adaptation
  rho = 0.05;
  tic;
  [X, histo, GAM1, GAM2, ppk, ddk, Cp, Cd, histo_obj] = LRTV_pds_ada(Tms/dom(2),Q,al,tv,be,lr,delta2,dom/dom(2),Gam,rho,maxiter,tol,verb);
  ADAP.X    = X;
  ADAP.time = toc;
  ADAP.histo= histo;
  ADAP.GAM1 = GAM1;
  ADAP.GAM2 = GAM2;
  ADAP.ppk  = ppk;
  ADAP.ddk  = ddk;
  ADAP.Cp   = Cp;
  ADAP.Cd   = Cd;
  ADAP.obj  = histo_obj;


% 1) Execution time comparison
figure(1); clf;
times = [FIX.time, NIPS.time, ADAP.time];
bar(times);
set(gca,'XTickLabel',{'FIX','NIPS','ADAP'});
ylabel('time (s)');
title('Execution time');

% 2) Reconstructed image comparison
figure(2); clf;
subplot(1,3,1); imagesc(FIX.X); axis image off; colormap gray; title('FIX');
subplot(1,3,2); imagesc(NIPS.X); axis image off; colormap gray; title('NIPS');
subplot(1,3,3); imagesc(ADAP.X); axis image off; colormap gray; title('ADAP');

% 3) Objective function history (log scale)
figure(3); clf; hold on;
plot(FIX.obj,'-','LineWidth',1.4);
plot(NIPS.obj,'-','LineWidth',1.4);
plot(ADAP.obj,'-','LineWidth',1.4);
legend('FIX','NIPS','ADAP');
set(gca,'yscale','log');
xlabel('iteration'); ylabel('objective (log)');
title('Objective history'); grid on;
hold off;
axis([0, length(FIX.obj), FIX.obj(end)/1.5, FIX.obj(1)]);


% 4) condition / ppk / ddk comparison
figure(4); clf;
subplot(3,1,1); hold on;
plot(FIX.histo,'-'); plot(NIPS.histo,'-'); plot(ADAP.histo,'-');
set(gca,'yscale','log');
legend('FIX','NIPS','ADAP'); title('Condition'); xlabel('iteration'); ylabel('condition'); grid on; hold off;

subplot(3,1,2); hold on;
plot(FIX.ppk,'-'); plot(NIPS.ppk,'-'); plot(ADAP.ppk,'-');
set(gca,'yscale','log');
legend('FIX','NIPS','ADAP'); title('Primal residual'); xlabel('iteration'); ylabel('primal'); grid on; hold off;

subplot(3,1,3); hold on;
plot(FIX.ddk,'-'); plot(NIPS.ddk,'-'); plot(ADAP.ddk,'-');
set(gca,'yscale','log');
legend('FIX','NIPS','ADAP'); title('Dual residual'); xlabel('iteration'); ylabel('dual'); grid on; hold off;
