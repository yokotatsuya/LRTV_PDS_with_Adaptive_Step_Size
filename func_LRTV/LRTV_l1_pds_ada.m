function [X, histo, GAM1, GAM2, PPk, DDk, Cp, Cd, histo_obj] = LRTV_l1_pds_ada(T,Q,al,w,be,lam,delta,dom,Gam,rho,maxiter,tol,verb)

  ep = 1e-16;

  gam1 = Gam(1);
  gam2 = Gam(2);

  II = size(T);
  NN = prod(II);
  N  = length(II);

  D  = @(z) z([2:end, end-1],:) - z;
  Dt = @(z) [-z(1,:); z(1:end-3,:) - z(2:end-2,:); z(end-2,:) - z(end-1,:) + z(end,:); z(end-1,:)-z(end,:)];

  vect = @(X) X(:);
  vect_cell = @(Z) vect(cat(ndims(Z{1})+1,Z{:}));
  Dual_para = @(A,B,C) [vect(A); vect(B); vect_cell(C)];
  csum = @(Z) sum(cat(ndims(Z{1})+1,Z{:}),ndims(Z{1})+1);
  inp    = @(Z1,Z2) Z1(:)'*Z2(:);

  % box constraint for gamma
  Pg = @(gam) min(1e3,max(1e-5,gam));

  % initialize
  X = T;
  U = X;
  Y = zeros(NN,N);
  Yp = zeros(NN,N);
  sumLtY = zeros(II);
  sumZ   = zeros(II);
  nu_norm= 0;
  for n = 1:N
    Y(:,n) = w(n)*vect( fold(D(unfold(X,n)),n,II) );
    sumLtY = sumLtY + w(n)*fold(Dt(unfold(reshape(Y(:,n),II),n)),n,II);
    Z{n} = X;
    sumZ = sumZ + Z{n};
    nu_norm = nu_norm + lam(n)*sum(sqrt(svd(unfold(X,n)*unfold(X,n)')));
  end

  obj = (al*sum(sqrt(sum(Y.^2,2)),1) + be*nu_norm)/prod(II);
  if verb
    fprintf('0 :: %f \n',obj);%fflush(1);
  end

  rho_t=rho;
  c   = 0.9;
  eta_p= 1.01;
  eta_d= 1.01;

  for iter = 1:maxiter

    % update X
    Xb = X;
    Xp = X;
    X  = prox_indicator_l1(X - gam1*(U + sumLtY + csum(Z)),T,Q,delta);
    dX = 2*X - Xp;

    % update U
    Ub= U;
    Up= U + gam2*dX;
    U = Up - gam2*prox_dom(Up/gam2,dom);

    % update Y
    Yb = Y;
    Yp(:) = 0;
    for n = 1:N
      Yp(:,n) = Y(:,n) + gam2*w(n)*vect( fold(D(unfold(dX,n)),n,II) );
    end
    Y = Yp - gam2*prox_l12norm(Yp/gam2,al/gam2);

    % update Z
    Zb = Z;
    for n = 1:N
      Zp = unfold(Z{n} + gam2*dX,n);
      Z{n} = fold(Zp - gam2*prox_nuc(Zp/gam2,be*lam(n)/gam2),n,II);
    end

    % calc_obj
    %if skip_obj == 0
      sumLtY(:) = 0;
      V      = zeros(II);
      nu_norm= 0;
      for n = 1:N
        sumLtY = sumLtY + w(n)*fold(Dt(unfold(reshape(Y(:,n),II),n)),n,II);
        V = V + (w(n)^2)*fold(D(unfold(X,n)),n,II).^2;
        nu_norm = nu_norm + lam(n)*sum(sqrt(svd(unfold(X,n)*unfold(X,n)')));
      end
      obj2 = (al*sum(sqrt(V(:)),1) + be*nu_norm)/prod(II);
      histo_obj(iter) = obj2;
    %end

    % check and adjust step size parameters
    Xdif = Xb - X;
    Udif = Ub - U;
    Ydif = Yb - Y;
    Pk = Xdif/gam1 - Udif;
    Dk_1 = Udif/gam2 - Xdif;
    for n = 1:N
      Pk = Pk - w(n)*fold(Dt(unfold(reshape(Ydif(:,n),II),n)),n,II);
      Pk = Pk - (Zb{n} - Z{n});
      Dk_2(:,n) = Ydif(:,n)/gam2 - w(n)*vect( fold(D(unfold(Xdif,n)),n,II) );
      Dk_3{n} = (Zb{n} - Z{n})/gam2 - Xdif;
    end
    Dk = Dual_para(Dk_1, Dk_2, Dk_3);
    Vdif = Dual_para(Ub,Yb,Zb) - Dual_para(U,Y,Z);

    PPk(iter)   = norm(Pk(:));
    DDk(iter)   = norm(Dk(:));
    pdr(iter)   = PPk(iter)/DDk(iter);
    GAM1(iter)  = gam1;
    GAM2(iter)  = gam2;
    condition   = (PPk(iter)^2 + DDk(iter)^2)/NN;
    histo(iter) = condition;

    if mod(iter,verb) == 0
      fprintf('%d :: condintion(%e) :: pdr(%e) :: gam(%e,%e), rho(%e) \n',iter,condition,pdr(iter),gam1,gam2,rho_t);
    end
    if condition < tol
      break;
    end

    Xdif = Xdif/norm(Xdif(:));
    Vdif = Vdif/norm(Vdif(:));
    Pk = Pk/PPk(iter);
    Dk = Dk/DDk(iter);
    Cp(iter) = inp(Xdif,Pk);
    Cd(iter) = inp(Vdif,Dk);

    if Cp(iter) < 0
      gam1 = gam1*0.9;
      eta_p= eta_p^0.1;
    elseif Cp(iter) > 0.9
      gam1 = gam1*eta_p;
    else
    end
    if  Cd(iter) < 0
      gam2 = gam2*0.9;
      eta_d= eta_d^0.1;
    elseif Cd(iter) > 0.9
      gam2 = gam2*eta_d;
    else
    end

    rnd = 1.0;
    gam1 = Pg((gam1)*( pdr(iter)^(rho_t*rnd) ));
    gam2 = Pg((gam2)/( pdr(iter)^(rho_t*rnd) ));

  end


















