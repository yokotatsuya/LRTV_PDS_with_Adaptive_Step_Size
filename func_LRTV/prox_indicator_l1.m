% calculate Proximal Projection
% x = argmin_x || x - a ||_2^2, s.t. || q .* (f - x) ||_1 <= delta
%
% This code was written by Tatsuya Yokota
function x = prox_indicator_l1(a,f,q,delta)

  II = size(a);
  x  = zeros(II);
  x(q==0) = a(q==0);

  v = a(q(:)==1) - f(q(:)==1);
  s = sign(v);
  p = abs(v);
  Psum = sum(p);
  if Psum <= delta
    x(q==1) = a(q==1);
  else
    N = length(p);
    [pii id] = sort(p);
    Pk = Psum - cumsum(diff([0; pii]) .* [N:-1:1]');
    k  = sum(Pk > delta);
    mu = pii(k+1) - (delta - Pk(k+1))/(N-k);
    x(q(:)==1) = f(q(:)==1) + s.*max(p-mu,0);
  end


