%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Random Numbers with Laplace distribution %  
%             with MATLAB Implementation                 %
%                                                        %
% Author: M.Sc. Eng. Hristo Zhivomirov          05/01/15 %
% Modified by T. Yokota                         17/10/19 %
%  (m, n) -- > II                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = randl(II)

% function: x  = randl(II)
%%% m - number of matrix rows
%%% n - number of matrix columns
% II - size of x
% x - tensor with Laplacian distributed numbers with mu = 0 and sigma = 1

% generation of a numbers with Uniform distribution
u1 = rand(II);
u2 = rand(II);

% generation of a numbers with Laplace distribution
x = log(u1./u2);

end

