function [x,w] = Gaulagwt(n, alf)
% Gaulagwt.m

% This script is for computing definite integrals using Laguerre-Gauss 
% Quadrature. Computes the Laguerre-Gauss nodes and weights  on an interval
% [0, inf] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [0, inf]
% which you can evaluate at any x in [0,inf]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using w*f;

% Input
% n   =  number of nodes for the integration
% alf = the power of x in the generalized form of gauss laguerre integral

% Output
% x   = absicsa at which the function is to be evaluated
% w   = weights corresponding to each abscisca
%
%
% Note: w is a row vector while x is a column vector. This choice is made
% so that if f is vectorized, we can compute the integration by w*f(x)

% Written by Lateef Kareem - 17/10/2018

if nargin == 1
    alf = 0;
end

MAXIT = 20;
p2 = 0;
pp = 0;
z = 0; 
eps = 3e-14;
x = zeros(n, 1); 
w = zeros(1, n);
for i = 1:n
    % Loop over the desired roots.
    if (i == 1)
        % Initial guess for the smallest root.
        z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
    elseif (i == 2)
        %Initial guess for the second root.
        z = z + (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
    else
        % Initial guess for the other roots.
        ai = i - 2;
        z = z + ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf / ....
                (1.0 + 3.5 * ai)) * (z - x(i - 2)) / (1.0 + 0.3 * alf);
    end
    for its = 1: MAXIT
        % Refinement by Newton?s method.
        p1 = 1.0;
        p2 = 0.0;
        for j = 1:n
            % Loop up the recurrence relation to get the
            p3 = p2; % Laguerre polynomial evaluated at z.
            p2 = p1;
            p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
        end
        % p1 is now the desired Laguerre polynomial. We next compute pp, its derivative,
        % by a standard relation involving also p2, the polynomial of one lower order.
        pp = (n * p1 - (n + alf) * p2) / z;
        z1 = z;
        z = z1 - p1 / pp; % Newton?s formula.
        if (abs(z - z1) <= eps); break; end
    end
    x(i) = z; % Store the root and the weight.
    w(i) = -exp(gammaln(alf + n) - gammaln(n)) / (pp * n * p2);
end
