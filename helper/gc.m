function [xn, iter, err] = gc (A, b, x0, nmax, tol)
%
%
% Metodo del gradiente coniugato
%
% Parametri di ingresso:
%   A      Matrice del sistema
%   b      Termine noto
%   x0     Vettore iniziale
%   nmax   Numero massimo di iterazioni
%   tol    Tolleranza sul test d'arresto
%
% Parametri in uscita
%   xn     Vettore soluzione
%   iter   Iterazioni effettuate
%   err    Vettore contenente gli errori relativi sul residuo

[n, m] = size (A);
if n~=m
  error ('matrice non quadrata')
end

iter = 0;
xn = zeros(n,1);

% Iterazioni
bnrm2 = norm (b);
r = b - A * x0;
d = r;
errk = norm (r) / bnrm2;
err = errk;
xv = x0;

while (errk > tol) && (iter < nmax)
    alpha = d'*r / (d'*A*d);
    xn = xv + alpha * d;
    r = r - alpha*A*d;
    beta = (A*d)'*r / ((A*d)'*d);
    d = r - beta*d;
    errk = norm (r) / bnrm2;
    err = [err errk];
    xv = xn;
    iter = iter + 1;
end

if (iter == nmax)
    fprintf('Il metodo gc non converge in %d iterazioni \n', iter);
end