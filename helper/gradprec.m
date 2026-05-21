function [xn, iter, err] = gradprec (A, b, P, x0, nmax, tol)
%
%        [xn, iter, err] = gradprec (A, b, P, x0, nmax, tol)
% Metodo del gradiente precondizionato
%
% Parametri di ingresso:
%   A      Matrice del sistema
%   b      Termine noto
%   P      Precondizionatore
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
errk = norm (r) / bnrm2;
err = errk;
xv = x0;

while (errk > tol) && (iter < nmax)
    z = P\r;
    alpha = z'*r / (z'*A*z);
    xn = xv + alpha * z;
    r = r-alpha*A*z;
    errk = norm (r) / bnrm2;
    err = [err errk];
    xv = xn;
    iter = iter + 1;
end

if (iter == nmax)
    fprintf('Il metodo gradprec non converge in %d iterazioni \n', iter);
end