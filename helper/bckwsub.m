function x = bckwsub(A,b)

% Algoritmo di sostituzione all'indietro
% A: matrice quadrata triangolare superiore
% b: termine noto
% x: soluzione del sistema Ax=b

% La matrice deve essere quadrata ed essere compatibile con
% le dimensioni del termine noto
n = length(b);

if ((size(A,1) ~= n) || (size(A,2) ~= n))
  error('Dimensioni incompatibili')
end

% La matrice deve essere inoltre triangolare superiore
if (A == triu(A))

    % La matrice deve essere non singolare.
	% Essendo triangolare, i suoi autovalori si trovano sulla diagonale principale
    if (prod(diag(A)) == 0)  % almeno un elemento diagonale nullo
        error('Matrice singolare')
    end

    x = zeros(n,1);
    x(n) = b(n) / A(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
else
    error('Matrice non triangolare superiore')
end
