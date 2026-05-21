function x = fwsub(A,b)

% Algoritmo di sostituzione in avanti
% A: matrice quadrata triangolare inferiore
% b: termine noto
% x: soluzione del sistema Ax=b

% La matrice deve essere quadrata ed essere compatibile con
% le dimensioni del termine noto
n = length(b);

if ((size(A,1) ~= n) || (size(A,2) ~= n))
  error('Dimensioni incompatibili')
end

% La matrice deve essere inoltre triangolare inferiore
if (A == tril(A))

    % La matrice deve essere non singolare.
	% Essendo triangolare, i suoi autovalori si trovano sulla diagonale principale
    if (prod(diag(A)) == 0)  % almeno un elemento diagonale nullo
        error('Matrice singolare')
    end

    x = zeros(n,1);
    x(1) = b(1) / A(1,1);
    for i = 2:n
        x(i) = (b(i) - A(i,1:i-1) * x(1:i-1)) / A(i,i);
    end
else
    error('Matrice non triangolare inferiore')
end
