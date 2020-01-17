
function [A]=cholesky_fact(A)
//-----------------------
// factorise en place
//----------------------
n = size(A,"r");     // A est supposee de taille n*n.
A = tril(A);
A(1, 1) = sqrt(A(1, 1));
A(2:n, 1) = (1 / A(1, 1)) * A(2:n, 1);
for column=1:(n-1),
    for line=(column + 1):n,
        A(line:n, line) = A(line:n, line) - A(line, column)*A(line:n, column);
    end
    A(column + 1, column + 1) = sqrt(A(column + 1, column + 1));
    A(column + 2:n, column + 1) = (1 / A(column + 1, column + 1)) * A(column + 2:n, column + 1);
end
endfunction

function [y]=up_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
    y = []
    y(n) = x(n) / A(n, n)
    for line=(n - 1):-1:1,
        y(line) = (x(line) - sum(A(line, (line + 1):n) * y((line + 1):n))) / A(line, line);
    end
  end
endfunction

function [y]=down_sweep_cholesky(A,x)
  [m,n]=size(A);
  if (m~=n) then
    print(%io(2), "error, not a square matrix");
  else
    y = []
    y(1) = x(1) / A(1, 1)
    for line=2:n,
        y(line) = (x(line) - sum(A(line, 1:(line - 1)) * y(1:(line - 1)))) / A(line, line);
    end
end
endfunction

function [U]=my_cholesky(N,S)
//---------------
// On decompose N = LL'.
// On cherche y tel que Ly = S, puis on cherche x tel que L'x = y.
// Dans ce cas, x = U recherche.
//--------------
L = cholesky_fact(N)
y = down_sweep_cholesky(L,S)
U = up_sweep_cholesky(L',y)
endfunction

