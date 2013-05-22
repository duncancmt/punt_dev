
L[n_] := gamma * Exp[(64/9)^(1/3) * (n * Log[2])^(1/3) * (Log[n * Log[2]])^(2/3)]
delta = (2^j - 1)^(-1) * M^(-1) * epsilon
expr = L[n] / (36*n * Log[n,2] * delta^(-2)) - 2^(2*j + 9) * n * delta^(-4)
lagrange = lambda * (expr - 2^security) + n^2/j
Assuming[n \[Element] Integers && j \[Element] Integers && security \[Element] Integers && n > 0 && j > 0 && security > 0 && epsilon > 0 && gamma > 0, Solve[D[lagrange, n] == 0 && D[lagrange, j] == 0 && D[lagrange, lambda] == 0, {n, j}, Reals]]
