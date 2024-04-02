n = 32;
M = multipleidentityfactory(n);
problem.M = M;
problem.cost = @(X) n * abs(X - 1).^2;
problem.grad = @(X) 2 * (X - 1);

checkgradient(problem)