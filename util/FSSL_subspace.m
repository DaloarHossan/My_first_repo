function A = FSSL_subspace(X, Y, regu)
[d, ~] = size(X);
[n, nClass] = size(Y);
% Check the solutions
nSolutionCheck = 0;
r1 = rank(X');
r2 = rank([X', Y]);
if r1 == r2 && r1 < d
    % X'*A = Y has many solution == rank(X') == rank([X', Y]) < d
    nSolutionCheck = 1;
end

A = zeros(d, nClass);
% Step 2: Find A satisfies the linear system
nIter = 20;
if nSolutionCheck
    % Situation 1, Infinitely many solutions
    G = eye(d);
    for iter = 1:nIter
        Gi = inv(G);
        A = Gi*X*inv(X'*Gi*X)*Y; %#ok
        normG = sqrt(sum(A.^2,2));
        nzIdx = (normG ~= 0);
        dd = zeros(d, 1);
        dd(nzIdx) = 1./normG;
        G = diag(dd);
    end
else
    % Situation 1, Single or No solution
    G = eye(d);
    for iter = 1:nIter
        Gi = inv(G);
        A = Gi*X*inv(X'*Gi*X + 0.5/regu*eye(n))*Y; %#ok
        normG = sqrt(sum(A.^2,2));
        nzIdx = (normG ~= 0);
        dd = zeros(d, 1);
        dd(nzIdx) = 1./normG;
        G = diag(dd);
    end
end
end

