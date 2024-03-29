%% Question 3

%%Setup
syms k1 k2 lambda
A = [0 1; 1 0];
B = [0;1];

%% Make M Matrix
% M matrix is based on linear state feedback control
M = A - B*[k1 k2]

%% Make Characteristic Equation and Subsitute estimator poles
charEq = det(M - lambda*[1 0;0 1]) == 0
eq1 = subs(charEq, lambda, -2 )
eq2 = subs(charEq, lambda, -3 )


%% Solve for Ks
eqs = [eq1 eq2];
ks = [k1 k2]
[solvK1, solvK2] = solve(eqs, ks)

K = [solvK1; solvK2]

% The controller is a linear state feedback controller
% The output state is used to affect input 
% input to the system becomes U = -K*x

