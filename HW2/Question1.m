%% Problem 1

%% Define Symbols
syms m1 m2 q1 q2 q3 a b c g
syms dq1 dq2 dq3 ddq1 ddq2 ddq3
syms tau1 force2 force3

%% Define State Vectors
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

%% Put Torque vector into MCG Form
tau = [ddq1*m2*(b + q2)^2;-m2*(b + q2)*dq1^2 + ddq2*m2; m2*(ddq3 - g)]

M11 = simplify((tau(1) - subs(tau(1),ddq(1),0))/ddq(1));
M12 = simplify((tau(1) - subs(tau(1),ddq(2),0))/ddq(2));
M13 = simplify((tau(1) - subs(tau(1),ddq(3),0))/ddq(3));
M21 = simplify((tau(2) - subs(tau(2),ddq(1),0))/ddq(1));
M22 = simplify((tau(2) - subs(tau(2),ddq(2),0))/ddq(2));
M23 = simplify((tau(2) - subs(tau(2),ddq(3),0))/ddq(3));
M31 = simplify((tau(3) - subs(tau(3),ddq(1),0))/ddq(1));
M32 = simplify((tau(3) - subs(tau(3),ddq(2),0))/ddq(2));
M33 = simplify((tau(3) - subs(tau(3),ddq(3),0))/ddq(3));

M = [M11 M12 M13;
 M21 M22 M23;
 M31 M23 M33];
M = simplify(expand(M))

G = subs(tau, {ddq(1),ddq(2), ddq(3),dq(1),dq(2), dq(3)},{0,0,0,0,0,0})

C1 = simplify(expand(tau(1) - M(1,:)*[ddq1 ddq2 ddq3].' - G(1)));
C2 = simplify(expand(tau(2) - M(2,:)*[ddq1 ddq2 ddq3].' - G(2)));
C3 = simplify(expand(tau(3) - M(3,:)*[ddq1 ddq2 ddq3].' - G(3)));
C = [C1;C2;C3]

%% Torque and Force equations
torqEq1 = M(1,1)*ddq1+C(1)*dq1+G(1) == tau1;
forceEq2 = M(2,2)*ddq2+C(2)*dq2+G(2) == force2;
forceEq3 = M(3,3)*ddq3+C(3)*dq3+G(3) ==force3;

%% Change to State Space Format
% establish state vector
x1 = q1; x2 = dq1;
x3 = q2; x4 = dq2;
x5 = q3; x6 = dq3;

% find state space
dx1 = x2;
dx2 = solve(torqEq1,ddq1);
dx3 = x4;
dx4 = solve(forceEq2,ddq2);
dx5 = x6;
dx6 = solve(forceEq3,ddq3);

stateSpace = [dx1;dx2;dx3;dx4;dx5;dx6]