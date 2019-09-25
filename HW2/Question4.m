%% Question 4
%% Part A Obtain State space
% I took the solved dyanmics from 1.7.4 of Zaks textbook, but with a point
% mass and no friction adjustments

syms m M I l theta dtheta ddtheta ddx g fc u a

%a = 1/(m+M)

diffdiffx = -m*a*l*ddtheta*cos(theta)+m*a*l*dtheta^2*sin(theta)-a*fc+a*u 

eq = I*ddtheta == m*g*l*sin(theta)-m*l^2*ddtheta-m*diffdiffx*l*cos(theta)

diffdifftheta = solve(eq,ddtheta)
diffdifftheta = simplify(subs(diffdifftheta, I, m*l^2))
diffdifftheta = simplify(subs(diffdifftheta, fc, 0))

diffdiffx = subs(diffdiffx, ddtheta, diffdifftheta)
diffdiffx = simplify(subs(diffdiffx, I, m*l^2))
diffdiffx = simplify(subs(diffdiffx, fc, 0))

syms x1 x2 x3 x4

output = x1

dx1 = x2
dx2 = diffdiffx
dx2 = subs(dx2, [theta dtheta], [x1 x2])
pretty(dx2)
dx3 = x4
dx4 = diffdifftheta
dx4 = subs(dx4, [theta dtheta], [x1 x2])

g1 = 0
f1 = dx1

g2 = coeffs(dx2, u);
g2 = g2(2)
f2 = simplify(dx2-g2*u)

g3 = 0
f3 = dx3

g4 = coeffs(dx4, u);
g4 = g4(2)
f4 = simplify(dx4-g4*u)

%% Part B

syms x1e x3e u1e u3e

dx2 = subs(dx2, a, (1/(m+M)))
dx4 = subs(dx4, a, (1/(m+M)))

dx2 = subs(dx2, [l m M g], [1 0.1 1 10])
dx4 = subs(dx4, [l m M g], [1 0.1 1 10])

dxAll = [dx1;dx2;dx3;dx4]

A = simplify([diff(dxAll,x1), diff(dxAll,x2),diff(dxAll,x3),diff(dxAll,x4)])

B = simplify(diff(dxAll,u))







