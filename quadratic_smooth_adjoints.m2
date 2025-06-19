-- This file contains code used to compute the quadratic adjoint of a 4D polytope with no planes in the residual arrangement and to show that it is smooth
-- The same is done for two combinatorial types in 5 dimensions
R = QQ[x0,x1,x2,x3,x4]

-- The polytope is given by the inequalities x_i>=0 for i=0,...,4 and by h1>=0, h2>=0 with h1, h2 below
h1 = 2*x0 + 3*x1 + x2 - x3 - 3*x4
h2 = -2*x0 - x1 + x2 + 2*x3 + 2*x4

-- The residual arrangement consists of 7 lines. This can be verified e.g. in Wolfram Mathematica 
r1 = ideal(x0, x1, x2)
r2 = ideal(x0, x1, h2)
r3 = ideal(x0, h1, h2)
r4 = ideal(x1, x2, x3)
r5 = ideal(x2, x3, x4)
r6 = ideal(x3, x4, h1)
r7 = ideal(x4, h1, h2)

-- Ideal of the residual arrangement
I = intersect(r1, r2, r3, r4, r5, r6, r7)

-- the adjoint
adj = first first entries gens I 

-- Smoothness test
J = ideal(diff(x0, adj), diff(x1,adj), diff(x2,adj), diff(x3,adj), diff(x4,adj))
radical(J) == ideal(x0,x1,x2,x3,x4)
-- returns true, thus the adjoint is smooth and has no determinantal representation

-- The same for two combinatorial types of polytopes in P^5 with 8 facets having no codimension two subspaces in the residual arrangement 
R = QQ[x0, x1, x2, x3, x4, x5]

--Residual arrangement consists of 5 planes and 1 line
h1 = 3*x0 + 3*x1 + x2 - x3 - 2*x4 - 3*x5
h2 = -5*x0 - 4*x1 - x2 + 3*x3 + 3*x4 + 4*x5

s1 = ideal(x0, x1, x2)
s2 = ideal(x0, h1, h2)
s3 = ideal(x1, x2, x3)
s4 = ideal(x2, x3, x4)
s5 = ideal(x3, x4, x5)
s6 = ideal(x4, x5, h1, h2)
I = intersect(s1, s2, s3, s4, s5, s6)
adj = first first entries gens I 

J = ideal(diff(x0, adj), diff(x1,adj), diff(x2,adj), diff(x3,adj), diff(x4,adj), diff(x5, adj))
radical(J) == ideal(x0,x1,x2,x3,x4,x5)

--Residual arrangement consists of 4 planes and 3 lines
h1 = 4*x0 + 2*x1 + 5*x2 + x3 - x4 - 5*x5 
h2 = -5*x0 - x1 - x2 + x3 + 2*x4 + 4*x5

s1 = ideal(x0, h1, h2)
s2 = ideal(x3, x4, x5)
s3 = ideal(x4, x5, h1)
s4 = ideal(x5, h1, h2)
s5 = ideal(x0, x1, x2, x3)
s6 = ideal(x0, x1, x2, h2)
s7 = ideal(x1, x2, x3, x4)
I = intersect(s1, s2, s3, s4, s5, s6, s7)
adj = first first entries gens I 

J = ideal(diff(x0, adj), diff(x1,adj), diff(x2,adj), diff(x3,adj), diff(x4,adj), diff(x5, adj))
radical(J) == ideal(x0,x1,x2,x3,x4,x5)
