-- This code accompanies the example of the adjoint of a planar polygon 
-- with a determinantal representation having recursive structure
R = QQ[x1,x2]
-- Define the adjoint polynomials
a3 = -70 - 116*x1 + 46*x1^2 - 4*x1^3 - 36*x2 - 92*x1*x2 + 15*x1^2*x2 + 30*x2^2 + 24*x1*x2^2 - 4*x2^3
a2 = -22 - 40*x1 + 8*x1^2 - 16*x2 - 40*x1*x2 + 6*x2^2
a1 = 12*x1-x2+6
-- Define the two edges used in Dixon's process
l1 = x1-2*x2-2 
l2 = x2

-- Define the off-diagonal adntries of the adjugate matrix Madj by applying Max Noether's theorem
l1*a1*l1*a1 // gens ideal(a2, a3)
m22 = 312*x1^2+358*x1*x2-32*x2^2-1728*x1+349*x2-942

l1*l2*l1*l2 // gens ideal (a2, a3)
m33 = -6*x1^2+22*x1*x2+40*x2^2+72*x1-143*x2-210

l1*l2*l1*a1 // gens ideal (a2, a3)
m23 = 26*x1*x2+32*x2^2-157*x2

-- Construct the adjugate matrix
Madj = matrix{{a2, l1*a1, l1*l2},{l1*a1, m22, m23},{l1*l2, m23, m33}}

-- One then needs to take the adjugate of Madj and divide each entry by a3. We did this in Wolfram Mathematica. The resulting determinant is below. 
dr = -54454680 - 90239184*x1 + 35784504*x1^2 - 3111696*x1^3 - 28005264*x2 - 71569008*x1*x2 + 11668860*x1^2*x2 + 23337720*x2^2 + 18670176*x1*x2^2 - 3111696*x2^3

-- Verifying correctness
dr % ideal(a3)
subm2 = matrix{{-66 + 12*x1 - 59*x2, 49*x2}, {49*x2, -49*(6 + 12*x1 - x2)}}
det(subm2) % ideal(a2)