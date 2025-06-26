-- In this file we verify that for each of 10 combinatorial types of 3-dimensional polytopes
-- with simple facet hyperplane arrangement there is a representative with a smooth adjoint
R = QQ[x,y,z,w]
-- Quartic adjoints
-- Type 7
A = 3*y+3/2 * z + w 
B = 2*x + 2*y + 2/5*z + w 
C = 3*x-1/5 *z + w 
D = 2*x-2*y-1/2*z+w 
E = -3*y+z+w 
F = -2*x - 2*y-z + w 
G = -3*x-3/5*y-z + w 
H = -2*x+2*y-z+w 
s1 = ideal(A, D)
s2 = ideal(A, F)
s3 = ideal(B, D)
s4 = ideal(B, E)
s5 = ideal(B, F)
s6 = ideal(B, G)
s7 = ideal(C, F)
s8 = ideal(D, H)
s9 = ideal(E, H)
s10 = ideal(G, D, E)
s11 = ideal(A, C, G)
s12 = ideal(A, C, H)
s13 = ideal(C, E, G)
I = intersect(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13)
p1 = first first entries gens I

--Smoothness test
J = ideal(diff(x, p1), diff(y,p1), diff(z,p1), diff(w,p1))
radical(J) == ideal(x,y,z,w)

--Type 10
A = -x + 3*y - 1/2*z + w
B = x + 3*y + 1/5*z + w 
C = 3*x+y-z+w 
D = 3*x-y+w 
E = x-3*y+2*z + w 
F = -x - 3*y + 3/2*z + w 
G  = -3*x - y - 2/3*z + w 
H = -3*x+y-3/2 * z + w 

s1 = ideal(A, D)
s2 = ideal(A, E)
s3 = ideal(B, H)
s4 = ideal(B, G)
s5 = ideal(C, F)
s6 = ideal(D, H)
s7 = ideal(E, G)
s8 = ideal(E, H)
s9 = ideal(F, H)
s10 = ideal(A, C, G)
s11 = ideal(B, D, F)
s12 = ideal(C, E)
I = intersect(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12)
p2 = first first entries gens gb I

J = ideal(diff(x, p2), diff(y,p2), diff(z,p2), diff(w,p2))
radical(J) == ideal(x,y,z,w)

--Type 14 
A = -x + 3*y - z/2 + w 
B = x + 3*y -z/2 + w 
C = 3*x+y+3*z+ w 
D = 3*x-y+z+w 
E = x - 3*y - 3/2*z + w 
F = -x -3*y - z + w 
G = -3*x - y + 3*z + w
H = -3*x + y + 2*z + w
 
s1 = ideal(A, D)
s2 = ideal(A, G)
s3 = ideal(B, F)
s4 = ideal(B, G)
s5 = ideal(B, H)
s6 = ideal(C, E)
s7 = ideal(C, F)
s8 = ideal(D, H)
s9 = ideal(D, F) 
s10 = ideal(E, H)
I = intersect(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10) 
p3 = first first entries gens gb I

J = ideal(diff(x, p3), diff(y,p3), diff(z,p3), diff(w,p3))
radical(J) == ideal(x,y,z,w)

-- Cubic adjoints 
-- Type 3
a = -x+w 
b = -3*x+2*y-z+w 
c = 4*y+w-z/4
d = 5*y+2*z+w 
e = 3*x+4*y+w 
f = 3*x+2*y+w 
g = x+z+w
J = intersect(ideal(a,d,f), ideal(b,e,d), ideal(b,d,f),ideal(b,g), ideal(a,e), ideal(a,c), ideal(c,f), ideal(c,g), ideal(e,g))
p1 = first first entries gens J
I = ideal(diff(x, p1), diff(y,p1), diff(z,p1), diff(w,p1))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

-- Type 4
a = -x+w 
b = -3*x+2*y+z+w 
c = -3*x+4*y+3/2*z+w 
d = 5*y-7/2*z+w 
e = 3*x+4*y-z+w 
f = 3*x+2*y+4*z+w 
g = x+2/3*z+w 
J = intersect(ideal(a,c), ideal(b,e), ideal(b,g), ideal(c,g), ideal(d,f), ideal(d,g), ideal(a,e,f))
p2 = first first entries gens J 
I = ideal(diff(x, p2), diff(y,p2), diff(z,p2), diff(w,p2))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

-- Type 5
a = -x-2*z+w
b = -3*x+2*y+w
c = -3*x + 4*y + z/2 + w 
d = 5*y-z+w
e = 3*x+4*y+w
f = 3*x+2*y+2*z+w 
J = intersect(ideal(a,c), ideal(a,f), ideal(b,e), ideal(c,g), ideal(d,f), ideal(d,g)) 
p3 = first first entries gens J 
I = ideal(diff(x, p3), diff(y,p3), diff(z,p3), diff(w,p3))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

-- Quadratic adjoints 

-- Perturbed cube
J = intersect(ideal(x+y+w,x+y+4*w), ideal(2*x+y/2+z+w,2*y+2*z+w), ideal(z+w,3/2*x+5/2*y+3/2*z+w))
p1 = first first entries gens J
I = ideal(diff(x, p1), diff(y,p1), diff(z,p1), diff(w,p1))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

-- The other combinatorial type, also perturbed to get a simple facet hyperplane arrangement
a = 4*w-3*y+2*z
b =  4*w-y+2*z
c = -4*x-5*y-2*z+20*w
d = 4*x-5*y-2*z+20*w
e  = -26*x-5*y+12*z+50*w
f = 2*w+2*x+y
J = intersect(ideal(c,f), ideal(d,e), ideal(e,f), ideal(a,b,d), ideal(a,b,c))
p2 = first first entries gens J 
I = ideal(diff(x, p2), diff(y,p2), diff(z,p2), diff(w,p2))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)


