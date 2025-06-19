-- In this file we verify that for each of 10 combinatorial types of 3-dimensional polytopes
-- with simple facet hyperplane arrangement there is a representative with a smooth adjoint
R = QQ[x,y,z,w]
-- Quartic adjoints
-- Type 7
p1 = 525600*x^4+284760*x^3*y-428400*x^2*y^2-284760*x*y^3-97200*y^4+408300*x^3*z+63372*x^2*y*z-130920*x*y^2*z-45552*y^3*z-138454*x^2*z^2+46836*x*y*z^2+38206*y^2*z^2-76636*x*z^3-18956*y*z^3-51000*x^3*w+40920*x^2*y*w+204000*x*y^2*w+79320*y^3*w-699760*x^2*z*w-84756*x*y*z*w+51046*y^2*z*w-282469*x*z^2*w-158619*y*z^2*w+38318*z^3*w-414760*x^2*w^2-166110*x*y*w^2-18260*y^2*w^2-166747*x*z*w^2-171187*y*z*w^2+175848*z^2*w^2+12750*x*w^3-37230*y*w^3+207276*z*w^3+70840*w^4
I = ideal(diff(x, p1), diff(y,p1), diff(z,p1), diff(w,p1))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

--Type 10
p2 = 4524*x^4+2080*x^3*y-49592*x^2*y^2+6880*x*y^3+3084*y^4+360*x^3*z+31376*x^2*y*z-18280*x*y^2*z-2704*y^3*z-4587*x^2*z^2+9410*x*y*z^2+969*y^2*z^2-1344*x*z^3-152*y*z^3+21*z^4-3888*x^3*w-560*x^2*y*w+112*x*y^2*w-1040*y^3*w+7036*x^2*z*w-856*x*y*z*w-7284*y^2*z*w+1888*x*z^2*w+4716*y*z^2*w-734*z^3*w+8696*x^2*w^2-2240*x*y*w^2+8088*y^2*w^2+2616*x*z*w^2-5408*y*z*w^2+1819*z^2*w^2+944*x*w^3+400*y*w^3+196*z*w^3-1572*w^4
I = ideal(diff(x, p2), diff(y,p2), diff(z,p2), diff(w,p2))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

--Type 14 
p3 = 1368*x^4-2400*x^3*y-16464*x^2*y^2+8800*x*y^3-1032*y^4-780*x^3*z-3432*x^2*y*z+5740*x*y^2*z-2072*y^3*z-1803*x^2*z^2+2908*x*y*z^2+5983*y^2*z^2+681*x*z^3+1285*y*z^3+336*z^4-80*x^3*w+176*x^2*y*w-240*x*y^2*w-624*y^3*w-1376*x^2*z*w+712*x*y*z*w+10048*y^2*z*w+232*x*z^2*w+2386*y*z^2*w+1093*z^3*w+2592*x^2*w^2-1600*x*y*w^2+3968*y^2*w^2-580*x*z*w^2+1256*y*z*w^2-213*z^2*w^2+80*x*w^3+112*y*w^3-1488*z*w^3-632*w^4
I = ideal(diff(x, p3), diff(y,p3), diff(z,p3), diff(w,p3))
dim I 
degree I 
decompose I
radical(I) == ideal(x,y,z,w)

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


