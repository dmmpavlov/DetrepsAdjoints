restart
R = QQ[x0,x1,x2];

-- compute the adjoint of a subpolytope given by a subset of V
    -- if the second argument is a positive integer n, then compute the adjoint of V_{0..n-1}
    -- if it is a list of natural numbers, compute the subpolytope with corresponding vertex labels
subadj = (V,subVertices) -> (
    if class subVertices === ZZ then (
    -- subpolygon which vertices V_{0..subVertices-1})
        -- subVertices: integer between 3 and length(V)=:n
    k := subVertices-1;
    Vn = V_{0..k};
    Ln = (toList(0..k))/(i -> (mingens intersect(Vn_{i,(i+1)%(k+1)}))_(0,0));
    ResVn = subsets(0..k,2)/(pair -> ideal(Ln_pair));
    Resn = (intersect ResVn) : (intersect Vn);
    adjn = (mingens Resn)_(0,0);
    return adjn;
    --
    )
    else if class subVertices === List then (
    -- subpolygon with vertices as given in subVertices
        -- subVertices: list of labels of length between 3 and length(V)=:n
    Vn = V_subVertices;
    k = (length subVertices)-1;
    Ln = (toList(0..k))/(i -> (mingens intersect(Vn_{i,(i+1)%(k+1)}))_(0,0));
    ResVn = subsets(0..k,2)/(pair -> ideal(Ln_pair));
    Resn = (intersect ResVn) : (intersect Vn);
    adjn = (mingens Resn)_(0,0);
    return adjn;
    --
    )
    else return print "Expected integer or list";
);

----- exmample 5.13.
-- vertex data of polygon (should be given in cw or ccw order)
V0 = ideal(x1-3*x0, x2-3*x0); -- (3,3)
V1 = ideal(x1-1*x0, x2-7*x0); -- (1,7)
V2 = ideal(x1+2*x0, x2-6*x0); -- (-2,6)
V3 = ideal(x1+3*x0,x2-5*x0);  -- (-3,5)
V4 = ideal(x1+3*x0, x2-3*x0); -- (-3,3)
V5 = ideal(x1+x0, x2-x0);     -- (-1,1)
V6= ideal(x1-x0, x2-x0);      -- (1,1)


V = {V0,V1,V2,V3,V4,V5,V6};
L = toList(0..length(V)-1) / (i -> (mingens intersect(V#((i-1)%length(V)), V#(i%length(V))))_(0,0));


-- computing a DetRep
M = mutableMatrix(R,length(V)-3, length(V)-3)
M_(0,0) = subadj(V,{0,1,2,3});

for i in 0..length(V)-5 do (
    M_(i,i+1) = L#(i+3); M_(i+1,i) = L#(i+3);
    S := R/ideal(subadj(V, {0,i+1,i+2,i+3,i+4}));
    M_(i+1,i+1) = sub(sub(M_(i,i+1)^2, S)//sub(M_(i,i),S),R);
)

M = matrix(M);
print(ideal det M == ideal subadj(V,length(V)));
