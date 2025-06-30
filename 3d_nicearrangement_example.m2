restart

S = QQ[x0,x1,x2,x3];

-- Example 7.9 (combinatorial type 10 from Figure 5 in Britton-Dunitz)
F0 = x1 - 3*x2 + 2*x3 + x0;         -- adjacent to 1,2,3
F1 = x1 + 3*x2 + 1/5*x3 + x0;       -- adjacent to 0,2,3,5,6
F2 = 3*x1 - x2 + x0;                -- adjacent to 0,1,3,4,6
F3 = -x1 - 3*x2 + 3/2*x3 + x0;      -- adjacent to 0,1,2,4,5
F4 = -3*x1 - x2 - 2/3*x3 + x0;      -- adjacent to 2,3,5,6,7
F5 = -x1 + 3*x2 - 1/2*x3 + x0;      -- adjacent to 1,3,4,6,7
F6 = 3*x1 + x2 - x3 + x0;           -- adjacent to 1,2,4,5,7
F7 = -3*x1 + x2 - 3/2*x3 + x0;      -- adjacent to 4,5,6

-- residual arrangement
R04 = ideal(F0,F4);
R05 = ideal(F0,F5);
R07 = ideal(F0,F7);
R17 = ideal(F1,F7);
R27 = ideal(F2,F7);
R14 = ideal(F1,F4);
R25 = ideal(F2,F5);
R36 = ideal(F3,F6);
R06 = ideal(F0,F6);
R37 = ideal(F3,F7);
 -- residual points
P123 = ideal(F1,F2,F3);
P456 = ideal(F4,F5,F6);

-- adjoint hypersurface
a = (mingens intersect(P123,P456,R04,R05,R06,R07,R17,R27,R37,R25,R36,R14))_(0,0);

------------------------
-- computation of a determinantal representation
-- D ... vanishing ideal of nice arrangement
D = intersect{R04, R05, R07, R17, R27, R14};

Ma = D*(S^1/(ideal(a)*S^1));
rso = (res Ma).dd;
DetRep = rso_1

print(ideal(det(DetRep)) == ideal(a));
-- returns true
