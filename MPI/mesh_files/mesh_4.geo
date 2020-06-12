//cylinder

lc1 = 1.0;

Point(1) = {4, -6, 0, lc1};
Point(2) = {30, -6,  0,lc1} ;
Point(3) = {30, 1.2, 0, lc1} ;
Point(4) = {4, 1.2, 0, lc1} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {4,1,2,3} ;

lc2 = lc1/4.0;
Point(5) = {-1, 0, 0, lc2};
Point(6) = {0, 1, 0, lc2};
Point(7) = {1, 0, 0, lc2};
Point(8) = {0, -1, 0, lc2};
Point(9) = {0, 0, 0, lc2};





Plane Surface(11) = {5} ;

inflow = 12;
outflow = 13;
noslip = 14;
inside = 15;



Physical Line(outflow) = {2} ;
Physical Line(noslip) = {3};


Physical Surface(1) = 11; 

//Physical Surface("My fancy surface label") = {11} ;


