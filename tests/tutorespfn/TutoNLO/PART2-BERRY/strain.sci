clear
format(20)

// Lattice parameters of equilibrum structure
// ******************************************
acell = [ 7.5232751513E+00  7.5232751513E+00  7.5232751513E+00];
rprim = [     0.0000000000E+00  7.0710678119E-01  7.0710678119E-01
              7.0710678119E-01  0.0000000000E+00  7.0710678119E-01
              7.0710678119E-01  7.0710678119E-01  0.0000000000E+00];
          
          

rprimd(1,1:3) = rprim(1,1:3)*acell(1);
rprimd(2,1:3) = rprim(2,1:3)*acell(2);
rprimd(3,1:3) = rprim(3,1:3)*acell(3);

// Define strain tensor
// ********************
e1 = 0;
e2 = 0;
e3 = 0;
e4 = -0.01;
e5 = 0;
e6 = 0;

e = [ e1    e6/2    e5/2 
      e6/2  e2      e4/2
      e5/2  e4/2    e3 ];

rd(1:3,1:3) = 0;
for i = 1:3
for j = 1:3
 for k = 1:3
  rd(i,j) = rd(i,j) + rprimd(i,k)*e(k,j);
 end
end
end
rd = rd + rprimd

for i = 1:3
 ac(i) = sqrt(rd(i,1)**2 + rd(i,2)**2 + rd(i,3)**2);
 rp(i,1:3) = rd(i,1:3)/ac(i);
end
rp
ac
 
//c0 = rprimd(1,:)+rprimd(2,:)+rprimd(3,:)
//c = rd(1,:) + rd(2,:) + rd(3,:)
//(c(3) - c0(3))/c0(3)

