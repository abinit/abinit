clear

pol_eq = [  -0.151086089E+01  -0.151086089E+01  -0.151086089E+01];
pol_pos = [ -0.151780475E+01  -0.151844972E+01  -0.151844972E+01];
pol_neg = [ -0.150398574E+01  -0.150334077E+01  -0.150334077E+01];

piezo_imp = (pol_pos - pol_neg)/0.02;
disp(piezo_imp)

piezo_prop = [0.0 0.0 0.0];
piezo_prop(1) = piezo_imp(1)
piezo_prop(2) = piezo_imp(2) - pol_eq(3)/2
piezo_prop(3) = piezo_imp(3) - pol_eq(3)/2
disp(piezo_prop)

quit
