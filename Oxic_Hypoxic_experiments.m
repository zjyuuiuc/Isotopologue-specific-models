function dydt = NFmodel(t,y,theta) % theta(1)=net mineralization rate; theta(2)=alpha_min; 
                                   % theta(3)=nitrification rate; 
                                   % theta(4)=alpha_amo;
                                   % theta(5)=nitrate consumption
                                   % rate;theta(6)=alpha_nim
                                   
%global alpha_nim

R15_N2 = 0.0036765;

d15N_orgN = 5.30;
F15_OrgN = (d15N_orgN+1000)/(d15N_orgN+1000+(1000/R15_N2));
F14_OrgN = 1 - F15_OrgN;

F16O_H2O = 0.9976205751661710; %d18O-H2O = 0;
F17O_H2O = 0.0003789960565056;
F18O_H2O = 0.0020004287773232;

F16O_O2 = 0.997569085053419; %d18O-O2 = 23.5
F17O_O2 = 0.000383581767291935;
F18O_O2 = 0.00204733317928882;

alpha_O2 = 1.019;
alpha_H2O1 = 1.018;
alpha_H2O2 = 1.018;

O_ex = 0.2;

alpha_ex = 1.0135;

dydt = zeros(9,1);

dydt(1) = theta(1)*(y(8)/(y(8)+y(9))) - theta(3)*(y(1)/(y(1)+y(2))); % [14N-NH4]

dydt(2) = theta(1)/theta(2)*(y(9)/(y(8)+y(9))) - theta(3)/theta(4)*(y(2)/(y(1)+y(2))); % [15N-NH4]

dydt(3) = theta(3)*(y(1)/(y(1)+y(2))) - theta(5)*(y(3)/(y(3)+y(4))); %[14N-NO3]

dydt(4) = theta(3)/theta(4)*(y(2)/(y(1)+y(2))) - theta(5)/theta(6)*(y(4)/(y(3)+y(4))); %[15N-NO3]

dydt(5) = 2*theta(3)*((0.5*F16O_O2 + 0.5*F16O_H2O)*(1 - O_ex) + O_ex*F16O_H2O) + theta(3)*F16O_H2O - 3*theta(5)*(y(5)/(y(5)+y(6)+y(7)));% [16O-NO3]

dydt(6) = 2*theta(3)*((0.5*F17O_O2/(alpha_O2^0.52) + 0.5*F17O_H2O/(alpha_H2O1^0.52))*(1 - O_ex) + O_ex*F17O_H2O/(alpha_ex^0.52)) + theta(3)*F17O_H2O/(alpha_H2O2^0.52) - 3*theta(5)/(theta(6))^0.52*(y(6)/(y(5)+y(6)+y(7))); % [17O-NO3]

dydt(7) = 2*theta(3)*((0.5*F18O_O2/alpha_O2 + 0.5*F18O_H2O/alpha_H2O1)*(1 - O_ex) + O_ex*F18O_H2O/alpha_ex) + theta(3)*F18O_H2O/alpha_H2O2 - 3*theta(5)/theta(6)*(y(7)/(y(5)+y(6)+y(7))); % [18O-NO3]

dydt(8) = - theta(1)*(y(8)/(y(8)+y(9)));

dydt(9) = - theta(1)/theta(2)*(y(9)/(y(8)+y(9)));
