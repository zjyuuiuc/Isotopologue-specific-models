function dydt = DNFmodel(t,y,theta)

global alpha_nar                            
global alpha_nxr
global Rnar;global Rnir;global Rnor;global knxr

alpha_MDF = 0.52;

F16O_H2O = 0.997642503169431; %d18O-H2O = -10;
F17O_H2O = 0.000377028810686926;
F18O_H2O = 0.00198046801988179;

F16O_NO2 = 0.9976129028294750; %d18O-NO2 = 3.5;
F17O_NO2 = 0.0003805908382095;
F18O_NO2 = 0.0020122182000000;

alpha_H2O2 = 1.018;

traj = 0.7;

alpha_nxr18 = 0.996;

dydt = zeros(7,1);

dydt(1) =  knxr*y(3) - Rnar*(y(1)/(y(1)+y(2))) - ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))*y(1);% [14N-NO3]

dydt(2) =  knxr/alpha_nxr*y(4) - Rnar/alpha_nar*(y(2)/(y(1)+y(2))) - ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))/alpha_nar*y(2); % [15N-NO3]

dydt(3) =  Rnar*(y(1)/(y(1)+y(2))) - Rnir*(y(3)/(y(3)+y(4))) - knxr*y(3) + ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))*y(1); %[14N-NO2]

dydt(4) =  Rnar/alpha_nar*(y(2)/(y(1)+y(2))) - Rnir/theta(1)*(y(4)/(y(3)+y(4))) - knxr/alpha_nxr*y(4) + ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))/alpha_nar*y(2);%[15N-NO2]

dydt(5) =  knxr*(2*(y(3)+y(4))*F16O_NO2) + knxr*y(3)*F16O_H2O - 3*Rnar*(y(5)/(y(5)+y(6)+y(7))) - ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))*y(5); %[16O-NO3]

dydt(6) = knxr/(alpha_nxr18^alpha_MDF)*(2*(y(3)+y(4))*F17O_NO2) + knxr*y(3)/(alpha_H2O2^alpha_MDF)*F17O_H2O - 3*Rnar/(((alpha_nar-1)*traj+1)^alpha_MDF)*(y(6)/(y(5)+y(6)+y(7))) - ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))/(((alpha_nar-1)*traj+1)^alpha_MDF)*y(6); %[17O-NO3]

dydt(7) = knxr/alpha_nxr18*(2*(y(3)+y(4))*F18O_NO2) + knxr*y(3)/alpha_H2O2*F18O_H2O - 3*Rnar/((alpha_nar-1)*traj+1)*(y(7)/(y(5)+y(6)+y(7))) - ((knxr*y(3)+knxr/alpha_nxr*y(4))/(y(1)+y(2)/alpha_nar))/((alpha_nar-1)*traj+1)*y(7); %[18O-NO3]
