function [G_deform] = DeformObject(G, RFP, angles_rad, DeformType)

N = size(G,1);
Tomorange_Xmin = -1;
Tomorange_Xmax = 1;
H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N); % step for X
Tomorange_X_Ar = linspace(Tomorange_Xmin, Tomorange_Xmax-H_x, N);
Tomorange_Y_Ar = Tomorange_X_Ar;



G_deform = zeros(N, N);
for i = 1:N
    for j = 1:N
        
           xx = Tomorange_X_Ar(i)*cos(angles_rad) + Tomorange_Y_Ar(j)*sin(angles_rad);
           yy = -Tomorange_X_Ar(i)*sin(angles_rad) + Tomorange_Y_Ar(j)*cos(angles_rad);
        
        if (DeformType == 0)
            % do forward transform
            xPersp1 = xx.*(1 - yy.*RFP);
            yPersp1 = yy;
            xPersp = xPersp1.*(1 - yPersp1.*RFP).*cos(angles_rad) - yPersp1.*sin(angles_rad);
            yPersp = xPersp1.*(1 - yPersp1.*RFP).*sin(angles_rad) + yPersp1.*cos(angles_rad);
        elseif (DeformType == 1)
            % do inverse transform
            xPersp1 = xx./(1 - yy.*RFP);
            yPersp1 = yy;
            xPersp = xPersp1./(1 - yPersp1.*RFP).*cos(angles_rad) - yPersp1.*sin(angles_rad);
            yPersp = xPersp1./(1 - yPersp1.*RFP).*sin(angles_rad) + yPersp1.*cos(angles_rad);
        end        
        
        G_deform(j,i) = BilinearInterpolation(G, Tomorange_Y_Ar, Tomorange_X_Ar, yPersp, xPersp);
    end
end

return
