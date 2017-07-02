function [Image_G] = Deform_Forward_IMG_F(G, angle,...
    SizeimageX, SizeimageY, Tomorange_X_Ar, Tomorange_Y_Ar, FocalP)

%Forward deformation of an object from any given angle

angles_rad = (pi/180)*angle;

SizeX = size(G,1);
SizeY = size(G,2);

Image_G = zeros(SizeX, SizeY);

for i = 1:SizeimageX
    for j = 1:SizeimageY
        
        
        xx = Tomorange_X_Ar(i)*cos(angles_rad) + Tomorange_Y_Ar(j)*sin(angles_rad);
        yy = -Tomorange_X_Ar(i)*sin(angles_rad) + Tomorange_Y_Ar(j)*cos(angles_rad);
        
        xPersp1 = xx*(1 - yy/FocalP);
        yPersp1 = yy;
        xPersp = xPersp1*(1 - yPersp1/FocalP)*cos(angles_rad) - yPersp1*sin(angles_rad);
        yPersp = xPersp1*(1 - yPersp1/FocalP)*sin(angles_rad) + yPersp1*cos(angles_rad);
        
        %xPersp = Tomorange_X_Ar(i)*(1 - Tomorange_Y_Ar(j)/FocalP)*cos(angles_rad) - Tomorange_Y_Ar(j)*sin(angles_rad);
        %yPersp = Tomorange_X_Ar(i)*(1 - Tomorange_Y_Ar(j)/FocalP)*sin(angles_rad) + Tomorange_Y_Ar(j)*cos(angles_rad);
        
        Image_G(j,i) = BilinearInterpolation(G, Tomorange_Y_Ar, Tomorange_X_Ar, yPersp, xPersp);
        
    end %for
end %for

return