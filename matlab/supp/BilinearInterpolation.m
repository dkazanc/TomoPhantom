function [Interp_P3] = BilinearInterpolation(F, AR1, AR2, xx, yy)

%Bilinear 2D Interpolation
SXmax = size(AR1,2);
SYmax = size(AR2,2);
SXel_min = AR1(1);
SYel_min = AR2(1);

Ar1step = AR1(2) - AR1(1);
Ar2step = AR2(2) - AR2(1);

Ar1step_inv = 1/Ar1step;
Ar2step_inv = 1/Ar2step;

ll = (xx-SXel_min)*Ar1step_inv; 
mm = (yy-SYel_min)*Ar2step_inv;

i0 = floor(ll);
j0 = floor(mm);
u = ll - i0;
v = mm - j0;

i=i0+1;
j=j0+1;


i1 = i+1;
j1 = j+1;

%two variant of interpolation code (both are correct!)

%         %Interpolation in first point
%         Interp_P1 = F(i,j) * ((AR1(i1) - xx)*Ar1step_inv) + F(i1,j) * ((xx - AR1(i))*Ar1step_inv);
%         %Interpolation in second point
%         Interp_P2 = F(i,j1) * ((AR1(i1) - xx)*Ar1step_inv) + F(i1,j1) * ((xx - AR1(i))*Ar1step_inv);
%         %Interpolation in third point(between prevoius two)
%         Interp_P3 = Interp_P1*((AR2(j1) - yy)*Ar2step_inv) + Interp_P2*((yy - AR2(j))*Ar2step_inv);
%         %

%Interpolation with border cases (where i = SXmax)
a = 0; b = 0; c = 0; d = 0;

if (i > 0) && (i <= SXmax)
    if (j > 0) && (j <= SYmax)
        a = F(i,j);
    end
end
if (i > 0) && (i < SXmax)
    if (j > 0) && (j <= SYmax)
        b = F(i1,j);
    end
end
if (i > 0) && (i <= SXmax)
    if (j > 0) && (j < SYmax)
        c = F(i,j1);
    end
end
if (i > 0) && (i < SXmax)
    if (j > 0) && (j < SYmax)
        d = F(i1,j1);
    end
end


Interp_P3 = (1 - u)*(1 - v)*a + u*(1 - v)*b + (1 - u)*v*c+...
    u*v*d;

return