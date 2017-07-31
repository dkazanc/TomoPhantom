function [F] = buildSino(ModelNo,G,P,angles)
% this function provides analytical sinograms to phantoms
% (parallel beam geometry)

% input includes:
% 1. model number
% 2. the phantom (N x N)
% 3. detector dimension
% 4. projection angles

name = 'PhantomLibrary.dat';

% read file with parameters
if (exist(name,'file') == 2)
    
    fid = fopen(name, 'r');
    
    if fid < 0
        error(['Could not open ',name,' for input']);
    else
        
        while feof(fid) == 0
            tline = fgetl(fid);
            k = strfind(tline,'#');
            if (isempty(k) == 1)
                matches = strfind(tline, 'Model No.');
                if (matches > 0)
                    [token] = strtok(tline, 'Model No.');
                    tt = textscan(token, '%d %d'); %where - ModNo (model number) and CompNo (components number)
                    ModNo = tt{1};
                    if (ModelNo == ModNo)
                        CompNo = tt{2};
                        C0 = zeros(CompNo,1);
                        IN = zeros(CompNo,1);
                        x0 = zeros(CompNo,1);
                        y0 = zeros(CompNo,1);
                        a = zeros(CompNo,1);
                        b = zeros(CompNo,1);
                        phi_rot = zeros(CompNo,1);
                        for i=1:CompNo
                            tline1 = fgetl(fid);
                            IN(i)=sscanf(tline1, '%d');  %tomo object
                            tline1 = fgetl(fid);
                            ParaM = sscanf(tline1, '%f');
                            C0(i) = ParaM(1);
                            x0(i) = ParaM(2); %x0(i) = x0(i) * H_x; %real step
                            y0(i) = ParaM(3); %y0(i) = y0(i) * H_y; %real step
                            a(i) = ParaM(4);
                            b(i) = ParaM(5);
                            phi_rot(i) = ParaM(6);
                        end % for
                    end % working_Model_no==ModNo
                end % matches == 1
            end % k~=1
        end
        
    end %if error
    fclose(fid);
else
    error('File PhantomLibrary.dat is NOT found');
end

N = size(G,1);
Sinorange_Pmin = -(P)/(N);
Sinorange_Pmax = (P)/(N);
Sinorange_P_Ar = linspace(Sinorange_Pmax, Sinorange_Pmin, P);
Tomorange_Xmin = -1;
Tomorange_Xmax = 1;
H_x = (Tomorange_Xmax - Tomorange_Xmin)/(N); % step for X

C1 = -4 *log(2);
AnglesTot = length(angles);
AnglesRad = angles*(pi/180);

F = zeros(P,AnglesTot);
for i = 1:CompNo   %number of models loop
    a1 = a(i);
    b1 = b(i);
    a22=(a1)^2;
    b22=(b1)^2;
    C00 = C0(i);
    mod = IN(i);
    phi_rot_radian = (phi_rot(i))*(pi/180);
    x00 = x0(i) + H_x;
    y00 = y0(i) + H_x;
    
    if (mod == 1)
        % gaussian
        AA5 = ((C00*a1*b1)/2.0)*sqrt((pi)/log(2));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta1 = 1.0/(a22*cos_2+b22*sin_2);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                under_exp = (C1*AA3)*delta1;
                F(j,ll) = F(j,ll) + first_dr*exp(under_exp);  % sinogramm computing
            end
        end
    elseif (mod == 2)
        % the object is a parabola Lambda = 1/2
        AA5 = ((pi/2.0)*C00*a(i)*b(i));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta1 = 1.0/(a22*cos_2+b22*sin_2);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                if (AA6 < 1)
                    F(j,ll) = F(j,ll) + first_dr*(1. - AA6);  % sinogramm computing
                end
            end
        end
    elseif (mod == 3)
        % the object is an elliptical disk
        AA5 = ((pi/2.0)*C00*a(i)*b(i));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta1 = 1.0/(a22*cos_2+b22*sin_2);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                if (AA6 < 1)
                    F(j,ll) = F(j,ll) + first_dr*(1. - AA6)^(1/2);  % sinogramm computing
                end
            end
        end
    elseif (mod == 12)
        % the object is a parabola Lambda = 1
        a1 = a(i);
        b1 = b(i);
        AA5 = 4.*((0.25*a1*b1*C00)/2.5);
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta = (0.25*a22)*cos_2 + (0.25*b22)*sin_2;
            delta1 = 1/(delta);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                if AA6 < 1
                    F(j,ll) = F(j,ll) + first_dr*(1 - AA6);  % sinogramm computing
                end
            end
        end
    elseif (mod == 13)
        % the object is a cone
        AA5 = C00*a(i)*b(i);
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta = (a22)*cos_2 + (b22)*sin_2;
            delta1 = 1/(delta);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                % might be a bug
                F(j,ll) = F(j,ll) + first_dr*(sqrt((1 - AA6)) - AA6*(1/2)*log2((1 + sqrt((1 - AA6)))./(1 - sqrt(1 - AA6))));  % sinogramm computing
            end
        end
    elseif (mod == 14)
        % the object is a parabola Lambda = 3/2
        AA5 = ((3./8.)*pi*C00*a(i)*b(i));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta = (a22)*cos_2 + (b22)*sin_2;
            delta1 = 1/(delta);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x00*cos(AnglesRad(ll))+y00*sin(AnglesRad(ll)); %p0
            for j = 1:P
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                if AA6 < 1
                    F(j,ll) = F(j,ll) + first_dr*(1 - AA6)^2;  % sinogramm computing
                end
            end
        end
        
    elseif (mod == 18)
        if (phi_rot_radian < 0)
            phi_rot_radian = pi + phi_rot_radian;
        end
        % a rectangle
        for ll = 1:AnglesTot
            for j = 1:P
                F(j,ll) = F(j,ll) + Rectang2DSino(Sinorange_P_Ar(j),AnglesRad(AnglesTot-ll+1),C00,-2*y0(i)-H_x,2*x0(i)+H_x,phi_rot_radian,b(i),a(i));
            end  
        end 
    end
end

return