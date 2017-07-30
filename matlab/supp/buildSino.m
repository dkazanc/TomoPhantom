function [F] = buildSino(ModelNo,G,P,angles)
% this function provides analytical sinograms to phantoms
% it uses parallel beam geometry

% input includes:
% 1. model number
%
% 2. detector dimension
% 3. projection angles

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
Sinorange_Pmin = -P/(N);
Sinorange_Pmax = P/(N);
% H_p = (Sinorange_Pmax - Sinorange_Pmin)/(P); % step for p-detectors
% Sinorange_P_Ar = linspace(Sinorange_Pmin, Sinorange_Pmax-H_p, P);
Sinorange_P_Ar = linspace(Sinorange_Pmin, Sinorange_Pmax, P+1);
% Sinorange_P_Ar = linspace(Sinorange_Pmin,Sinorange_Pmax,P);
C1 = -4 *log(2);
AnglesTot = length(angles);
AnglesRad = angles*(pi/180);

F = zeros(P,AnglesTot);
for i = 1:CompNo   %number of models loop
    a1 = a(i);
    b1 = b(i);
    a22=(a1)^2;
    a2=1/a22;
    b22=(b1)^2;
    b2=1/b22;
    C00 = C0(i);
    mod = IN(i);
    phi_rot_radian = (phi_rot(i))*(pi/180);
    

    if (mod == 1)
        % gaussian
        AA5 = ((C0(i)*a1*b1)/2.0)*sqrt((pi)/log(2));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta1 = 1.0/(a22*cos_2+b22*sin_2);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x0(i)*cos(AnglesRad(ll))+y0(i)*sin(AnglesRad(ll)); %p0
            for j = 1:P-1
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                under_exp = (C1*AA3)*delta1;
                F(P-j,ll) = F(P-j,ll) + first_dr*exp(under_exp);  % sinogramm computing                
            end
        end               
    elseif (mod == 2)
        % the object is a parabola Lambda = 1/2
        AA5 = ((pi/2.0)*C0(i)*a(i)*b(i));
        for ll = 1:AnglesTot
            sin_2 = (sin((AnglesRad(ll)) + phi_rot_radian)).^2;
            cos_2 = (cos((AnglesRad(ll)) + phi_rot_radian)).^2;
            delta1 = 1.0/(a22*cos_2+b22*sin_2);
            delta_sq = sqrt(delta1);
            first_dr = AA5*delta_sq;
            AA2 = -x0(i)*cos(AnglesRad(ll))+y0(i)*sin(AnglesRad(ll)); %p0
            for j = 1:P-1
                AA3 = (Sinorange_P_Ar(j) - AA2)^2; %(p-p0)^2
                AA6 = (AA3)*delta1;
                if (AA6 < 1)
                 F(P-j,ll) = F(P-j,ll) + first_dr*(1. - AA6);  % sinogramm computing
                end                
            end
        end            
    elseif (mod == 3)
        % the object is an elliptical disk
        T = (a2*(X*cos_phi+Y*sin_phi).^2 + b2*(-X*sin_phi+Y*cos_phi).^2);
        T(T <= 1) = C00;
        T(T>1) = 0;
        G = G + T;
    elseif (mod == 12)
        % the object is a parabola Lambda = 1
        a2x=4.*a2;
        b2y=4.*b2;
        T = (a2x*(X*cos_phi+Y*sin_phi).^2 + b2y*(-X*sin_phi+Y*cos_phi).^2);
        T(T <= 1) = C00.*sqrt(1.0 - T(T <= 1));
        T(T>1) = 0;
        G = G + T;
    elseif (mod == 13)
        % the object is a cone
        T = a2*(X*cos_phi+Y*sin_phi).^2 + b2*(-X*sin_phi+Y*cos_phi).^2;
        T(T <= 1) = C00.*(1.0 - sqrt(T(T <= 1)));
        T(T>1) = 0;
        G = G + T;
    elseif (mod == 14)
        % the object is a parabola Lambda = 3/2
        T = (a2*(X*cos_phi+Y*sin_phi).^2 + b2*(-X*sin_phi+Y*cos_phi).^2);
        T(T <= 1) = C00.*((1.0 - T(T <= 1)).^1.5);
        T(T>1) = 0;
        G = G + T;
    elseif (mod == 18)
        % a rectangle
        x0r=x0(i)*cos(0) + y0(i)*sin(0);
        y0r=-x0(i)*sin(0) + y0(i)*cos(0);
        if (phi_rot_radian < 0)
            phi_rot_radian = pi + phi_rot_radian;
            sin_phi=sin(phi_rot_radian);
            cos_phi=cos(phi_rot_radian);
        end
        a2=a(i)*0.5;
        b2=b(i)*0.5;
        
        T = zeros(N,N);
        HX = abs((X-x0r)*cos_phi + (Y-y0r)*sin_phi);
        HY = zeros(N,N);
        HY(HX <= a2) = abs((Y(HX <= a2) - y0r)*cos_phi -(X(HX <= a2) - x0r)*sin_phi);
        T((HY <= b2) & (HX <= a2)) = C00;
        G = G + T;
    end
end

return