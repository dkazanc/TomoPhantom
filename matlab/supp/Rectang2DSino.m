function [SS] = Rectang2DSino(p00,ksi00,C0,x0,y0,ksi1,xwid,ywid)
% !      Calculate line integral of phantom: Rectangle with parameters
% !      c0, x0, y0, ksi1, xwid, ywid; along line (p,ksi).
% !
% !    xwid, ywid - full width of rectangle;
% !
% !               ksi is in rad,
% !               ksi1 also is in rad (is angle of tilt)
% !
% !   Constructed - 14 March 2001
% !    Updated - 14.03.01
% !    number of revisions - 1
PI2 = pi*0.5;
p=p00;
ksi=ksi00;

if (ksi > pi)
    ksi=ksi-pi;
    p=-p00;
end %(ksi > pi)

C =cos(ksi); S=sin(ksi);
XSYC = -x0*S + y0*C;
A2 = xwid*0.5;
B2 = ywid*0.5;

%    .....FI - ���� ������-� � ���. ����.,0<=FI<=PI/2
if (ksi - ksi1) < 0.0
    FI = pi + ksi - ksi1;
else
    FI = ksi - ksi1;
end %(ksi - ksi1 < 0.0)

if (FI > PI2)
    FI = pi - FI;
end %(FI > PI2)

CF = cos(FI);
SF = sin(FI);
P0 = abs(p-XSYC);

%   ���. �����.-� ������ ������ P0(����.�����. � ���.����.)>=0
if (abs(CF) <= 1.E-12)
    if ((P0 - A2) > 1.E-12)
        %!C   ....��� �� ���������� ������
        SS=0.; return
    end %((P0 - A2) > 1.E-12)
    %    ............���������........
    SS = ywid*C0; return
end %(abs(CF) <= 1.E-12)


if (abs(SF) <= 1.E-12)
    if ((P0 - B2) > 1.E-12)
        %!C   ....��� �� ���������� ������
        SS=0.; return
    end %((P0 - B2) > 1.E-12)
    %      ..........�����������......
    SS = xwid*C0; return
end %(abs(SF) <= 1.E-12)

TF = SF/CF;	%! tangent
PC = P0/CF;
QP = B2+A2*TF;

if (PC >= QP)
    %!C   ....��� �� ���������� ������
    SS=0.; return
end %(PC >= QP)

QM = QP+PC;
if (QM > ywid)
    DEL = P0+B2*CF;
    if (DEL > A2*SF)
        %    ............���. �  �����. �������.....
        SS = (QP-PC)/SF*C0; return
    end % (DEL > A2*SF)
    
    %    ��� ���������� ����. � �����. �������
    SS = ywid/SF*C0; return
end %(QM > ywid)

%    ��� ���������� ���. � ����. �������
SS=xwid/CF*C0;

%Sin2DRect=SS

return