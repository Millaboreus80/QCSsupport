R0=6.6*10^17;%core density 
M=10^29;%mantle modulus
d=0.1;R1=R0*d;%mantle density
dw2=-1;dr=10;a=10000;b=9500;%core and mantle radius
G=6.67*10^(-11);
[UL,VS,GL]=boundary(R0,M,dr,d,a,b,dw2);x=0
E=NaN(a/dr,360);
for c=b:dr:a
    rad=c/dr;x=x+1
    for THD=0:1:359
        if THD<=90
            TH=90-THD;
        elseif THD<=270
            TH=THD-90;
        else
            TH=450-THD;
        end
        TH=TH/180*pi;ang=THD+1;
        E1=strain(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,3);
        %E3=strain(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,3);
        %T1=stress(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,1)/10^29;
        %T3=stress(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,3)/10^29;
        %BTH=45;
        %E(rad,ang)=(T3-T1)/sind(2*BTH)-(T3+T1)*cotd(2*BTH);
        E(rad,ang)=E1;
    end
end
pic(a,b,dr,E);
function [UL,VS,GL]=boundary(R0,M,dr,d,a,b,dw2)
    G=6.67*10^(-11);R1=R0*d;LP2=-2/3*dw2;LQ2=-1/3*dw2;
    syms ul gl vs;
    F=(gl-4*pi*G*R0*ul)/2*b;
    p=R0*g(R0,b,d,a,b)*ul+R0*F-R0*LQ2*b*b;
    Xs=[ul;vs;p;0;F;gl;];
    for i=b/dr:1:a/dr-1
        r=i*dr;k=i-b/dr+1
        a0=A(g(R0,r,d,a,b),R0,M,r,d);a1=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d);a2=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d);
        A1=a0;A2=a1+a1*A1*dr/2;A3=a1+a1*A2*dr/2;A4=a2+a2*A3*dr;
        B1=[0;0;-R1*LP2*r;-R1*LQ2*r;0;0];B2=a1*B1*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B3=a1*B2*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B4=a2*B3*dr+[0;0;-R1*LP2*(r+dr);-R1*LQ2*(r+dr);0;0];
        Xs=(eye(6)+(A1+2*A2+2*A3+A4)*dr/6)*Xs+(B1+2*B2+2*B3+B4)*dr/6;
    end
    eqn=[Xs(3)==0,Xs(4)==0,Xs(6)*a+3*Xs(5)==0];
    U=solve(eqn,[ul gl vs]);UL=double(U.ul);VS=double(U.vs);GL=double(U.gl);
end%solve the dislocation at core-mantle surface
function E=strain(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,Lam)
    G=6.67*10^(-11);R1=R0*d;LP2=-2/3*dw2;LQ2=-1/3*dw2;
    F=(GL-4*pi*G*R0*UL)/2*b;
    p=R0*g(R0,b,d,a,b)*UL+R0*F-R0*LQ2*b*b;
    XsN=[UL;VS;p;0;F;GL;];
    for i=b/dr:1:c/dr-1
        r=i*dr;
        a0=A(g(R0,r,d,a,b),R0,M,r,d);a1=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d);a2=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d);
        A1=a0;A2=a1+a1*A1*dr/2;A3=a1+a1*A2*dr/2;A4=a2+a2*A3*dr;
        B1=[0;0;-R1*LP2*r;-R1*LQ2*r;0;0];B2=a1*B1*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B3=a1*B2*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B4=a2*B3*dr+[0;0;-R1*LP2*(r+dr);-R1*LQ2*(r+dr);0;0];
        XsN=(eye(6)+(A1+2*A2+2*A3+A4)*dr/6)*XsN+(B1+2*B2+2*B3+B4)*dr/6;
    end
    r=c;U2=XsN(1);V2=XsN(2);Q2=XsN(4);
    EA=(-2*U2+6*V2)*1/2*(3*cos(TH)^2-1)/r;
    EB=(U2*(3*cos(TH)^2-1)/2-V2*3*cos(2*TH))/r;
    EC=(U2*(3*cos(TH)^2-1)/2-V2*3*cos(TH)^2)/r;
    ED=-Q2*3*cos(TH)*sin(TH)/2/M;
    EE(1)=EA;EE(2)=EB;EE(3)=EC;
    %EE=sort(EE);
    E=EE(Lam);
end%solve the strain, Lam means the subscript of strain component.
function T=stress(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH,Lam)
    G=6.67*10^(-11);R1=R0*d;LP2=-2/3*dw2;LQ2=-1/3*dw2;LP0=2/3*dw2;
    F=(GL-4*pi*G*R0*UL)/2*b;
    p=R0*g(R0,b,d,a,b)*UL+R0*F-R0*LQ2*b*b;
    XsN=[UL;VS;p;0;F;GL;];
    for i=b/dr:1:c/dr-1
        r=i*dr;
        a0=A(g(R0,r,d,a,b),R0,M,r,d);a1=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d);a2=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d);
        A1=a0;A2=a1+a1*A1*dr/2;A3=a1+a1*A2*dr/2;A4=a2+a2*A3*dr;
        B1=[0;0;-R1*LP2*r;-R1*LQ2*r;0;0];B2=a1*B1*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B3=a1*B2*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B4=a2*B3*dr+[0;0;-R1*LP2*(r+dr);-R1*LQ2*(r+dr);0;0];
        XsN=(eye(6)+(A1+2*A2+2*A3+A4)*dr/6)*XsN+(B1+2*B2+2*B3+B4)*dr/6;
    end
    r=c;U2=XsN(1);V2=XsN(2);P2=XsN(3);Q2=XsN(4);P0=-1/2*R1*LP0*(r*r-a*a);
    T1=P2*1/2*(3*cos(TH)^2-1)-2*M/r*(-2*U2+6*V2)*1/2*(3*cos(TH)^2-1);
    TA=T1+2*M/r*(-2*U2+6*V2)*1/2*(3*cos(TH)^2-1)+P0;
    TB=T1+2*M/r*(U2*(3*cos(TH)^2-1)/2-V2*3*cos(2*TH))+P0;
    TC=T1+2*M/r*(U2*(3*cos(TH)^2-1)/2-V2*3*cos(TH)^2)+P0;
    TD=-Q2*3*cos(TH)*sin(TH);
    TT(1)=TA;TT(2)=TB;TT(3)=TC;
    %TT=sort(TT);
    T=TT(Lam);
end%solve the stress, Lam means the subscript of strain component.
function A=A(g,R0,M,r,d)
    G=6.67*10^(-11);K=3*M;R1=R0*d;
    A=[ -2/r      6/r    0   0   0   0;...
        -1/r  1/r   0   1/M 0   0;...
        4*K/r/r-4*R1*g/r      -6*(2*K/r/r-R1*g/r)   0    6/r   0   R1;...
        R1*g/r-2*K/r/r        (6*K+4*M)/r/r    -1/r         -3/r  R1/r     0;...
        -4*pi*G*R1     0   0   0   0   1;...
        0   24*pi*G*R1/r     0   0   6/r/r   -2/r;];
end
function g=g(R0,r,d,a,b)
    G=6.67*10^(-11);
    R1=R0*d;
    if r<b
        g=4/3*pi*G*R0*r;
    elseif r>a
        g=(4/3*pi*(R0*b^3+R1*(a^3-b^3)))*G/r/r;
    else
        g=(4/3*pi*(R0*b^3+R1*(r^3-b^3)))*G/r/r;
    end
end
function pic(a,b,dr,T0)
    t = linspace(0, 2*pi, 360);
    r=linspace(0,a,a/dr); 
    for i=1:1:a/dr
        if r(i)<b
            r(i)=r(i)/3;
        else
            r(i)=(r(i)-b)/(a-b)*(a-b/3)+b/3;
        end
    end
    [tt, rr] = meshgrid(t, r);
    [x,y] = pol2cart(tt, rr);
    figure
    hold on
    surf(x,y,T0);
    axis([-a a -a a],"equal")
    axis off
    shading interp
    colormap jet
    colorbar
end