R0=6.6*10^17;M=10^29;dr=10;d=0.1;a=10000;b=9500;dw2=-1;G=6.67*10^(-11);R1=R0*d;
[UL,VS,GL]=boundary(R0,M,dr,d,a,b,dw2);x=0
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
        [E1,E2,E3,DTH]=strain(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH);
        [Al,D,L,so]=sort(E1,E2,E3,DTH,45);
        n=c/dr;k=(a-c)/dr+1;
        di(rad,ang)=-(2*sind(2*Al)*sind(D)*cosd(L)+cosd(2*Al)*sind(2*D)*sind(L))*sin(TH)^2*di1(k)...
            -2*sind(L)*sind(2*D)*(cos(TH)^2-1/3)*di2(k)+...
            (cosd(Al)*cosd(D)*cosd(L)-sind(Al)*cosd(2*D)*sind(L))*sin(2*TH)*di3(k);
        dii(rad,ang)=so;
    end
end
pic(a,b,dr,di);
function [UL,VS,GL]=boundary(R0,M,dr,d,a,b,dw2)
    G=6.67*10^(-11);R1=R0*d;LP2=-2/3*dw2;LQ2=-1/3*dw2;
    syms ul gl vs;
    F=(gl-4*pi*G*R0*ul)/2*b;
    p=R0*g(R0,b,d,a,b)*ul+R0*F-R0*LQ2*b*b;
    Xs=[ul;vs;p;0;F;gl;];
    for i=b/dr:1:a/dr-1
        r=i*dr;
        a0=A(g(R0,r,d,a,b),R0,M,r,d);a1=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d);a2=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d);
        A1=a0;A2=a1+a1*A1*dr/2;A3=a1+a1*A2*dr/2;A4=a2+a2*A3*dr;
        B1=[0;0;-R1*LP2*r;-R1*LQ2*r;0;0];B2=a1*B1*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B3=a1*B2*dr/2+[0;0;-R1*LP2*(r+dr/2);-R1*LQ2*(r+dr/2);0;0];B4=a2*B3*dr+[0;0;-R1*LP2*(r+dr);-R1*LQ2*(r+dr);0;0];
        Xs=(eye(6)+(A1+2*A2+2*A3+A4)*dr/6)*Xs+(B1+2*B2+2*B3+B4)*dr/6;
    end
    eqn=[Xs(3)==0,Xs(4)==0,Xs(6)*a+3*Xs(5)==0];
    U=solve(eqn,[ul gl vs]);UL=double(U.ul);VS=double(U.vs);GL=double(U.gl);
end
function [E1,E2,E3,DTH]=strain(UL,VS,GL,R0,M,dr,d,a,b,c,dw2,TH)
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
    r=c;U=XsN(1);V=XsN(2);Q=XsN(4);
    EA=(-2*U+6*V)*1/2*(3*cos(TH)^2-1)/r;
    EB=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(2*TH));
    EC=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(TH)^2);
    ED=-Q*3*cos(TH)*sin(TH)/2/M;
    %E1=1/2*(EA+EB+((EA-EB)^2+4*ED^2)^0.5);E2=1/2*(EA+EB-((EA-EB)^2+4*ED^2)^0.5);E3=EC;
    if ED==0
        E1=EA;E2=EB;E3=EC;DTH=0;
    else
        Aa=(EA-EB)/2/ED;
        E1=1/2*(EA+EB)+ED*(1+Aa*Aa)^0.5;E2=1/2*(EA+EB)-ED*(1+Aa*Aa)^0.5;E3=EC;
        si=1/(1+Aa*Aa)^0.5;co=Aa/(1+Aa*Aa)^0.5;
        DTH=angle(si,co)/2;
    end
end%solve the strain and seismic parameters
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
function [A,D,L,so]=sort(E1,E2,E3,DTH,BTH)
        if E1>=E2
            if E2>=E3
                si=(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=-sind(BTH)*cosd(DTH);
                D=angle(si,co);
                si=-sind(BTH)*sind(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=-cosd(BTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                A=angle(si,co);
                si=-cosd(BTH)*cosd(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=sind(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                L=angle(si,co);%123
                so=123;
            elseif E1>=E3
                si=(cosd(BTH+DTH)^2)^0.5;
                co=-sind(BTH+DTH);
                D=angle(si,co);
                if cosd(DTH+BTH)>=0
                    A=90;
                    L=270;
                else
                    A=270;
                    L=90;
                end%132
                so=132;
            else
                si=(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=-cosd(BTH)*sind(DTH);
                D=angle(si,co);
                si=cosd(BTH)*cosd(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=sind(BTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                A=angle(si,co);
                si=sind(BTH)*sind(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=-cosd(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                L=angle(si,co);%312
                so=312;
            end
        elseif E1>=E3
                si=(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=-cosd(BTH)*sind(DTH);
                D=angle(si,co);
                si=cosd(BTH)*cosd(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=sind(BTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                A=angle(si,co);
                si=sind(BTH)*sind(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                co=-cosd(DTH)/(sind(BTH)^2+cosd(BTH)^2*cosd(DTH)^2)^0.5;
                L=180+angle(si,co);%213
                so=213;
        else
            if E2>=E3
                si=(cosd(BTH+DTH)^2)^0.5;
                co=-sind(BTH+DTH);
                D=angle(si,co);
                if cosd(DTH+BTH)>=0
                    A=90;
                    L=90;
                else
                    A=270;
                    L=270;
                end%231
                so=231;
            else
                si=(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=-sind(BTH)*cosd(DTH);
                D=angle(si,co);
                si=-sind(BTH)*sind(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=-cosd(BTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                A=angle(si,co);
                si=-cosd(BTH)*cosd(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                co=sind(DTH)/(cosd(BTH)^2+sind(BTH)^2*sind(DTH)^2)^0.5;
                L=180+angle(si,co);%321
                so=321;
            end
        end
end%solve the strain and seismic parameters
function theta=angle(si,co)
    if si>=0
        if co>=0
            theta=asind(si);
        else
            theta=180-asind(si);
        end
    elseif co>=0
        theta=asind(si);
    else
        theta=180-asind(si);
    end
end