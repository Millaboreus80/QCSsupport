R0=6.6*10^17;M=10^34;a=10000;dw2=-1;
B=eye(6);G=6.67*10^(-11);
syms au1 au3 af2;
[u,v,p,q,f,g]=solution(au1,au3,af2,a,R0,M,dw2);
Xs=[u;v;p;q;f;g];
eqn=[Xs(3)==0,Xs(4)==0,Xs(6)*a+3*Xs(5)-5/3*dw2*a*a==0];
X=solve(eqn,[au1 au3 af2]);
AU1=double(X.au1);AU3=double(X.au3);AF2=double(X.af2);dr=100;
for c=dr:dr:10000
    rad=c/dr
    for THD=0:1:359
        if THD<=90
            TH=90-THD;
        elseif THD<=270
            TH=THD-90;
        else
            TH=450-THD;
        end
        TH=TH/180*pi;ang=THD+1;
        T1=stress(AU1,AU3,AF2,c,a,TH,R0,M,dw2,1)/M;
        T3=stress(AU1,AU3,AF2,c,a,TH,R0,M,dw2,3)/M;
        BTH=45;
        E(rad,ang)=(T3-T1)/sind(2*BTH)-(T3+T1)*cotd(2*BTH);
    end
end
pic(a,dr,E);
clear;
function [u,v,p,q,f,g]=solution(u1,u3,f2,r,R0,M,dw2)
    G=6.67*10^(-11);
    u=u1*r+u3*r^3;
    v=1/2*u1*r+5/6*u3*r^3;
    q=M*(u1+8/3*u3*r^2);
    f=f2*r^2;
    g=2*f2*r+4*pi*G*R0*u;
    p=4/3*pi*G*R0^2*u*r-6*M*u/r+22*M*v/r-3*q+R0*f-M*16/3*u3*r^2;
end%solve the dislocation at sphere centre
function E=strain(u1,u3,f2,r,a,TH,R0,M,dw2,Lam)
    [U,V,P,Q,~,~]=solution(u1,u3,f2,r,R0,M,dw2);
    EA=(-2*U+6*V)*1/2*(3*cos(TH)^2-1)/r;
    EB=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(2*TH));
    EC=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(TH)^2);
    ED=-Q*3*cos(TH)*sin(TH)/2/M;
    EE(1)=EA;EE(2)=EB;EE(3)=EC;EE(4)=ED;
    %EE(1)=1/2*(EA+EB+((EA-EB)^2+4*ED^2)^0.5);EE(2)=1/2*(EA+EB-((EA-EB)^2+4*ED^2)^0.5);EE(3)=EC; 
    %EE=sort(EE);
    E=EE(Lam);
end
function T=stress(u1,u3,f2,r,a,TH,R0,M,dw2,Lam)
    [U,V,P,Q,~,~]=solution(u1,u3,f2,r,R0,M,dw2);
    P0=-1/3*R0*dw2*(r*r-a*a);
    T1=P*1/2*(3*cos(TH)^2-1)-2*M/r*(-2*U+6*V)*1/2*(3*cos(TH)^2-1);
    TA=T1+2*M/r*(-2*U+6*V)*1/2*(3*cos(TH)^2-1)+P0;
    TB=T1+2*M/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(2*TH))+P0;
    TC=T1+2*M/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(TH)^2)+P0;
    TD=-Q*3*cos(TH)*sin(TH);
    TT(1)=1/2*(TA+TB+((TA-TB)^2+4*TD^2)^0.5);TT(2)=1/2*(TA+TB-((TA-TB)^2+4*TD^2)^0.5);TT(3)=TC; 
    %TT(1)=TA;TT(2)=TB;TT(3)=TC;TT(4)=TD;
    TT=sort(TT);
    T=TT(Lam);
end
function pic=pic(a,dr,T0)
    t = linspace(0, 2*pi, 360); 
    r=linspace(0,a,a/dr); 
    [tt, rr] = meshgrid(t, r);
    [x,y] = pol2cart(tt, rr); 
    figure
    hold on
    surf(x,y,T0);
    axis([-10000 10000 -10000 10000],"equal")
    shading interp
    colormap jet
    colorbar
end