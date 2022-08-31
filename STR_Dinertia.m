R0=6.6*10^17;M=10^34;a=10000;dw2=-1;
B=eye(6);G=6.67*10^(-11);
syms au1 au3 af2;
[u,v,p,q,f,g]=solution(au1,au3,af2,a,R0,M);
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
        [E1,E2,E3,DTH]=strain(AU1,AU3,AF2,c,TH,R0,M);
        [Al,D,L,so]=sort(E1,E2,E3,DTH,60);n=c/dr;k=c/dr;
        di(rad,ang)=-(2*sind(2*Al)*sind(D)*cosd(L)+cosd(2*Al)*sind(2*D)*sind(L))*sin(TH)^2*di1(k)...
            -2*sind(L)*sind(2*D)*(cos(TH)^2-1/3)*di2(k)+...
           (cosd(Al)*cosd(D)*cosd(L)-sind(Al)*cosd(2*D)*sind(L))*sin(2*TH)*di3(k);
    end
end
pic(a,dr,-di);
function [u,v,p,q,f,g]=solution(u1,u3,f2,r,R0,M)
    G=6.67*10^(-11);
    u=u1*r+u3*r^3;
    v=1/2*u1*r+5/6*u3*r^3;
    q=M*(u1+8/3*u3*r^2);
    f=f2*r^2;
    g=2*f2*r+4*pi*G*R0*u;
    p=4/3*pi*G*R0^2*u*r-6*M*u/r+22*M*v/r-3*q+R0*f-M*16/3*u3*r^2;
end
function [E1,E2,E3,DTH]=strain(u1,u3,f2,r,TH,R0,M)
    [U,V,~,Q,~,~]=solution(u1,u3,f2,r,R0,M);
    EA=(-2*U+6*V)*1/2*(3*cos(TH)^2-1)/r;
    EB=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(2*TH));
    EC=1/r*(U*(3*cos(TH)^2-1)/2-V*3*cos(TH)^2);
    ED=-Q*3*cos(TH)*sin(TH)/2/M;
    %E1=1/2*(EA+EB+((EA-EB)^2+4*ED^2)^0.5);E2=1/2*(EA+EB-((EA-EB)^2+4*ED^2)^0.5);E3=EC;
    if ED==0
        E1=EA;E2=EB;E3=EC;DTH=0;
    else
        A=(EA-EB)/2/ED;
        E1=1/2*(EA+EB)+ED*(1+A*A)^0.5;E2=1/2*(EA+EB)-ED*(1+A*A)^0.5;E3=EC;
        si=1/(1+A*A)^0.5;co=A/(1+A*A)^0.5;
        DTH=angle(si,co)/2;
    end
end
function pic(a,dr,T0)
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
end
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
    