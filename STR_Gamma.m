R0=6.6*10^17;M=10^34;a=10000;dr=100;
B=eye(6);G=6.67*10^(-11);
for c=dr:dr:a
    k=c/dr
    h(k)=a-c;
    di1(k)=dinertia(R0,M,1,a,c);
    di2(k)=dinertia(R0,M,2,a,c);
    di3(k)=dinertia(R0,M,3,a,c);
end
figure,plot(h,di1,h,di2,h,di3);
xlabel("h(m)");
ylabel("Γ(s^2)");
legend("Γ_1","Γ_2","Γ_3",'Location','Best');
function [u,v,p,q,f,g]=solution(u1,u3,f2,u02,u04,f03,r,R0,M)
    G=6.67*10^(-11);
    u=u1*r+u3*r^3+u02*r^(-2)+u04*r^(-4);
    v=1/2*u1*r+5/6*u3*r^3-1/3*u04*r^(-4);
    q=M*(u1+8/3*u3*r^2+u02*r^(-3)+8/3*u04*r^(-5));
    f=f2*r^2+f03*r^(-3);
    g=2*f2*r-3*f03*r^(-4)+4*pi*G*R0*u;
    p=4/3*pi*G*R0^2*u*r-6*M*u/r+22*M*v/r-3*q+R0*f-M*(16/3*u3*r^2-3*u02*r^(-3)-40/3*u04*r^(-5));
end
function di=dinertia(R0,M,L,a,c)
    syms au1 au3 af2 bu1 bu3 bf2 bu02 bu04 bf03;
    [u,v,p,q,f,g]=solution(au1,au3,af2,0,0,0,c,R0,M);
    Xs1=[u;v;p;q;f;g];
    if L==1
            Xs1=Xs1+[0;0;0;-1/6/c^3;0;0];
        elseif L==2
            Xs1=Xs1+[0;0;3/2/c^3;-3/4/c^3;0;0];
        else
            Xs1=Xs1+[0;1/6/M/c^2;0;0;0;0];
    end
    [u,v,p,q,f,g]=solution(bu1,bu3,bf2,bu02,bu04,bf03,c,R0,M);
    Xs2=[u;v;p;q;f;g];
    [u,v,p,q,f,g]=solution(bu1,bu3,bf2,bu02,bu04,bf03,a,R0,M);
    Xs3=[u;v;p;q;f;g];
    eqn=[Xs1(1)-Xs2(1)==0,Xs1(2)-Xs2(2)==0,Xs1(3)-Xs2(3)==0,...
         Xs1(4)-Xs2(4)==0,Xs1(5)-Xs2(5)==0,Xs1(6)-Xs2(6)==0,...
         Xs3(3)==0,Xs3(4)==0,Xs3(6)*a+3*Xs3(5)==0];
    X=solve(eqn,[au1 au3 af2 bu1 bu3 bf2 bu02 bu04 bf03;]);
    AU1=double(X.au1);AU3=double(X.au3);AF2=double(X.af2);
    BU1=double(X.bu1);BU3=double(X.bu3);BF2=double(X.bf2);
    BU02=double(X.bu02);BU04=double(X.bu04);BF03=double(X.bf03);
    di=R0*(BU1*a^5+BU3*a^7+BU02*a^2+BU04);
end