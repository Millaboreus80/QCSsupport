R0=6.6*10^17;M=10^29;dr=10;d=0.1;a=10000;b=9500;c=b;
for n=c/dr:1:a/dr
    L=1;k=n-c/dr+1
    [UL,VS,GL]=boundary(R0,M,dr,d,n,L,a,b);
    di1(k)=dinertia(UL,VS,GL,R0,M,dr,d,n,L,a,b);
    h(k)=a-n*dr;
    L=2;
    [UL,VS,GL]=boundary(R0,M,dr,d,n,L,a,b);
    di2(k)=dinertia(UL,VS,GL,R0,M,dr,d,n,L,a,b);
    L=3;
    [UL,VS,GL]=boundary(R0,M,dr,d,n,L,a,b);
    di3(k)=dinertia(UL,VS,GL,R0,M,dr,d,n,L,a,b);
end
figure,plot(h,di1,h,di2,h,di3);
xlabel("h(m)");
ylabel("Γ(s^2)");
legend("Γ_1","Γ_2","Γ_3",'Location','Best');

function [UL,VS,GL]=boundary(R0,M,dr,d,n,L,a,b)
    B=eye(6);G=6.67*10^(-11);
    syms ul gl vs;
    F=(gl-4*pi*G*R0*ul)/2*b;
    p=R0*g(R0,b,d,a,b)*ul+R0*F;
    Xs=[ul;vs;p;0;F;gl;];
    for i=b/dr:1:n-1
        r=i*dr;
        A1=A(g(R0,r,d,a,b),R0,M,r,d)*B;
        A2=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(B+dr/2*A1);
        A3=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(B+dr/2*A2);
        A4=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d)*(B+dr*A3);
        B=B+(A1+2*A2+2*A3+A4)*dr/6;
    end
    i=n;
    if L==1
        Xs=B*Xs+[0;0;0;-1/6/(i*dr)^3;0;0];
    elseif L==2
        Xs=B*Xs+[0;0;3/2/(i*dr)^3;-3/4/(i*dr)^3;0;0];
    else
        Xs=B*Xs+[0;1/6/M/(i*dr)^2;0;0;0;0];
    end
    B=eye(6);
    for i=n:1:a/dr-1
        r=i*dr;
        A1=A(g(R0,r,d,a,b),R0,M,r,d)*B;
        A2=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(B+dr/2*A1);
        A3=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(B+dr/2*A2);
        A4=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d)*(B+dr*A3);
        B=B+(A1+2*A2+2*A3+A4)*dr/6;
    end
    Xs=B*Xs;
    eqn=[Xs(3)==0,Xs(4)==0,Xs(6)*a+3*Xs(5)==0];
    U=solve(eqn,[ul gl vs]);UL=double(U.ul);VS=double(U.vs);GL=double(U.gl);
end
function di=dinertia(UL,VS,GL,R0,M,dr,d,n,L,a,b)
    di=0;G=6.67*10^(-11);R1=R0*d;
    F=(GL-4*pi*G*R0*UL)/2*b;
    p=R0*g(R0,b,d,a,b)*UL+R0*F;
    XsN=[UL;VS;p;0;F;GL;];
    di=di+(b)^4*R0*XsN(1);
    di=di-(b)^4*R1*XsN(1);
    for i=b/dr:1:n-1
        r=i*dr;
        A1=A(g(R0,r,d,a,b),R0,M,r,d)*XsN;
        A2=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(XsN+dr/2*A1);
        A3=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(XsN+dr/2*A2);
        A4=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d)*(XsN+dr*A3);
        XsN=XsN+(A1+2*A2+2*A3+A4)*dr/6;
        %di=di+r^3*R1*(2*XsN(1)+6*XsN(2))*dr;
    end
    i=n;
    if L==1
        XsN=XsN+[0;0;0;-1/6/(i*dr)^3;0;0];
    elseif L==2
        XsN=XsN+[0;0;3/2/(i*dr)^3;-3/4/(i*dr)^3;0;0];
    else
        XsN=XsN+[0;1/6/M/(i*dr)^2;0;0;0;0];
    end
    for i=n:1:a/dr-1
        r=i*dr;
        A1=A(g(R0,r,d,a,b),R0,M,r,d)*XsN;
        A2=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(XsN+dr/2*A1);
        A3=A(g(R0,r+dr/2,d,a,b),R0,M,r+dr/2,d)*(XsN+dr/2*A2);
        A4=A(g(R0,r+dr,d,a,b),R0,M,r+dr,d)*(XsN+dr*A3);
        XsN=XsN+(A1+2*A2+2*A3+A4)*dr/6;
        %di=di+r^3*R1*(2*XsN(1)+6*XsN(2))*dr;
    end
    di=di+a^4*R1*XsN(1);
end%solve the Γ(h)
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
