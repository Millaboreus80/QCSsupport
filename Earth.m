M=zeros(1,199);L=zeros(1,199);g=zeros(1,199);m=0;G=6.67*10^(-11);PG=pi*G;
for i=1:1:198
    M(i)=EM(i,4)*EM(i,4)*EM(i,2);
    L(i)=EM(i,3)*EM(i,3)*EM(i,2)-2*M(i);
    if i==1
        m=m+EM(i,2)*EM(i,1)^3*4*pi/3;
    else
        m=m+EM(i,2)*(EM(i,1)^3-EM(i-1,1)^3)*4*pi/3;
    end
    g(i)=G*m/EM(i,1)/EM(i,1);
end%
syms ul gl vs;
Xl(:,1)=[ul;0;0;gl];B=eye(4);
for i=1:1:13
    dr=EM(i+1,1)-EM(i,1);r=EM(i,1);l=L(i);ro=EM(i,2);RO=EM(i+1,2);G0=g(i);G1=g(i+1);
    A=[     (2*G0*ro*log(1+dr/r)+G0*ro-4*PG*RO*dr*ro)       (-1+G0*dr*ro/l)      (RO)        (RO*dr-dr*ro);... 
                (G0*G1*RO*ro+2*G0*G1*RO*log(1+dr/r)*ro-4*G0*PG*RO*dr*ro^2)      (-G1*RO+G0*G1*RO*dr/l*ro)      (G0*RO*ro)         (G0*RO*dr*ro-G1*RO*dr*ro);...
                (-4*PG*dr*ro*(G0*ro-G1*RO))       0       (G0*ro-G1*RO)       (dr*(G0*ro-G1*RO));...
                (-16*PG^2*RO*dr*ro^2+4*G1*PG*RO*ro+16*G0*PG*log(1+dr/r)*ro^2-8*G1*PG*RO*log(1+dr/r)*ro)     (-4*PG*ro+4*G1*PG*RO*dr/l*ro)       (6*G0*dr/r/(r+dr)*ro-6*G1*RO*dr/r/(r+dr)+4*PG*RO*ro)      (-G1*RO+G0*ro-4*PG*dr*ro^2+2*G1*RO*log(1+dr/r)-2*G0*log(1+dr/r)*ro+4*PG*RO*dr*ro);];
    A=A/(G0*ro-G1*RO);
    B=A*B;
end
for i=15:1:38
    dr=EM(i+1,1)-EM(i,1);r=EM(i,1);l=L(i);ro=EM(i,2);RO=EM(i+1,2);G0=g(i);G1=g(i+1);
    A=[     (2*G0*ro*log(1+dr/r)+G0*ro-4*PG*RO*dr*ro)       (-1+G0*dr*ro/l)      (RO)        (RO*dr-dr*ro);... 
                (G0*G1*RO*ro+2*G0*G1*RO*log(1+dr/r)*ro-4*G0*PG*RO*dr*ro^2)      (-G1*RO+G0*G1*RO*dr/l*ro)      (G0*RO*ro)         (G0*RO*dr*ro-G1*RO*dr*ro);...
                (-4*PG*dr*ro*(G0*ro-G1*RO))       0       (G0*ro-G1*RO)       (dr*(G0*ro-G1*RO));...
                (-16*PG^2*RO*dr*ro^2+4*G1*PG*RO*ro+16*G0*PG*log(1+dr/r)*ro^2-8*G1*PG*RO*log(1+dr/r)*ro)     (-4*PG*ro+4*G1*PG*RO*dr/l*ro)       (6*G0*dr/r/(r+dr)*ro-6*G1*RO*dr/r/(r+dr)+4*PG*RO*ro)      (-G1*RO+G0*ro-4*PG*dr*ro^2+2*G1*RO*log(1+dr/r)-2*G0*log(1+dr/r)*ro+4*PG*RO*dr*ro);];
    A=A/(G0*ro-G1*RO);
    B=A*B;
end
Xl(:,2)=B*Xl(:,1);
Xs=[Xl(1,2);vs;Xl(2,2);0;Xl(3,2);Xl(4,2)];
for n=87:1:197 
    C=eye(6);n
    for i=40:1:n
        B=L(i)+2*M(i);K=L(i)+M(i)-L(i)*L(i)/B;
        A=[ -2*L(i)/EM(i,1)/B       6*L(i)/EM(i,1)/B    1/B     0   0   0;...
            -1/EM(i,1)  1/EM(i,1)   0   1/M(i)  0   0;...
            4*K/EM(i,1)/EM(i,1)-4*EM(i,2)*g(i)/EM(i,1)      -6*(2*K/EM(i,1)/EM(i,1)-EM(i,2)*g(i)/EM(i,1))   2*(L(i)/B-1)/EM(i,1)    6/EM(i,1)   0   EM(i,2);...
            EM(i,2)*g(i)/EM(i,1)-2*K/EM(i,1)/EM(i,1)        (6*K+4*M(i))/EM(i,1)/EM(i,1)    -L(i)/B/EM(i,1)         -3/EM(i,1)      EM(i,2)/EM(i,1)     0;...
            -4*pi*G*EM(i,2)     0   0   0   0   1;...
            0   24*pi*G*EM(i,2)/EM(i,1)     0   0   6/EM(i,1)/EM(i,1)   -2/EM(i,1)];
        dr=EM(i+1,1)-EM(i,1);
        C=(eye(6)+A*dr)*C;
    end
    i=n+1;B=L(i)+2*M(i);
    Xs1=C*Xs+[0;0;0;-1/6/EM(i,1)^3;0;0];
    Xs2=C*Xs+[0;0;3/2/EM(i,1)^3;-1/4/EM(i,1)^3;0;0]+[1/B/(EM(i,1))^2;0;2*(L(i)/B-1)/(EM(i,1))^3;-L(i)/B/(EM(i,1))^3;0;0]/2;
    Xs3=C*Xs+[0;0;-1/(EM(i,1))^3;1/2/EM(i,1)^3;0;0]+[0;1/M(i)/(EM(i,1))^2;6/(EM(i,1))^3;-3/(EM(i,1))^3;0;0]/6;
    C=eye(6);
    for i=n+1:1:197
        B=L(i)+2*M(i);K=L(i)+M(i)-L(i)*L(i)/B;
        A=[ -2*L(i)/EM(i,1)/B       6*L(i)/EM(i,1)/B    1/B     0   0   0;...
            -1/EM(i,1)  1/EM(i,1)   0   1/M(i)  0   0;...
            4*K/EM(i,1)/EM(i,1)-4*EM(i,2)*g(i)/EM(i,1)      -6*(2*K/EM(i,1)/EM(i,1)-EM(i,2)*g(i)/EM(i,1))   2*(L(i)/B-1)/EM(i,1)    6/EM(i,1)   0   EM(i,2);...
            EM(i,2)*g(i)/EM(i,1)-2*K/EM(i,1)/EM(i,1)        (6*K+4*M(i))/EM(i,1)/EM(i,1)    -L(i)/B/EM(i,1)         -3/EM(i,1)      EM(i,2)/EM(i,1)     0;...
            -4*pi*G*EM(i,2)     0   0   0   0   1;...
            0   24*pi*G*EM(i,2)/EM(i,1)     0   0   6/EM(i,1)/EM(i,1)   -2/EM(i,1)];
        dr=EM(i+1,1)-EM(i,1);
        C=(eye(6)+A*dr)*C;
    end
    i=198;k=n-86;
    Xs1=C*Xs1;
    eqn=[Xs1(3)==0,Xs1(4)==0,Xs1(6)*EM(i,1)+3*Xs1(5)==0];
    U=solve(eqn,[ul gl vs]);
    UL=double(U.ul);
    VS=double(U.vs);
    GL=double(U.gl);
    dc1(k)=DC(UL,VS,GL,g,n,EM,L,M,1);
    eqn=[Xs2(3)==0,Xs2(4)==0,Xs2(6)*EM(i,1)+3*Xs2(5)==0];
    U=solve(eqn,[ul gl vs]);
    UL=double(U.ul);
    VS=double(U.vs);
    GL=double(U.gl);
    dc2(k)=DC(UL,VS,GL,g,n,EM,L,M,2);
    eqn=[Xs3(3)==0,Xs3(4)==0,Xs3(6)*EM(i,1)+3*Xs3(5)==0];
    U=solve(eqn,[ul gl vs]);
    UL=double(U.ul);
    VS=double(U.vs);
    GL=double(U.gl);
    dc3(k)=DC(UL,VS,GL,g,n,EM,L,M,3);
    h(k)=6371-EM(n,1)/1000;
end
figure,plot(h,dc1,h,-dc2,h,dc3);
xlabel("h(km)");
ylabel("Γ(s^2)");
legend("Γ_1","-Γ_2","Γ_3");
function DC=DC(UL,VS,GL,g,n,EM,L,M,Lam)
    XN=[UL;0;0;GL;];G=6.67*10^(-11);DC=0;PG=pi*G;
    for i=1:1:13
        dr=EM(i+1,1)-EM(i,1);r=EM(i,1);l=L(i);ro=EM(i,2);RO=EM(i+1,2);G0=g(i);G1=g(i+1);
        A=[     (2*G0*ro*log(1+dr/r)+G0*ro-4*PG*RO*dr*ro)       (-1+G0*dr*ro/l)      (RO)        (RO*dr-dr*ro);... 
                    (G0*G1*RO*ro+2*G0*G1*RO*log(1+dr/r)*ro-4*G0*PG*RO*dr*ro^2)      (-G1*RO+G0*G1*RO*dr/l*ro)      (G0*RO*ro)         (G0*RO*dr*ro-G1*RO*dr*ro);...
                    (-4*PG*dr*ro*(G0*ro-G1*RO))       0       (G0*ro-G1*RO)       (dr*(G0*ro-G1*RO));...
                    (-16*PG^2*RO*dr*ro^2+4*G1*PG*RO*ro+16*G0*PG*log(1+dr/r)*ro^2-8*G1*PG*RO*log(1+dr/r)*ro)     (-4*PG*ro+4*G1*PG*RO*dr/l*ro)       (6*G0*dr/r/(r+dr)*ro-6*G1*RO*dr/r/(r+dr)+4*PG*RO*ro)      (-G1*RO+G0*ro-4*PG*dr*ro^2+2*G1*RO*log(1+dr/r)-2*G0*log(1+dr/r)*ro+4*PG*RO*dr*ro);];
        A=A/(G0*ro-G1*RO);
        XN=A*XN;
        DC=DC-dr*EM(i+1,1)^4*EM(i+1,2)^2*XN(3)/L(i+1);
    end
    for i=15:1:38
        dr=EM(i+1,1)-EM(i,1);r=EM(i,1);l=L(i);ro=EM(i,2);RO=EM(i+1,2);G0=g(i);G1=g(i+1);
        A=[     (2*G0*ro*log(1+dr/r)+G0*ro-4*PG*RO*dr*ro)       (-1+G0*dr*ro/l)      (RO)        (RO*dr-dr*ro);... 
                    (G0*G1*RO*ro+2*G0*G1*RO*log(1+dr/r)*ro-4*G0*PG*RO*dr*ro^2)      (-G1*RO+G0*G1*RO*dr/l*ro)      (G0*RO*ro)         (G0*RO*dr*ro-G1*RO*dr*ro);...
                    (-4*PG*dr*ro*(G0*ro-G1*RO))       0       (G0*ro-G1*RO)       (dr*(G0*ro-G1*RO));...
                    (-16*PG^2*RO*dr*ro^2+4*G1*PG*RO*ro+16*G0*PG*log(1+dr/r)*ro^2-8*G1*PG*RO*log(1+dr/r)*ro)     (-4*PG*ro+4*G1*PG*RO*dr/l*ro)       (6*G0*dr/r/(r+dr)*ro-6*G1*RO*dr/r/(r+dr)+4*PG*RO*ro)      (-G1*RO+G0*ro-4*PG*dr*ro^2+2*G1*RO*log(1+dr/r)-2*G0*log(1+dr/r)*ro+4*PG*RO*dr*ro);];
        A=A/(G0*ro-G1*RO);
        XN=A*XN;
        DC=DC-dr*EM(i+1,1)^4*EM(i+1,2)^2*XN(3)/L(i+1);
    end
    i=39;
    DC=DC+EM(i,1)^4*EM(i,2)*XN(1);
    XsN=[XN(1);VS;XN(2);0;XN(3);XN(4)];
    for i=40:1:n
        B=L(i)+2*M(i);K=L(i)+M(i)-L(i)*L(i)/B;
        A=[ -2*L(i)/EM(i,1)/B       6*L(i)/EM(i,1)/B    1/B     0   0   0;...
            -1/EM(i,1)  1/EM(i,1)   0   1/M(i)  0   0;...
            4*K/EM(i,1)/EM(i,1)-4*EM(i,2)*g(i)/EM(i,1)      -6*(2*K/EM(i,1)/EM(i,1)-EM(i,2)*g(i)/EM(i,1))   2*(L(i)/B-1)/EM(i,1)    6/EM(i,1)   0   EM(i,2);...
            EM(i,2)*g(i)/EM(i,1)-2*K/EM(i,1)/EM(i,1)        (6*K+4*M(i))/EM(i,1)/EM(i,1)    -L(i)/B/EM(i,1)         -3/EM(i,1)      EM(i,2)/EM(i,1)     0;...
            -4*pi*G*EM(i,2)     0   0   0   0   1;...
            0   24*pi*G*EM(i,2)/EM(i,1)     0   0   6/EM(i,1)/EM(i,1)   -2/EM(i,1)];
        dr=EM(i+1,1)-EM(i,1);
        U1=XsN(1);
        XsN=(eye(6)+A*dr)*XsN;
        U2=XsN(1);
        DC=DC+dr*EM(i+1,1)^3*EM(i+1,2)*(2*XsN(1)+6*XsN(2));
    end
    i=n+1;B=L(i)+2*M(i);
    if Lam==1
        XsN=XsN+[0;0;0;-1/6/EM(i,1)^3;0;0];
    elseif Lam==2
        XsN=XsN+[0;0;3/2/EM(i,1)^3;-1/4/EM(i,1)^3;0;0]+[1/B/(EM(i,1))^2;0;2*(L(i)/B-1)/(EM(i,1))^3;-L(i)/B/(EM(i,1))^3;0;0]/2;
    else
        XsN=[0;0;-1/(EM(i,1))^3;1/2/EM(i,1)^3;0;0]+[0;1/M(i)/(EM(i,1))^2;6/(EM(i,1))^3;-3/(EM(i,1))^3;0;0]/6;
    end
    for i=n+1:1:197
        B=L(i)+2*M(i);K=L(i)+M(i)-L(i)*L(i)/B;
        A=[ -2*L(i)/EM(i,1)/B       6*L(i)/EM(i,1)/B    1/B     0   0   0;...
            -1/EM(i,1)  1/EM(i,1)   0   1/M(i)  0   0;...
            4*K/EM(i,1)/EM(i,1)-4*EM(i,2)*g(i)/EM(i,1)      -6*(2*K/EM(i,1)/EM(i,1)-EM(i,2)*g(i)/EM(i,1))   2*(L(i)/B-1)/EM(i,1)    6/EM(i,1)   0   EM(i,2);...
            EM(i,2)*g(i)/EM(i,1)-2*K/EM(i,1)/EM(i,1)        (6*K+4*M(i))/EM(i,1)/EM(i,1)    -L(i)/B/EM(i,1)         -3/EM(i,1)      EM(i,2)/EM(i,1)     0;...
            -4*pi*G*EM(i,2)     0   0   0   0   1;...
            0   24*pi*G*EM(i,2)/EM(i,1)     0   0   6/EM(i,1)/EM(i,1)   -2/EM(i,1)];
        dr=EM(i+1,1)-EM(i,1);
        U1=XsN(1);
        XsN=(eye(6)+A*dr)*XsN;
        U2=XsN(1);
        DC=DC+dr*EM(i+1,1)^3*EM(i+1,2)*(2*XsN(1)+6*XsN(2));
    end
    i=198;
end