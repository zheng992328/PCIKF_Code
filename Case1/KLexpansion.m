% KL expansion
function Fi = KLexpansion(Nrow,Ncol,Lx,Lz,CorLengthx,CorLengthz,VarF)
% Nrow=31;Ncol=31;
Npd=500;% the number in fun.dat
% Np=100;
Nroot=40;%超越方程的个数
% Lx=200;
% Lz=200;
% CorLengthx=100;
% CorLengthz=50;
% VarF=1.0;

dx=Lx/(Ncol-1)*ones(1,Ncol);
dz=Lz/(Nrow-1)*ones(1,Nrow);

% eigenpairAna_2D
X0RLn=0.1;
dXRLn=1.0e-5;
z0RLn=0.1;
dzRLn=1.0e-5;
Rx=CorLengthx/Lx;
Rz=CorLengthz/Lz;
wx=Findwxz(Nroot,X0RLn,dXRLn,1.0,Rx);
wz=Findwxz(Nroot,z0RLn,dzRLn,1.0,Rz);

for i=1:Nroot
    wx(i)=wx(i)/Lx;
    wz(i)=wz(i)/Lz;
end

RLn=zeros(3,Nroot^2);
for i=1:Nroot
    for j=1:Nroot
        k=(i-1)*Nroot+j;
        RLn(1,k)=wx(i);
        RLn(2,k)=wz(j);
        RLn(3,k)=4.0*CorLengthx*CorLengthz*VarF/((CorLengthx*wx(i))^2+1.)...
            /((CorLengthz*wz(j))^2+1.);
        if(k>1)
            for kk=k:-1:2
                if(RLn(3,kk)>RLn(3,kk-1))
                  r1=RLn(1,kk);
                  r2=RLn(2,kk);
                  r3=RLn(3,kk);
                  RLn(1,kk)=RLn(1,kk-1);
                  RLn(2,kk)=RLn(2,kk-1);
                  RLn(3,kk)=RLn(3,kk-1);
                  RLn(1,kk-1)=r1;
                  RLn(2,kk-1)=r2;
                  RLn(3,kk-1)=r3;
                end
            end
        end
    end
end

Fi=zeros(Npd,Ncol*Nrow);
for k=1:Npd
    z=0;
    for j=1:Nrow
        x=0;
        for i=1:Ncol
            N=(j-1)*Ncol+i;
            r1=CorLengthx*RLn(1,k);
            r2=CorLengthz*RLn(2,k);
            r4=RLn(1,k)*x;
            r5=RLn(2,k)*z;
            r3=(r1*cos(r4)+sin(r4))*(r2*cos(r5)+sin(r5))/...
                sqrt(((r1*r1+1.)*Lx/2.+CorLengthx)*...
                ((r2*r2+1.)*Lz/2.+CorLengthz));
            Fi(k,N)=sqrt(RLn(3,k))*r3;
            x=x+dx(i);
        end
        z=z+dz(j);
    end
end
    
Eigenv=zeros(1,Npd);
for i=1:Npd
     Eigenv(i)=RLn(3,i);
end

Fi=Fi';
            
    


