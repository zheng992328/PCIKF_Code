function x=Findwxz(Nroot,x0,dx,L,CorLength)
accuricywn=1e-5;
Nr=1;
x1=x0;
f0=fout(x0,CorLength,L);

while(Nr<=Nroot)
    x1=x1+dx;
    f1=fout(x1,CorLength,L);
    if (f1==0)
        x(Nr)=x1;
        y(Nr)=fout(x(Nr),CorLength,L);
        Nr=Nr+1;
        x0=x1+0.5*dx;
        f0=fout(x0,CorLength,L);
    elseif(f0*f1<0)
            x3=x1;
            while(abs(x0-x1)>accuricywn)
                x2=(x0+x1)*0.5;
                f2=fout(x2,CorLength,L);
                if(f2==0)
                    x(Nr)=x2;
                    y(Nr)=f2;
                    Nr=Nr+1;
                    x0=x1+0.5*dx;
                    f0=fout(x0,CorLength,L);
                    % go to 10?
                elseif(f1*f2<0)
                        x0=x2;
                        f0=f2;
                    else
                        x1=x2;
                        f1=f2;
                end
            end
            x(Nr)=x1-0.5*abs(x1-x0);
            y(Nr)=fout(x(Nr),CorLength,L);
            Nr=Nr+1;
            % 10 continue?
            x0=x3;
            f0=fout(x3,CorLength,L);
    end
end

end

