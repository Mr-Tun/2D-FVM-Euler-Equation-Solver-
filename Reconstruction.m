function [UL,UR] = Reconstruction(U_1,U,U1,U2)
% Reconstruction function is used to reconstruct the solution vector
% UL =U;
% UR =U1;
%% define the U left
rl=zeros(1,4);
phil=zeros(1,4);
for i=1:4
    if U1(i)==U(i)  && U_1(i)==U(i)
        rl(i)=1;
    else
        rl(i)=(U1(i)-U(i))/(U(i)-U_1(i));
    end
    if rl(i)<0
        phil(i)=0;
    else
        phil(i)=min([1,rl(i)]);
    end 
end
UL = U + 0.5*phil.*(U-U_1);
%% define the U right
rr=zeros(1,4);
phir=zeros(1,4);
for i=1:4
    if U2(i)==U1(i)  && U1(i)==U(i)
        rr(i)=1;
    else
        rr(i)=(U1(i)-U(i))/(U2(i)-U1(i));
    end
    if rr(i)<0
        phir(i)=0;
    else
        phir(i)=min([1,rr(i)]);
    end 
end
UR = U1 - 0.5*phir.*(U2-U1);
end

