clear
clc
nu=.48;  %Poisson's ratio
E=2e3;% youngs modulus

%%%%%%%% thickness in microns
t=6e-6;
%%%%%%%%%
Ai=2401.91e-12;%initial area of pattern in m^2
Aff=[2252.47];%final area of pattern in Um^2

for i=1:length(Aff)
    
Af=Aff(i)*1e-12;
l=sqrt(Ai);


syms p;
%eq=-nu*(p/E)^2+(1-nu)*p/E+(1-(Af/Ai));
eq=-nu*(p/E)^2+(1-nu)*p/E+(1-(Af/Ai));

P=solve(eq,p);
P=eval(P);
f1=P(1)*l*t;
%f2=P(2)*l*t;

Dl1_1=l*P(1)/E;
%Dl1_2=l*P(2)/E;

Dl2_1=-nu*l*P(1)/E;
%Dl2_2=-nu*l*P(2)/E;


U1=f1*Dl1_1*1e12;
%U2=f2*Dl1_2*1e12;

U11(i)=U1;
%U22(i)=U2;
end

%%%TFM values
% U_tfm=[0.319501817];

Ar=Aff/Ai;
Strainenergy=U11;
% scatter(Aff,U_tfm)
% hold on
% scatter(Aff,U11)
%scatter(Aff,U22)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Af=2232.059*1e-12;
% l=sqrt(Ai);
% t=300e-6;
% 
% syms p;
% eq=-nu/E^2*p^2+(1-nu)*p/E+(1-Af/Ai);
% P=solve(eq,p);
% P=eval(P);
% f1=P(1)*l*t;
% f2=P(2)*l*t;
% 
% Dl1_1=l*P(1)/E;
% Dl1_2=l*P(2)/E;
% 
% Dl2_1=-nu*l*P(1)/E;
% Dl2_2=-nu*l*P(2)/E;
% 
% bound=sqrt(Ai)/2;
% int=1e-6;
% 
% X=-bound:int:bound;
% X(length(X))=bound;
% Y=-bound:int:bound;
% Y(length(X))=bound;
% [x,y]=meshgrid(X,Y);
% 
% u=Dl1_2*2.*x/l;
% v=Dl2_2*2.*y/l;
% 
% sx=P(2);
% sy=-nu*P(2);
% 
% ubar1=u*sx;
% ubar2=v*sy;
% U1=trapz(Y,trapz(X,ubar1));
% U2=trapz(Y,trapz(X,ubar2));
% U=0.5*(U1+U2)*1e12
% 
% 
% 
% 
% 
% 
