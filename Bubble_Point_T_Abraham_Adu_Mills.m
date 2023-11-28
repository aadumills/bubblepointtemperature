% Equation of state for vapor phase, Gamma model for liquid phase
% Bubble point temperature calculation
clear
clc
P=101325; %Pa
% component 1 cyclohexane(C6H12), component 2 tert-butanol(C4H10O)
TC= [554 506.2]; %critical temperature in Kelvin
PC= [40.7e5 39.7e5]; %critical pressure in Pascal
omg= [0.20926 0.6252];   %accentric factors
x1=linspace(0.001,0.999,100);
x2=1-x1;
x=zeros(2,100);
x(1,:)=x1;
x(2,:)=x2;
T_bp=linspace(0,0,100);     %bubble point temperature
y=zeros(2,100);


for i=1:100
    T=323.15;
    y(:,i)=x(:,i);
    it=0;
    for iter=1:50
        it=it+1;
        kk= valik(x(:,i), y(:,i),T,P,TC,PC,omg); 
        w=kk'.*x(:,i);
        S=sum(w);
        FO=log(S);
        if abs(FO)<0.001
            y(:,i)=w./S;
            break
        else
            TT=T+1;
        kk= valik(x(:,i), y(:,i),TT,P,TC,PC,omg); 
        w=kk'.*x(:,i);
        SS=sum(w);
            F1=log(SS);
            
            T=(F1-FO)*T/(F1-(T*FO)/(T+1));
            y(:,i)=w./SS;
        end
    end
    T_bp(i)=T;
    
end
format short e;
T_bp;x;y;
figure(1)
hold on
plot(x(1,:),T_bp-273.15,y(1,:),T_bp-273.15),xlabel('x1/y1'),ylabel('Temperature (K)')
hold off
figure(2)
hold on
plot(x(1,:),y(1,:),x(1,:),x(1,:)),xlabel('x1'),ylabel('y1')
hold off


function [output]=valik(x,y,T,P,TC,PC,omg)

R=8.314;     %Gas constant in J/molK
Rg=1.98721;   %Gas constant in cal/molK
A12=1286.0948;
A21=-18.2840;   %parameters in cal/mol

alpha12=0.3011;   %parameters in cal/mol
alpha21=alpha12;  %parameters in cal/mol
%NRTL solved parameters
tao12=(A12)/(Rg*T);
tao21=(A21)/(Rg*T);
G12=exp(-alpha12*tao12);
G21=exp(-alpha21*tao21);
%Activity coefficients based on NRTL
gamma1=exp((x(2)).^2*((tao21*(G21./(x(1)+(x(2)).*G21)).^2)+(tao12*G12./((x(2))+x(1)*G12).^2)));
gamma2=exp((x(1)).^2*((tao12*(G12./(x(2)+(x(1))*G12)).^2)+(tao21*G21./((x(1))+x(2)*G21).^2)));
gamma=[gamma1 gamma2];

%fugacity coefficient of gas phase phiv using Peng-Robinson eos
Tr=T./TC;
f=0.37464 + 1.54226*omg - 0.26992*omg.^2;   %parameter in Peng-Robinson eos
alfa=(1+f.*(1-sqrt(Tr))).^2;                  
aii=(0.45724*(R*TC).^2./(PC)).*alfa;        %vapor vallues of a11 and a22
bii=0.0778*R*TC./PC;                       %vapor values of b11 and b22
kij=-0.0101;
a12=((aii(1)*aii(2)).^0.5)*(1-kij);                %vapor value of a12
a21=a12;
av=(y(1).^2)*aii(1)+2.*y(1).*y(2)*a12+(y(2).^2)*aii(2);  %a for vapor
bv=y(1).*bii(1)+y(2).*bii(2);     %b for vapor
Av=av.*P./(R*T).^2;     %A for vapor
Bv=bv.*P./(R*T);        %B for vapor
Zvroots=roots([1 -(1-Bv) (Av-3*Bv^2-2*Bv) -(Av*Bv-Bv^2-Bv^3)]);

 if isreal(Zvroots)
    Zvp=Zvroots(Zvroots>0);
   	Zv = max(Zvp);  
   else
       disp('imaginary value of Zv')
 end

  %compressibility factors for each component
Fv=Zv-1;
Gv=Zv-Bv;
Hv=Av/((2*sqrt(2))*Bv);
Iv1=(2*(y(1)*aii(1)+y(2)*a21)/av)-bii(1)/bv;
Iv2=(2*(y(1)*a21+y(2)*aii(2))/av)-bii(2)/bv;
Jv=(Zv+2.414*Bv)/(Zv-0.414*Bv);
phiv1=exp((bii(1)/bv)*Fv-log(Gv)-Hv*Iv1*log(Jv));
phiv2=exp((bii(2)/bv)*Fv-log(Gv)-Hv*Iv2*log(Jv));

phiv=[phiv1 phiv2];
%*****************************************************88


%fugacity coefficient of liquid phase phil using Peng-Robinson eos
al=(x(1).^2)*aii(1)+2.*x(1).*x(2)*a12+(x(2).^2)*aii(2);  %a for liquid
bl=x(1).*bii(1)+x(2).*bii(2);     %b for liquid
Al=al.*P./(R*T).^2;     %A for liquid
Bl=bl.*P./(R*T);        %B for liquid
Zlroots=roots([1 -(1-Bl) (Al-3*Bl^2-2*Bl) -(Al*Bl-Bl^2-Bl^3)]);

 if isreal(Zlroots)
    Zlp=Zlroots(Zlroots>0);
   	Zl = min(Zlp);  
   else
       disp('alfmaginary value of Zl')
 end

  %compressibility factors for each component
Fl=Zl-1;
Gl=Zl-Bl;
Hl=Al/((2*sqrt(2))*Bl);
Il1=(2*(x(1)*aii(1)+x(2)*a21)/al)-bii(1)/bl;
Il2=(2*(x(1)*a21+x(2)*aii(2))/al)-bii(2)/bl;
Jl=(Zl+2.414*Bl)/(Zl-0.414*Bl);
phil1=exp((bii(1)/bl)*Fl-log(Gl)-Hl*Il1*log(Jl));
phil2=exp((bii(2)/bl)*Fl-log(Gl)-Hl*Il2*log(Jl));

phil=[phil1 phil2];
Vl=Zl*R*T/P;
%Pisat using the antoine equation
TdegC=T-273.15;   %Temperature in degrees Celcius
A1=6.85146;
B1=1206.47;
C1=223.136;      %Antoine constant for cylcohexane
A2=7.36168;
B2=1180.93;
C2=180.476;      %Antoine constant for tert-butanol
PmmHg1=10^(A1-(B1/(TdegC+C1)));  %saturation pressure in mmHg
P1sat=(PmmHg1/760)*101325;   %saturation Pressure in Pa
PmmHg2=10^(A2-(B2/(TdegC+C2)));  %saturation pressure in mmHg
P2sat=(PmmHg2/760)*101325;   %saturation Pressure in Pa
Pisat=[P1sat P2sat];

%fugacity of pure liquid fi0
fio=phil.*Pisat.*exp((Vl/(R*T))*(P-Pisat));

%calculate K=yi/xi

output=(gamma.*(fio))./(phiv.*Pisat);

end