function[x,M1_out,M2_out] = IzuCPE(rates,conversion,chain_length)
%Numerically integrates the Izu copolymer equation
%Inputs: rates, a vector containing [f1_0, M0, k1, k2, ..., k7, k8] 
%f1_0 = initial monomer 1 mole fraction
%M0 = initial monomer concentration (M)
%k1, k2, k3, k4, k5, k6, k7, k8 are rate constants as defined in Scheme 1
%of the manuscript. 

%Output is a vector of total monomer conversion (x), and individual monomer
%concentrations (M1_out and M2_out, in M). 

global CPE_params
syms a b ep eta 

CPE_params = zeros(7,1);
CPE_params(1,1) = rates(1);
CPE_params(2,1) = 1-rates(1);
CPE_params(3,1) = CPE_params(1,1);
CPE_params(4,1) = CPE_params(2,1);


f1_0 = rates(1);
M1_0 = f1_0*1;
M2_0 = 1-M1_0;

%the bounds of integration in 'tspan' begin assuming 
tspan = linspace((1/chain_length), conversion,250);
options2 = odeset('MaxStep',0.001,'NonNegative',1);
[x, f] = ode23(@(x,f) odefunIzu(x,f,rates,a,b,ep,eta),tspan,f1_0,options2);
x=x';
M1_out = zeros(1,length(x));
M2_out = zeros(1,length(x));

for i=1:length(x)
    M1_out(i) = (M1_0-(1-x(i))*(M1_0+M2_0)*f(i))/M1_0;
    M2_out(i) = (M2_0-(1-x(i))*(M1_0+M2_0)*(1-f(i)))/M2_0;
end

%Option to plot 
%hold on
%plot(x,M1_out,'ko');
%plot(x,M2_out,'ro');

end

function dfdx = odefunIzu(x,f,rates,a,b,ep,eta) 
global CPE_params

options = optimoptions('fsolve','Display','off','MaxFunEvals', 1e6, 'MaxIter', 1e6);
CPE_variables= fsolve(@(y) Izu_algebra(y,rates,f,x),CPE_params(1:4,end),options);

a = CPE_variables(1);
b = CPE_variables(2);
ep = CPE_variables(3);
eta = CPE_variables(4);

CPE_params(1,end+1)=a;
CPE_params(2,end)=b;
CPE_params(3,end)=ep;
CPE_params(4,end)=eta;
CPE_params(5,end)=x;

k1 = 1; %Monomer 1 homo-propagation
k2 = rates(4); %Monomer 1 homo-depropagation 
k3 = 1/rates(2); %Monomer 1 cross-propagation 
k4 = rates(7); %Monomer 2 cross-depropagation 
k5 = 1/rates(3); %Monomer 2 cross-propagation
k6 = rates(5); %Monomer 1 cross-depropagation 
k7 = 1; %Monomer 2 homo-propagation
k8 = rates(6); %monomer 2 homo-depropagation 

M1_0 = rates(1); %Monomer 1 initial concentration
M2_0 = 1-rates(1); %Monomer 2 initial concentration

M1 = (1-x)*(M1_0+M2_0)*f;
M2 = (1-x)*(M1_0+M2_0) - M1;  

An = M1*(a*k1 + b*k5) - a*((1-ep)*k6+ep*k2);
Bn = M2*(b*k7 + a*k3) - b*((1-eta)*k4+eta*k8);
F1 = An/(An+Bn);
CPE_params(6,end)=F1;
CPE_params(7,end)=f;
dfdx = (f - F1)/(1-x);

end

function[F] =  Izu_algebra(y,rates,f,x)

a = y(1);
b = y(2);
ep = y(3);
eta = y(4);

k1 = 1; %Monomer 1 homo-propagation
k2 = rates(4); %Monomer 1 homo-depropagation 
k3 = 1/rates(2); %Monomer 1 cross-propagation 
k4 = rates(7); %Monomer 2 cross-depropagation 
k5 = 1/rates(3); %Monomer 2 cross-propagation
k6 = rates(5); %Monomer 1 cross-depropagation 
k7 = 1; %Monomer 2 homo-propagation
k8 = rates(6); %monomer 2 homo-depropagation 
M1_0 = rates(1); %Monomer 1 initial concentration
M2_0 = 1-rates(1); %Monomer 2 initial concentration

M1 = (1-x)*(M1_0+M2_0)*f;
M2 = (1-x)*(M1_0+M2_0) - M1;   

F(1) = a + b - 1;
F(2) =  a*(k3*M2 + (1-ep)*k6) - b*(k5*M1 + (1-eta)*k4);
F(3) = b*ep*(1-eta)*k4 -a*(ep*(k1*M1+k2+k3*M2)-(k1*M1+ep*ep*k2));
F(4) = a*eta*(1-ep)*k6 - b*(eta*(k7*M2+k8+k5*M1) - (k7*M2 + eta*eta*k8));

end