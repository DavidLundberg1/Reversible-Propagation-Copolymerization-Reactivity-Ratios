function[rxn_conversion,M1_conversion,M2_conversion,c_out] = PopulationBalanceODEs(rates,conv,chain_length)
%Numerically integrates a series of coupled ODEs describing the population
%balance of different dyad species
%Inputs: rates, a vector containing [f1_0, M0, k1, k2, ..., k7, k8] 
%f1_0 = initial monomer 1 mole fraction
%M0 = initial monomer concentration (M)
%k1, k2, k3, k4, k5, k6, k7, k8 are rate constants as defined in Scheme 1
%of the manuscript. 


global ode_conv
ode_conv = conv;
comparison_params = rates; %[f1, R1, R2, beta1, gamma1, beta2, gamma2]

M1_0 = comparison_params(1);
M2_0 = 1-comparison_params(1);
init_frac = 1/chain_length;
c0 = zeros(1,(8*chain_length+24));
c0(1) = M1_0*(1-init_frac); %Initiate chains according to chain length
c0(2) = M2_0*(1-init_frac);
c0(3) = M2_0*(init_frac);
c0(4) = M1_0*(init_frac);
c0(5:end) = 1e-18;

options2 = odeset('MaxStep',0.1,'Events',@myEvent,'NonNegative',1);
tspan = linspace(0, 5000,500); %NOTE: This integration proceeds in terms of unit time. Integration bounds may need to be manually adjusted depending on rate constant magnitudes. 
[t_out, c_out] = ode45(@(t,c) odefunODE(t,c,comparison_params,chain_length,ode_conv),tspan,c0,options2);

rxn_conversion = zeros(1,length(t_out));
M1_conversion = zeros(1,length(t_out));
M2_conversion = zeros(1,length(t_out));
for i = 1:length(t_out)
    rxn_conversion(i)= 1-(c_out(i,1)+c_out(i,2))/(M1_0+M2_0); 
    M1_conversion(i)= 1 - c_out(i,1)/M1_0;
    M2_conversion(i)= 1 - c_out(i,2)/M2_0;
end
%Potential to plot
%hold on
% plot(rxn_conversion,M1_conversion,'r.');
% plot(rxn_conversion,M2_conversion,'k.');

%Now compare results to CPE integration...
end


function dMdt = odefunODE(t, c,comparison_params,chain_length,ode_conv)
k1 = 1; %Monomer 1 homo-propagation
k2 = 1*comparison_params(4); %Monomer 1 homo-depropagation 
k3 = 1/comparison_params(2); %Monomer 1 cross-propagation 
k4 = comparison_params(7); %Monomer 2 cross-depropagation 
k5 = 1/comparison_params(3); %Monomer 2 cross-propagation
k6 = 1*comparison_params(5); %Monomer 1 cross-depropagation 
k7 = 1; %Monomer 2 homo-propagation
k8 = comparison_params(6); %monomer 2 homo-depropagation 
M1_0 = comparison_params(1); %Monomer 1 initial concentration
M2_0 = 1-comparison_params(1); %Monomer 2 initial concentration

k11 = k1;
k12 = k3;
k22 = k7;
k21 = k5;
k11r = k2;
k12r = k4;
k22r = k8;
k21r = k6;

%Define concentrtaions of living dyad species. 
active11 = c(5);
active21 = c(6);
active22 = c(7);
active12 = c(8);

for i = 0:(chain_length)
    j = (i*8)+13;
    active11 = active11 + c(j) ;
    active21 = active21 + c(j+1);
    active22 = active22 + c(j+2);
    active12 = active12 + c(j+3);
end

dMdt(1) = -k11*c(1)*(c(4) + active11 + active21) - k21*c(1)* (c(3) + active22 + active12) + k11r*active11 + k21r*active21; %M1
dMdt(2) = -k22*c(2)*(c(3) + active22 + active12) - k12*c(2)*(c(4)+active11 + active21) + k22r*active22 + k12r*active12; %M2

%M*2,1
dMdt(3) = -c(3)*c(2)*k22 -c(3)*c(1)*k21 + c(6)*k21r + c(7)*k22r; %Active M2 end 1 unit
%M*1,1 
dMdt(4) = -c(4)*c(1)*k11 - c(4)*c(2)*k12 + c(5)*k11r + c(8)*k12r;
%c(1) to c(4) are monomer and length 1 active species, respectively.
%Define length 2 active species here
dMdt(5) = +c(4)*c(1)*k11     -c(5)*(k11r + k11*c(1) + k12*c(2))    + c(13)*k11r*(c(9)/(c(9)+c(10))) + c(16)*k12r*(c(9)/(c(9)+c(10))); %11*
dMdt(6) = +c(3)*c(1)*k21     -c(6)*(k21r + k11*c(1) + k12*c(2))    + c(13)*k11r*(c(10)/(c(9)+c(10))) + c(16)*k12r*(c(10)/(c(9)+c(10))); %21*
dMdt(7) = +c(3)*c(2)*k22     -c(7)*(k22r + k21*c(1) + k22*c(2))    + c(14)*k21r*(c(11)/(c(11)+c(12))) + c(15)*k22r*(c(11)/(c(11)+c(12))) ; %22*
dMdt(8) = +c(4)*c(2)*k12     -c(8)*(k12r + k21*c(1) + k22*c(2))    + c(14)*k21r*(c(12)/(c(11)+c(12))) + c(15)*k22r*(c(12)/(c(11)+c(12))) ; %12*
%Define dormant dyads at length 2
dMdt(9) = +c(5)*(k11*c(1) + k12*c(2))  - c(13)*k11r*(c(9)/(c(9)+c(10))) - c(16)*k12r*(c(9)/(c(9)+c(10))); %11
dMdt(10) = +c(6)*(k11*c(1) + k12*c(2)) - c(13)*k11r*(c(10)/(c(9)+c(10))) - c(16)*k12r*(c(10)/(c(9)+c(10))); %21
dMdt(11) = +c(7)*(k21*c(1) + k22*c(2)) - c(14)*k21r*(c(11)/(c(11)+c(12))) - c(15)*k22r*(c(11)/(c(11)+c(12))) ; %22
dMdt(12) = +c(8)*(k21*c(1) + k22*c(2)) - c(14)*k21r*(c(12)/(c(11)+c(12))) - c(15)*k22r*(c(12)/(c(11)+c(12))) ; %12
%Define legnth 3 and beyond active species based on number of chains to
%simulate.  
for i = 1:(chain_length-1)
    %Starting index to define new equations
    j = (i-1)*8 + 13;
    %Next index values for depropagation
    k = j + 8;
%Define active chain ends
dMdt(j)   = c(1)*k11*(c(j-8)+c(j-7))     -c(j)*(k11r + k11*c(1) + k12*c(2))    + c(k)*k11r*(c((j+4))/(c((j+4))+c((j+5)))) + c(k+3)*k12r*(c((j+4))/(c((j+4))+c((j+5)))); %11*
dMdt(j+1) = c(1)*k21*(c(j-6)+c(j-5))     -c(j+1)*(k21r + k11*c(1) + k12*c(2))    + c(k)*k11r*(c((j+5))/(c((j+4))+c((j+5)))) + c(k+3)*k12r*(c((j+5))/(c((j+4))+c((j+5)))); %21*
dMdt(j+2) = c(2)*k22*(c(j-6)+c(j-5))     -c(j+2)*(k22r + k21*c(1) + k22*c(2))    + c(k+1)*k21r*(c(j+6)/(c(j+6)+c(j+7))) + c(k+2)*k22r*(c(j+6)/(c(j+6)+c(j+7))) ; %22*
dMdt(j+3) = c(2)*k12*(c(j-8)+c(j-7))     -c(j+3)*(k12r + k21*c(1) + k22*c(2))    + c(k+1)*k21r*(c(j+7)/(c(j+6)+c(j+7))) + c(k+2)*k22r*(c(j+7)/(c(j+6)+c(j+7))) ; %12*
%Define dormant dyads at length 2
dMdt(j+4) = +c(j)*(k11*c(1) + k12*c(2))  - c(k)*k11r*(c((j+4))/(c((j+4))+c((j+5)))) - c(k+3)*k12r*(c((j+4))/(c((j+4))+c((j+5)))); %11
dMdt(j+5) = +c(j+1)*(k11*c(1) + k12*c(2)) - c(k)*k11r*(c((j+5))/(c((j+4))+c((j+5)))) - c(k+3)*k12r*(c((j+5))/(c((j+4))+c((j+5)))); %21
dMdt(j+6) = +c(j+2)*(k21*c(1) + k22*c(2)) - c(k+1)*k21r*(c(j+6)/(c(j+6)+c(j+7))) - c(k+2)*k22r*(c(j+6)/(c(j+6)+c(j+7))) ; %22
dMdt(j+7) = +c(j+3)*(k21*c(1) + k22*c(2)) - c(k+1)*k21r*(c(j+7)/(c(j+6)+c(j+7))) - c(k+2)*k22r*(c(j+7)/(c(j+6)+c(j+7))) ; %12

end
j = chain_length*8 + 13;
k = j+8;
 
dMdt(k) = 0;
dMdt(k+1) = 0;
dMdt(k+2) = 0;
dMdt(k+3) = 0;

dMdt = dMdt'; %Change dimension of derivative vector.

end

function [value, isterminal, direction] = myEvent(t, c,ode_conv)
global ode_conv
ode_conv;
%Define integration stopping criteria
%**********************************************************************************spec*
value      = ((c(1)+c(2))/1 < (1-ode_conv)) ;
isterminal = 1;   % Stop the integration
direction  = 0;
end
