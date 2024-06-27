function[] = Example_iPrSi8_EndoNB_Fitting()
%This functions is an example of fitting experimental reversible
%copolymerization data using numerical integration of the Izu copolymer
%equation. Data is for an iPrSi8/endoNB copolymerization following the
%experimental procedure in the Supporting Information
%For this example, we have previously measured the equilibrium
%concentration of iPrSi8 and are usign the average fit value as a constant.
%Therefore, only reactivity ratios for each comonomer and the heterodyad
%reversible propagation rate (gamma) are fit parameters. 

global fit_data mole_fraction_irreversible

total_conv = [0 
0.114531351
0.218601784
0.28241979
0.340307054
0.397784524
0.446512912
0.489586046
0.531282343
0.566291897
0.60017465
0.630041463
0.657349136
0.705293393
0.745977988
0.77881917
0.833871424
0.873984959
0.902120023];

nbe_conv = [0 
0.140487462
0.252187528
0.324250027
0.389350816
0.451892226
0.505882449
0.554539308
0.604951268
0.640420721
0.677954775
0.711861035
0.741476794
0.797666784
0.839796066
0.873655303
0.930832537
0.971541533
0.994899832];

iPrSi8_conv_data = [0 
0.08857524
0.185016039
0.240589553
0.291263293
0.343676822
0.387143375
0.424632785
0.457613417
0.492163073
0.522394525
0.548221891
0.573221478
0.612920002
0.652159909
0.683983037
0.73691031
0.776428386
0.809340214];

fit_data = zeros(3,length(total_conv));
fit_data(1,:) = total_conv;
fit_data(2,:) = nbe_conv; %second row monomer 1 (irreversible) conversion
fit_data(3,:) = iPrSi8_conv_data; %third row monomer 2 (reversible) conversion

mole_fraction_irreversible = 0.5; 
hold on
%Plots data points
plot(fit_data(1,:),fit_data(2,:),'o');
plot(fit_data(1,:),fit_data(3,:),'o');

options = optimoptions("lsqnonlin",'Display','iter','MaxIterations',50);
[fitted_parameters,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@(fit_params) comparison_fcn(fit_params),[0.4, 0.08, 0.047],[0.001, 0.001, 0.001],[10,10, 1],options);
ci = nlparci(fitted_parameters,residual,'jacobian',jacobian);

end

function ssr = comparison_fcn(fit_params)
global fit_data mole_fraction_irreversible mole_fraction_irreversible; 
f0 = 0.5;
tspan = [fit_data(1,:)];
options2 = odeset('MaxStep',0.01,'NonNegative',1);
%Run simulation with input parameters
[x, f] = ode45(@(x,f) odefunCPE(x,f,fit_params),tspan,f0,options2);
%calculate total monomer conversion from output of this function. 
CPE_M1_out = zeros(1,length(x));
CPE_M2_out = zeros(1,length(x)); 
M1_0 = 0.1; %Monomer 1 initial concentration
    M2_0 = 0.1; %Monomer 2 initial concentration
for i=1:length(x)
    CPE_M1_out(i) = (M1_0-(1-x(i))*(M1_0+M2_0)*f(i))/M1_0;
    CPE_M2_out(i) = (M2_0-(1-x(i))*(M1_0+M2_0)*(1-f(i)))/M2_0;
end
%calculate ssr value
diff = [(CPE_M1_out-fit_data(2,:)), CPE_M2_out-fit_data(3,:)];
ssr=diff;
plot(x,CPE_M1_out);
plot(x,CPE_M2_out);
end

function dfdx = odefunCPE(x,f,fit_params) 
global fit_data
%solves algebraic equations. 
options = optimoptions('lsqnonlin','Display','off');
CPE_variables= lsqnonlin(@(y) CPE_algebra(y,fit_params,f,x),[0.5,0.5,0.5,0.5],[0,0,0,0],[1,1,1,1],options);
a = CPE_variables(1);
b = CPE_variables(2);
ep = CPE_variables(3);
eta = CPE_variables(4);
k1 = 1; %Monomer 1 homo-propagation
k2 = 0; %Monomer 1 homo-depropagation *zero for this system
k3 = 1/fit_params(1); %Monomer 1 cross-propagation 
k4 = fit_params(3); %Monomer 2 cross-depropagation 
k5 = 1/fit_params(2); %Monomer 2 cross-propagation
k6 = 0; %Monomer 1 cross-depropagation *zero for this system
k7 = 1; %Monomer 2 homo-propagation
k8 =0.047; %monomer 2 homo-depropagation (measured previously)

M1_0 = 0.1; %Monomer 1 initial concentration
M2_0 = 0.1; %Monomer 2 initial concentration

M1 = (1-x)*(M1_0+M2_0)*f;
M2 = (1-x)*(M1_0+M2_0) - M1;  

%Calculate F1 value now
An = M1*(a*k1 + b*k5) - a*((1-ep)*k6+ep*k2);
Bn = M2*(b*k7 + a*k3) - b*((1-eta)*k4+eta*k8);
F1 = An/(An+Bn);
%Define value of df/dx given the above
dfdx = (f - F1)/(1-x);
end

function[F] =  CPE_algebra(y,fit_params,f,x)
a = y(1);
b = y(2);
ep = y(3);
eta = y(4);

k1 = 1; %Monomer 1 homo-propagation
k2 = 0; %Monomer 1 homo-depropagation *zero for this system
k3 = 1/fit_params(1); %Monomer 1 cross-propagation 
k4 = fit_params(3); %Monomer 2 cross-depropagation 
k5 = 1/fit_params(2); %Monomer 2 cross-propagation
k6 = 0; %Monomer 1 cross-depropagation *zero for this system
k7 = 1; %Monomer 2 homo-propagation
k8 =0.047; %monomer 2 homo-depropagation (measured previously)

M1_0 = 0.1; %Monomer 1 initial concentration
M2_0 = 0.1; %Monomer 2 initial concentration

M1 = (1-x)*(M1_0+M2_0)*f;
M2 = (1-x)*(M1_0+M2_0) - M1;   

F(1) = a + b - 1;
F(2) =  a*(k3*M2 + (1-ep)*k6) - b*(k5*M1 + (1-eta)*k4);
F(3) = b*ep*(1-eta)*k4 -a*(ep*(k1*M1+k2+k3*M2)-(k1*M1+ep*ep*k2));
F(4) = a*eta*(1-ep)*k6 - b*(eta*(k7*M2+k8+k5*M1) - (k7*M2 + eta*eta*k8));

end
