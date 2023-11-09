clear
close all
clc

%% Simualtion part: Anodic current from high field model
%Define high feld parameters
A = 1.2E-8; % cm^2, surface area
D = 9E-5; % cm^2/s, diffusion coefficient of O2
i0 = 6.5E-12; % A/cm^2, high field constant
beta = 3.2E-6; % cm/V, high field cosntant
E0 = -1.75; % V vs. SCE, flat band potential
d0 = 2.73E-7; % cm, initial oxide thickness
M = 102; % g/mol, molar mass of Al2O3
rho = 3; % g/cm^3, density of Al2O3
R = 1E-12; % mol/[cm^2*s], dissolution rate of Al2O3
K = 1.7E6; % adjusting constant of pH effect
F = 96485; % C/mol, Faraday's constant
initial_c = 1E-10; % mol/cm^3, initial concentration of H+ in pH=7 solution
z = 6; % number of electrons transferred to produce 1 mol of Al2O3
delta_t = 0.01; % s
i_c0 = 1E-12; % A, cathodic exchange current
beta_c = -0.75; % V, cathodic Tafel slope
E_c = -0.4; % V, OCP or cathodic equilibrium potential


% Cathodic current
% 200 mV/s
scan_rate1 = 0.2; % V/s
t1 = [0.01:0.01:2/scan_rate1];
E1 = scan_rate1*t1-1;

for k=1
    d1(k) = d0;
    i1(k) = i0*exp(beta*(E1(k)-E0)/d1(k))*A;
    V_diff1(k) = sqrt(2*D*t1(k))*A; % diffusion volume at time t=k
    delta_c1(k) = i1(k)*k*delta_t/(V_diff1(k)*F); % concentation change of H+
    C1(k) = 1E-20/initial_c-initial_c-delta_c1(k); % constant to calculate final H+ concentration
    final_c1(k) = (sqrt(C1(k)^2+4E-20)-C1(k))/2; % final H+ concentration
    epsilon_O1(k) = 1/(1+K*final_c1(k)); % oxide formation efficiency
    pH1(k) = -log10(final_c1(k)*1000); % pH at time t=k
end
    
for k=[2:200/scan_rate1]
    d1(k) = d1(k-1)+epsilon_O1(k-1)*M*i1(k-1)*delta_t/(rho*z*F*A)-R*M*delta_t/rho; 
    i1(k) = i0*exp(beta*(E1(k)-E0)/d1(k))*A;
    V_diff1(k) = sqrt(2*D*t1(k))*A;
    delta_c1(k) = 0;
    for n = 1:k
        delta_c1(k) = delta_c1(k)+i1(n)*delta_t/(V_diff1(k+1-n)*F);
    end
    C1(k) = 1E-20/initial_c-initial_c-delta_c1(k);
    final_c1(k) = (sqrt(C1(k)^2+4E-20)-C1(k))/2;
    epsilon_O1(k) = 1/(1+K*final_c1(k));
    pH1(k) = -log10(final_c1(k)*1000);

end

i_c1 = -i_c0*10.^((E1-E_c)/beta_c);

plot(E1,1E12*(i1+i_c1),'--','linewidth',2,'Color',[0, 0.4470, 0.7410])
hold on

% 100 mV/s
scan_rate2 = 0.1; % V/s
t2 = [0.01:0.01:2/scan_rate2];
E2 = scan_rate2*t2-1;

for k=1
    d2(k) = d0;
    i2(k) = i0*exp(beta*(E2(k)-E0)/d2(k))*A;
    V_diff2(k) = sqrt(2*D*t2(k))*A; % diffusion volume at time t=k
    delta_c2(k) = i2(k)*k*delta_t/(V_diff2(k)*F); % concentation change of H+
    C2(k) = 1E-20/initial_c-initial_c-delta_c2(k); % constant to calculate final H+ concentration
    final_c2(k) = (sqrt(C2(k)^2+4E-20)-C2(k))/2; % final H+ concentration
    epsilon_O2(k) = 1/(1+K*final_c2(k)); % oxide formation efficiency
    pH2(k) = -log10(final_c2(k)*1000); % pH at time t=k
end
    
for k=[2:200/scan_rate2]
    d2(k) = d2(k-1)+epsilon_O2(k-1)*M*i2(k-1)*delta_t/(rho*z*F*A)-R*M*delta_t/rho; 
    i2(k) = i0*exp(beta*(E2(k)-E0)/d2(k))*A;
    V_diff2(k) = sqrt(2*D*t2(k))*A;
    delta_c2(k) = 0;
    for n = 1:k
        delta_c2(k) = delta_c2(k)+i2(n)*delta_t/(V_diff2(k+1-n)*F);
    end
    C2(k) = 1E-20/initial_c-initial_c-delta_c2(k);
    final_c2(k) = (sqrt(C2(k)^2+4E-20)-C2(k))/2;
    epsilon_O2(k) = 1/(1+K*final_c2(k));
    pH2(k) = -log10(final_c2(k)*1000);

end

i_c2 = -i_c0*10.^((E2-E_c)/beta_c);

plot(E2,1E12*(i2+i_c2),'--','linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
hold on

% 50 mV/s
scan_rate3 = 0.05; % V/s
t3 = [0.01:0.01:2/scan_rate3];
E3 = scan_rate3*t3-1;

for k=1
    d3(k) = d0;
    i3(k) = i0*exp(beta*(E3(k)-E0)/d3(k))*A;
    V_diff3(k) = sqrt(2*D*t3(k))*A; % diffusion volume at time t=k
    delta_c3(k) = i3(k)*k*delta_t/(V_diff3(k)*F); % concentation change of H+
    C3(k) = 1E-20/initial_c-initial_c-delta_c3(k); % constant to calculate final H+ concentration
    final_c3(k) = (sqrt(C3(k)^2+4E-20)-C3(k))/2; % final H+ concentration
    epsilon_O3(k) = 1/(1+K*final_c3(k)); % oxide formation efficiency
    pH3(k) = -log10(final_c3(k)*1000); % pH at time t=k
end
    
for k=[2:200/scan_rate3]
    d3(k) = d3(k-1)+epsilon_O3(k-1)*M*i3(k-1)*delta_t/(rho*z*F*A)-R*M*delta_t/rho; 
    i3(k) = i0*exp(beta*(E3(k)-E0)/d3(k))*A;
    V_diff3(k) = sqrt(2*D*t3(k))*A;
    delta_c3(k) = 0;
    for n = 1:k
        delta_c3(k) = delta_c3(k)+i3(n)*delta_t/(V_diff3(k+1-n)*F);
    end
    C3(k) = 1E-20/initial_c-initial_c-delta_c3(k);
    final_c3(k) = (sqrt(C3(k)^2+4E-20)-C3(k))/2;
    epsilon_O3(k) = 1/(1+K*final_c3(k));
    pH3(k) = -log10(final_c3(k)*1000);

end

i_c3 = -i_c0*10.^((E3-E_c)/beta_c);

plot(E3,1E12*(i3+i_c3),'--','linewidth',2,'Color',[0.9290, 0.6940, 0.125])
hold on

% 25 mV/s
scan_rate4 = 0.025; % V/s
t4 = [0.01:0.01:2/scan_rate4];
E4 = scan_rate4*t4-1;

for k=1
    d4(k) = d0;
    i4(k) = i0*exp(beta*(E4(k)-E0)/d4(k))*A;
    V_diff4(k) = sqrt(2*D*t4(k))*A; % diffusion volume at time t=k
    delta_c4(k) = i4(k)*k*delta_t/(V_diff4(k)*F); % concentation change of H+
    C4(k) = 1E-20/initial_c-initial_c-delta_c4(k); % constant to calculate final H+ concentration
    final_c4(k) = (sqrt(C4(k)^2+4E-20)-C4(k))/2; % final H+ concentration
    epsilon_O4(k) = 1/(1+K*final_c4(k)); % oxide formation efficiency
    pH4(k) = -log10(final_c4(k)*1000); % pH at time t=k
end
    
for k=[2:200/scan_rate4]
    d4(k) = d4(k-1)+epsilon_O4(k-1)*M*i4(k-1)*delta_t/(rho*z*F*A)-R*M*delta_t/rho; 
    i4(k) = i0*exp(beta*(E4(k)-E0)/d4(k))*A;
    V_diff4(k) = sqrt(2*D*t4(k))*A;
    delta_c4(k) = 0;
    for n = 1:k
        delta_c4(k) = delta_c4(k)+i4(n)*delta_t/(V_diff4(k+1-n)*F);
    end
    C4(k) = 1E-20/initial_c-initial_c-delta_c4(k);
    final_c4(k) = (sqrt(C4(k)^2+4E-20)-C4(k))/2;
    epsilon_O4(k) = 1/(1+K*final_c4(k));
    pH4(k) = -log10(final_c4(k)*1000);

end

i_c4 = -i_c0*10.^((E4-E_c)/beta_c);

plot(E4,1E12*(i4+i_c4),'--','linewidth',2,'Color',[0.4940, 0.1840, 0.5560])


xlabel('E vs. SCE (V)')
ylabel('i (pA)')
xlim([-1 1])
ylim([-10 300])
set(gca,'FontSize',18,'linewidth',2)


figure
plot(E1, pH1,'linewidth',2)
hold on;
plot(E2, pH2,'linewidth',2)
hold on;
plot(E3, pH3,'linewidth',2)
hold on;
plot(E4, pH4,'linewidth',2)
xlabel('E vs. SCE (V)')
ylabel('pH')
xlim([-1 1])
set(gca,'FontSize',18,'linewidth',2)


