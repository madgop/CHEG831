clear all
close all
%% Constants and Function Definition
% hours = y; % Day cycle
% time = t;
% v_s= cos(2 * 3.1416 / hours * time) * 0.3785 + 0.3785;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;      % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hvr^-1; first-order rate constant for PER synthesis
v_d=0.95;  % um/hr; maximum degradation rate of biphosphorylated PER (P2)
k_1=1.9;    % hr^-1; first-order rate constant for P2 transport into nucleus
k_2=1.3;    % hr^-1; first-order rate constant for PN transport into cytosol
K_I=1;      % um; threshold constant for repression
K_d=0.2;    % um; Michaelis constant for P2 degradation 
K_14=2;     % um; Michaelis constants for kinases and phosphatases involved 
            %     in reverse phosphorylation of P0 into P1 and P1 into P2
V_1=3.2;    % um/hr; maximum rate for kinases and phosphatases
V_2=1.58;   % um/hr
V_3=5;      % um/hr
V_4=2.5;    % um/hr
n = 4;      % degree of cooperativity

%f is an anonymous function for the concentrations as defined in the paper
f = @(P,y, t) [(cos(2 * 3.1416 / y * t) * 0.3785 + 0.3785)*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
    k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
    V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
    V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
    k_1*P(4)-k_2*P(5)];
%jac is an anonymous Jacobian for the concentrations as defined in the paper
jac = @(PConc,hours,time) [-v_m.*K_m./((K_m+PConc(1)).^2) 0 0 0 -(cos(2 * 3.1416 / hours * time) * 0.37 + 0.3785).*K_I^n*n*(PConc(5)).^(n-1)./((K_I.^n+PConc(5).^n).^2);
    k_s -(V_1*K_14)./((K_14+PConc(2)).^2) V_2*K_14./((K_14+PConc(3)).^2) 0 0;
    0 (V_1.*K_14)/((K_14+PConc(2)).^2) (-V_2.*K_14)/((K_14+PConc(3)).^2) - (V_3.*K_14)./((K_14+PConc(3)).^2) (V_4.*K_14)./((K_14+PConc(4)).^2) 0;
    0 0 V_3.*K_14./((K_14+PConc(3)).^2) (-V_4.*K_14./((K_14+PConc(4)).^2)-k_1 - v_d.*K_d./((K_d+PConc(4)).^2)) k_2;
    0 0 0 k_1 -k_2];
%this is the df/dp vector that we solve for in parametric continuation
dfdhours = @(P_conc, y, t) [-(2*3.1516/y)*sin(2 * 3.1416 / y * t)*K_I^n/(P_conc(5)^n+K_I^n);0;0;0;0];

%% Bifurcation Algorithm

time = 308;

cycle = [44:-5:4];
k = 1;
i = 1; 
[t,P_Conc] = ode45(@(t,P)getC(t,P,cycle(i)),[0,time],[0.6;0.5;1.8;0.65;1.2]); 

P_Conc = P_Conc(t<24, :);
t = t(t<24);

psis = cell(length(cycle),1);
ts = cell(length(cycle),1);
Ps = cell(length(cycle),1);
mdromy = cell(length(cycle),1);

init = [0.6;0.5;1.8;0.65;1.2];

counters = [];

for j = 1:length(cycle)
    poincare = cycle(j);
    [t,P_Conc] = ode45(@(t,P)getC(t,P,poincare),[0,time], init);
    
    Ps{j} = P_Conc;
    ts{j} = t;
    
    stopping = t(end);
    a = 0.0;
    counter = 1;
    while (a + poincare) < stopping
        
        [T,index] = (min(abs(t - a - poincare)));
        mndrmy = (P_Conc(index,:)./P_Conc(counter,:)).';
        mdromy{j} = [mdromy{j},mndrmy];
        
        counter = counter +1;
        a = t(counter);

    end
    
    counters=[counters, counter-1];
    
    [T,index] = (min(abs(t - (fix(time/poincare)-1) * poincare)));
    init = Ps{j}(index,:).';
    

end

figure(1)
hold on

for j = 1:length(cycle)
    step = 1/length(cycle);
    colorspec = [1-step*j, 0.5*step*j,1 ];
    p1 = plot(t(1:(counters(j))), sum(mdromy{j},1), 'LineWidth', 8-8/(length(cycle)+2)*j);
    p1.Color(4) = 0.3;
end
ylabel("Trace of Monodromy Matrix")
xlabel("Time(hours)")
title("Trace of Monodromy Matrix (one cycle ahead) vs Time (308 hrs)")
legendCell = cellstr(num2str(cycle', 'Cycle (hrs) = %-d'));
legend(legendCell, 'location', 'southoutside');

ylim([3,10])
hold off

figure(2)
hold on

for j = 1:length(cycle)
    step = 1/length(cycle);
    colorspec = [1-step*j, 0.5*step*j,1 ];
    x = (cos(2*pi / cycle(j) * ts{j})*0.3875)+0.3785;
    y = sum(Ps{j},2);
    p2 = plot(x, y, 'LineWidth', 4);
    p2.Color(4) = 0.5;
end
legendCell = cellstr(num2str(cycle', 'Cycle (hrs) = %-d'));
legend(legendCell, 'location', 'southoutside');
ylabel("Total mRNA(M)")
xlabel("v_s (cyclic)")
title("Sustained Oscillations between Total Protein and per mRNA (308 hrs)")

hold off

figure(3)
hold on

for j = 1:length(cycle)
    poincare = cycle(j);
    
    stopping = ts{j}(end);
    a = 0.0;
    counter = 1;
    
    ptot = sum(Ps{j},2);
    y = [];
    x = [];
    while (a + poincare) < stopping
        
        [T,index] = (min(abs(ts{j} - poincare*(counter-1))));
        y = [y, ptot(index)];
        x = [x,(cos(2*pi / cycle(j) * ts{j}(index))*0.3785)+0.3785];
        counter = counter +1;
        a = ts{j}(index);

    end

    p3 = scatter(x, y, 50, 'filled');
    alpha(0.2)
end
legendCell = cellstr(num2str(cycle', 'Cycle (hrs) = %-d'));
legend(legendCell, 'location', 'southoutside');
ylabel("Total mRNA (M)")
xlabel("v_s (cyclic)")
title("Poincare Diagram: Total Protein per mRNA vs v_s (308 hrs)")

hold off

% cycle_stable = cycle(1:hopfBif(end));
% mConc_stable = P_array(1,(1:hopfBif(end)));
% cycle_unstable = cycle(hopfBif(end):end);
% mConc_unstable = P_array(1,(hopfBif(end):end));
% pTot_stable = P_t(1,(1:hopfBif(end)));
% pTot_unstable = P_t(hopfBif(end):end);
% 
% v_s = [0.65:0.005:3];
% vs_hopf = v_s;
% for i=1:length(v_s)
%     subM=[]; %empty vector that extract the values of M for each iteration around v_s
%     subP=[]; %empty vector that extract the values of P_tot for each iteration around v_s
%     [t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]); % get vals
%     subM = P_Conc(:,1); %subset only the mRNA data
%     subM = subM(round(length(subM)*0.9):end); % find the steady state vals by only extracting the last 10% of timepoints
%     max_stable_mrna(i) = max(subM); % find the min mRNA val
%     min_stable_mrna(i) = min(subM); %find the max mRNA val
%     
%     subP = (P_Conc(:,2:5)); %subset only the P0,P1,P2,PN data
%     subP = sum(subP,2); %sum along rows to calculate Pn and store in subP
%     subP = subP(round(length(subP)*0.9):end); % extract only the last 10% of timepoints
%     max_stable_totProt(i) = max(subP); %find max total protein at v_s(i)
%     min_stable_totProt(i) = min(subP); %find min total protein at v_s(i)
% end
% 
% %% Plot mRNA vs v_s
% figure(1)
% hold on
% plot(vs_stable, mConc_stable, '-')
% plot(cycle_unstable, mConc_unstable, '--')
% plot(vs_hopf,max_stable_mrna,'-')
% plot(vs_hopf,min_stable_mrna,'-')
% legend({'Stable SS','Unstable', 'Oscillatory Stable', 'Oscillatory Stable'},'Location','northwest')
% ylabel('mRNA Concentration')
% xlabel('v_s')
% title('mRNA Bifurcation Diagram for v_s')
% hold off
% 
% %% Plot Total Protein vs v_s
% figure(2)
% hold on
% plot(vs_stable, pTot_stable, '-')
% plot(cycle_unstable, pTot_unstable, '--')
% plot(vs_hopf,max_stable_totProt,'-')
% plot(vs_hopf,min_stable_totProt,'-')
% legend({'Stable SS','Unstable', 'Oscillatory Stable', 'Oscillatory Stable'},'Location','northwest')
% ylabel('Total Protein Concentration')
% xlabel('v_s')
% title('Total Protein Bifurcation Diagram for v_s')
% hold off

function dP = Jacobian(PConc,cycle,t)
% Function to calculate Jacobian of system
% Inputs: vector of protein concentrations
% (necessarily at a specific time point)
% and y (v_s parameter)
% PConc vector is defined as follows, from original code: 
% PConc(1): Cytosolic concentration, M
% PConc(2): Unphosphorylated PER, P0
% PConc(3): Monophosphorylated PER, P1
% PConc(4): Biphosphorylated PER, P2
% PConc(5): Nuclear PER, PN

PConc = PConc(:)';

% v_s=cos(2 * 3.1416 / cycle * t);   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;      % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d=0.95;  % um/hr; maximum degradation rate of biphosphorylated PER (P2)
k_1=1.9;    % hr^-1; first-order rate constant for P2 transport into nucleus
k_2=1.3;    % hr^-1; first-order rate constant for PN transport into cytosol
K_I=1;      % um; threshold constant for repression
K_d=0.2;    % um; Michaelis constant for P2 degradation 
K_14=2;     % um; Michaelis constants for kinases and phosphatases involved 
            %     in reverse phosphorylation of P0 into P1 and P1 into P2
V_1=3.2;    % um/hr; maximum rate for kinases and phosphatases
V_2=1.58;   % um/hr
V_3=5;      % um/hr
V_4=2.5;    % um/hr
n=4;        % [-]; degree of cooperativity

dP = zeros(5,5);
    for i = 1:5
        if i == 1
            dP(i,1) = -v_m*K_m/((K_m+PConc(1))^2);
            dP(i,5) = -((cos(2 * 3.1416 / cycle * t) * 0.3785 + 0.3785)*K_I^n)*n*(PConc(5))^(n-1)/((K_I^n+PConc(5)^n)^2);
        elseif i == 2
            dP(i,1) = k_s;
            dP(i,2) = -(V_1*K_14)/((K_14+PConc(2))^2);
            dP(i,3) = V_2*K_14/((K_14+PConc(3))^2);
        elseif i == 3
            dP(i,2) = (V_1*K_14)/((K_14+PConc(2))^2);
            dP(i,3) = -(V_2*K_14)/((K_14+PConc(3))^2) - (V_3*K_14)/((K_14+PConc(3))^2);
            dP(i,4) = (V_4*K_14)/((K_14+PConc(4))^2);
        elseif i == 4
            dP(i,3) = V_3*K_14/((K_14+PConc(3))^2);
            dP(i,4) = -V_4*K_14/((K_14+PConc(4))^2) - k_1 - v_d*K_d/((K_d+PConc(4))^2);
            dP(i,5) = k_2;
        else
            dP(i,4) = k_1;
            dP(i,5) = -k_2;
        end
    end
end

function P_Conc=getC(t,P,y)
% Protein Function
% Outputs: Concentrations of 5 proteins of interest
% Inputs: Time, concentration, PER maximum degradation rate

v_s= 0.3785* cos(2 * pi / y * t) + 0.3785   ;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;   % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d= 0.95;   % um/hr; maximum degradation rate of biphosphorylated PER (P2)
k_1=1.9;    % hr^-1; first-order rate constant for P2 transport into nucleus
k_2=1.3;    % hr^-1; first-order rate constant for PN transport into cytosol
K_I=1;      % um; threshold constant for repression
K_d=0.2;    % um; Michaelis constant for P2 degradation 
K_14=2;     % um; Michaelis constants for kinases and phosphatases involved 
            %     in reverse phosphorylation of P0 into P1 and P1 into P2
V_1=3.2;    % um/hr; maximum rate for kinases and phosphatases
V_2=1.58;   % um/hr
V_3=5;      % um/hr
V_4=2.5;    % um/hr
n=4;        % [-]; degree of cooperativity

% Equations 1a through 1e of paper
P_Conc=zeros(5,1);

% Cytosolic concentration, M
P_Conc(1)=v_s*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
% Unphosphorylated PER, P0
P_Conc(2)=k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
% Monophosphorylated PER, P1
P_Conc(3)=V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
% Biphosphorylated PER, P2
P_Conc(4)=V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
% Nuclear PER, PN
P_Conc(5)=k_1*P(4)-k_2*P(5);
%Total Protein
end

function [x_sol iterations] = Newton(fun,jacobian,x_init,N_max,tolerance)
% Newton implements Newton's iterative method for the solution of fun(x)=0

% INPUT:
% fun: user supplied function of x, the function that the root(s) of which are sought
% jacobian: user supplied function of x, the derivative of fun with respect to x
% x: is the independent variable
% x_init: is the initial guess for the independent variable
% N_max: is the maximum number of iterations (10 recommended)
% tolerance: is the absolute tolerance for the solution

% OUTPUT:
% x_sol: is the solution (if converged) last guess (if not converged) 
% iterations: number of iterations (k) required for convergence within given tolerance;
% is equal to -(N_max+1) if no solution is found within N_max iterations;
% it takes a negative value -k if the procedure diverged after k iterations

iterations = 0;

x_sol = x_init;

diff = jacobian(x_sol)\fun(x_sol);
diff_mag_old = norm(diff);
x_sol = x_sol - diff;

    for k=1:N_max
        diff = jacobian(x_sol)\fun(x_sol);
        diff_mag = norm(diff);
            if(diff_mag < tolerance)
                iterations = k;
                x_sol = x_sol - diff;
                break
            elseif(diff_mag > 10*diff_mag_old)
                iterations = -k;
                disp('Warning! Algorithm diverged!')
                break
            end
        diff_mag_old=diff_mag;
        x_sol = x_sol - diff;
    end
    if(iterations == 0) 
        iterations = -(N_max+1);
        disp('Warning! There is no convergence within N__max iterations')
    end
end