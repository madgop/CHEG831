clear all
%% Constants and Function Definition
v_s=0.76;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;      % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hvr^-1; first-order rate constant for PER synthesis
%v_d=y;  % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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
f = @(P,y) [v_s*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
    k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
    V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
    V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-y*(P(4)/(K_d+P(4)));
    k_1*P(4)-k_2*P(5)];
%jac is an anonymous Jacobian for the concentrations as defined in the paper
jac = @(PConc,y) [-v_m.*K_m./((K_m+PConc(1)).^2) 0 0 0 -(v_s.*K_I^n)*n*(PConc(5)).^(n-1)./((K_I.^n+PConc(5).^n).^2);
    k_s -(V_1*K_14)./((K_14+PConc(2)).^2) V_2*K_14./((K_14+PConc(3)).^2) 0 0;
    0 (V_1.*K_14)/((K_14+PConc(2)).^2) (-V_2.*K_14)/((K_14+PConc(3)).^2) - (V_3.*K_14)./((K_14+PConc(3)).^2) (V_4.*K_14)./((K_14+PConc(4)).^2) 0;
    0 0 V_3.*K_14./((K_14+PConc(3)).^2) (-V_4.*K_14./((K_14+PConc(4)).^2)-k_1 - y.*K_d./((K_d+PConc(4)).^2)) k_2;
    0 0 0 k_1 -k_2];
%this is the df/dp vector that we solve for in parametric continuation
dfdv_d = @(P_conc,v_d) [0;0;0;(-P_conc(4)./(K_d+P_conc(4)));0];

%% Bifurcation Algorithm
v_d = [0.01:0.0001:3];
k = 1;
i = 1;
[t,P_Conc] = ode23s(@(t,P)getC(t,P,v_d(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]);
P_array(:,1) = P_Conc(end,:);
for j=2:length(v_d)
    Jacob = Jacobian(P_array(:,j-1),v_d(j-1));
    eigJacob(:,j) = eig(Jacob);
    detJacob = abs(det(Jacob));
    if detJacob < 1e-4
        tp(k) = j;
        k = k+1;
        break
    end
    if abs(real(eigJacob(1,j)))<1e-4 || abs(real(eigJacob(2,j)))<1e-4 || abs(real(eigJacob(3,j)))<1e-4 || abs(real(eigJacob(4,j)))<1e-4 || abs(real(eigJacob(5,j)))<1e-4
        hopfBif(i) = j;
        i = i + 1;
    end
    delF = dfdv_d(P_array(:,j-1));
    dzdp = -inv(Jacob)*delF;
    P_array(:,j) = P_array(:,j-1)+dzdp*(v_d(j)-v_d(j-1));
    Jacob1 = det(Jacobian(P_array(:,j),v_d(j)));
    [xsol, iter] = Newton(@(P,y) f(P,v_d(j)), @(PConc,y) jac(PConc,v_d(j)), P_array(:,j), 10000, 1e-5);
    P_array(:,j) = xsol;
    k = k + 1;
end

[row,col]=size(P_array);
for j=1:col
    P_t(j)=(sum(P_array(2:5,j)));  % Total PER protein, eq. 2 in paper
end

vd_stable = v_d(1:hopfBif(1));
vd_unstable = v_d(hopfBif(1):hopfBif(end));
vd_stable2 = v_d(hopfBif(end):end);
mConc_stable = P_array(1,(1:hopfBif(1)));
mConc_unstable = P_array(1,(hopfBif(1):hopfBif(end)));
mConc_stable2 = P_array(1,(hopfBif(end)):end);
pTot_stable = P_t(1,(1:hopfBif(1)));
pTot_unstable = P_t(hopfBif(1):hopfBif(end));
pTot_stable2 = P_t(1,(hopfBif(end)):end);

vd_hopf = [v_d(hopfBif(1)):0.01:v_d(hopfBif(end))];
for i=1:length(vd_hopf)
    subM=[]; %empty vector that extract the values of M for each iteration around v_d
    subP=[]; %empty vector that extract the values of P_tot for each iteration around v_d
    [t,P_Conc] = ode23s(@(t,P)getC(t,P,vd_hopf(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]); % get vals
    subM = P_Conc(:,1); %subset only the mRNA data
    subM = subM(round(length(subM)*0.9):end); % find the steady state vals by only extracting the last 10% of timepoints
    max_stable_mrna(i) = max(subM); % find the min mRNA val
    min_stable_mrna(i) = min(subM); %find the max mRNA val
    
    subP = (P_Conc(:,2:5)); %subset only the P0,P1,P2,PN data
    subP = sum(subP,2); %sum along rows to calculate Pn and store in subP
    subP = subP(round(length(subP)*0.9):end); % extract only the last 10% of timepoints
    max_stable_totProt(i) = max(subP); %find max total protein at v_d(i)
    min_stable_totProt(i) = min(subP); %find min total protein at v_d(i)
end

%% Plot mRNA vs v_d
figure(1)
hold on
plot(vd_stable, mConc_stable, '-')
plot(vd_unstable, mConc_unstable, '--')
plot(vd_stable2, mConc_stable2, '-')
plot(vd_hopf,max_stable_mrna,'-')
plot(vd_hopf,min_stable_mrna,'-')
legend({'Stable SS','Unstable', 'Stable SS', 'Oscillatory Stable', 'Oscillatory Stable'},'Location','northwest')
ylabel('mRNA Concentration')
xlabel('v_d')
title('mRNA Bifurcation Diagram for v_d')
hold off

%% Plot Total Protein vs v_d
figure(2)
hold on
plot(vd_stable, pTot_stable, '-')
plot(vd_unstable, pTot_unstable, '--')
plot(vd_stable2, pTot_stable2, '-')
plot(vd_hopf,max_stable_totProt,'-')
plot(vd_hopf,min_stable_totProt,'-')
legend({'Stable SS','Unstable', 'Stable SS', 'Oscillatory Stable', 'Oscillatory Stable'},'Location','northwest')
ylabel('Total Protein Concentration')
xlabel('v_d')
title('Total Protein Bifurcation Diagram for v_d')
hold off

function dP = Jacobian(PConc,y)
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

v_s=0.76;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;      % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d=y;  % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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
            dP(i,5) = -(v_s*K_I^n)*n*(PConc(5))^(n-1)/((K_I^n+PConc(5)^n)^2);
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
            dP(i,4) = -V_4*K_14/((K_14+PConc(4))^2) - k_1 - y*K_d/((K_d+PConc(4))^2);
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

v_s=0.76   ;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;   % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d= y;   % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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