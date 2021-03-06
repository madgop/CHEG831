clear all
%% Constants and Function Definition
%v_s=y;   % um/hr; accumulation rate of per mRNA in cytosol
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
f = @(P,y) [y*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
    k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
    V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
    V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
    k_1*P(4)-k_2*P(5)];
%jac is an anonymous Jacobian for the concentrations as defined in the paper
jac = @(PConc,y) [-v_m.*K_m./((K_m+PConc(1)).^2) 0 0 0 -(y.*K_I^n)*n*(PConc(5)).^(n-1)./((K_I.^n+PConc(5).^n).^2);
    k_s -(V_1*K_14)./((K_14+PConc(2)).^2) V_2*K_14./((K_14+PConc(3)).^2) 0 0;
    0 (V_1.*K_14)/((K_14+PConc(2)).^2) (-V_2.*K_14)/((K_14+PConc(3)).^2) - (V_3.*K_14)./((K_14+PConc(3)).^2) (V_4.*K_14)./((K_14+PConc(4)).^2) 0;
    0 0 V_3.*K_14./((K_14+PConc(3)).^2) (-V_4.*K_14./((K_14+PConc(4)).^2)-k_1 - v_d.*K_d./((K_d+PConc(4)).^2)) k_2;
    0 0 0 k_1 -k_2];
%this is the df/dp vector that we solve for in parametric continuation
dfdv_s = @(P_conc,v_s) [K_I^n/(P_conc(5)^n+K_I^n);0;0;0;0];

%% Bifurcation Algorithm
v_s = [0.01:0.0001:3];
k = 1;
i = 1;
[t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]);
P_array(:,1) = P_Conc(end,:);
for j=2:length(v_s)
    Jacob = Jacobian(P_array(:,j-1),v_s(j-1));
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
    delF = dfdv_s(P_array(:,j-1),v_s(j-1));
    dzdp = -inv(Jacob)*delF;
    P_array(:,j) = P_array(:,j-1)+dzdp*(v_s(j)-v_s(j-1));
    Jacob1 = det(Jacobian(P_array(:,j),v_s(j)));
    [xsol, iter] = Newton(@(P,y) f(P,v_s(j)), @(PConc,y) jac(PConc,v_s(j)), P_array(:,j), 10000, 1e-5);
    P_array(:,j) = xsol;
    k = k + 1;
end

[row,col]=size(P_array);
for j=1:col
    P_t(j)=(sum(P_array(2:5,j)));  % Total PER protein, eq. 2 in paper
end
vs_stable = v_s(1:hopfBif(end));
mConc_stable = P_array(1,(1:hopfBif(end)));
vs_unstable = v_s(hopfBif(end):end);
mConc_unstable = P_array(1,(hopfBif(end):end));
pTot_stable = P_t(1,(1:hopfBif(end)));
pTot_unstable = P_t(hopfBif(end):end);

v_s = [0.65:0.005:3];
vs_hopf = v_s;
for i=1:length(v_s)
    subM=[]; %empty vector that extract the values of M for each iteration around v_s
    subP=[]; %empty vector that extract the values of P_tot for each iteration around v_s
    [t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]); % get vals
    subM = P_Conc(:,1); %subset only the mRNA data
    subM = subM(round(length(subM)*0.9):end); % find the steady state vals by only extracting the last 10% of timepoints
    max_stable_mrna(i) = max(subM); % find the min mRNA val
    min_stable_mrna(i) = min(subM); %find the max mRNA val
    
    subP = (P_Conc(:,2:5)); %subset only the P0,P1,P2,PN data
    subP = sum(subP,2); %sum along rows to calculate Pn and store in subP
    subP = subP(round(length(subP)*0.9):end); % extract only the last 10% of timepoints
    max_stable_totProt(i) = max(subP); %find max total protein at v_s(i)
    min_stable_totProt(i) = min(subP); %find min total protein at v_s(i)
end

%% Plot mRNA vs v_s
figure(1)
hold on
plot(vs_stable, mConc_stable, '-')
plot(vs_unstable, mConc_unstable, '--')
plot(vs_hopf,max_stable_mrna,'-')
plot(vs_hopf,min_stable_mrna,'-')
legend('Stable SS','Unstable', 'Oscillatory Stable', 'Oscillatory Stable')
ylabel('mRNA Concentration')
xlabel('v_s')
title('mRNA Bifurcation Diagram for v_s')
hold off

%% Plot Total Protein vs v_s
figure(2)
hold on
plot(vs_stable, pTot_stable, '-')
plot(vs_unstable, pTot_unstable, '--')
plot(vs_hopf,max_stable_totProt,'-')
plot(vs_hopf,min_stable_totProt,'-')
legend('Stable SS','Unstable', 'Oscillatory Stable', 'Oscillatory Stable')
ylabel('Total Protein Concentration')
xlabel('v_s')
title('Total Protein Bifurcation Diagram for v_s')
hold off