% CHEG 831 Project

clc
clear all
close all
% Initial concentrations of the proteins evaluated in equations 1a through 1e
init_1=[0.1 0.25 0.25 0.25 0.25;...    % 
                0.7 0.5 0.5 0.5 0.5;...        %
1.9 0.8 0.8 0.8 0.8];          %
hours = 24;
%
% Figure 1 (Figure 2 of paper)
% 
% Format for ode solver is ode45(function, time_span, initial_cond)
[t,P_Conc] = ode45(@(t,P)getC(t,P,0.95, hours),[0,100],[0.6;0.5;1.8;0.65;1.2]);
[rows,columns]=size(P_Conc);
for j=1:rows
    P_t(j,1)=(sum(P_Conc(j,2:5)));  % Total PER protein, eq. 2 in paper
end

P_Conc(:,7)=P_t;    
figure(1)
hold on

for i = 1:5
    plot(t, P_Conc(:,i))
end

plot(t, cos(2 * 3.1416 / hours * t)* 0.3785 + 0.3785);
plot(t,P_Conc(:,7)) % Plotting all protein concentrations vs time




legend("M","Po","P1","P2","PN", 'v_s', "Pt")
ylabel("Per forms of M")
xlabel("Time (hrs)")

title("Sustained oscillations of the model equations (24 hr cycle)")
hold off
%%
% Figure 2 (Figure 3 of paper)
totPhours12 = cell(3,1);
vss = cell(3,1);
t12s = cell(3,1);
tss = cell(3,1);
pt2 = cell(3,1);
P_t2=[]     % total PER protein 
P_M=[]      % PER mRNA
for k=1:(size(init_1))[2] 
% size(init_1)[2] is number of rows in init_1 (dimension 2)
P_Conc_2=[] % New vector of protein concentrations
[t2,P_Conc_2] = ode45(@(t,P)getC(t, P, 0.95, hours),[0,120],init_1(k,:));
[rows_2,columns]=size(P_Conc_2);
size(P_Conc_2)
    for j=1:rows_2
        P_t2(j,k)=(sum(P_Conc_2(j,2:5)));   % Total PER protein
        P_M(j,k)=P_Conc_2(j,1);             % PER mRNA
    end
figure(2)

hold on
% 
% for i = 3:5
%    p1 = plot((cos(2 * 3.1416 / hours * t2)* 0.3785 + 0.3785), P_Conc_2(:,i), 'LineWidth', 5);
% end

color = [0.9 - 0.4*(k-1), 0.4*(k-1), 0.9];

p1 = plot((cos(2 * 3.1416 / hours * t2)* 0.3785 + 0.3785),P_t2(1:length(P_Conc_2),k), 'LineWidth', (10 - 3*k), 'Color', color);
p1.Color(4) = 0.25;
ylabel("Total mRNA(M)")
xlabel("v_s (cyclic)")
title("Sustained Oscillations between Total Protein and per mRNA")

poincare = hours;
tss{k} = t2;
pt2{k} = P_t2(1:length(P_Conc_2), k);
    for i = 0:(fix(t2(end)/poincare))
    
        [T,index] = (min(abs(t2 - poincare * i)));
        vss{k}(i+1) = cos(2 * 3.1416 / poincare * t2(index))* 0.3785 + 0.3785;
        totPhours12{k}(i+1) = P_t2(index,k);
        t12s{k}(i+1) = t2(index);
        
    end

end

for k=1:(size(init_1))[2] 
color = [0.9 - 0.4*(k-1), 0.4*(k-1), 0.9];
    p2 = scatter(vss{k}, totPhours12{k}, 60, color, 'filled');
    alpha(0.4)
end
ylabel("Total mRNA(M)")
xlabel("v_s constant")
title("Sustained Oscillations between Total Protein and per mRNA (24 hr cycle)")

legend('IC 1 Solution', 'IC 2 Solution', 'IC 3 Solution', 'Poincare Diagram IC 1', 'Poincare Diagram IC 2', 'Poincare Diagram IC 3', 'Location', 'SouthOutside')

hold off





figure(3)
hold on
for k=1:(size(init_1))[2]    
    
    color = [0.9 - 0.4*(k-1), 0.4*(k-1), 0.9];
    p3 = plot(tss{k}, pt2{k},'LineWidth', (10 - 2.5*k));
    p3.Color(4) = 0.25;

end
for k=1:(size(init_1))[2]    
  color = [0.9 - 0.4*(k-1), 0.4*(k-1), 0.9];  
p3 = scatter(t12s{k}, totPhours12{k}, 80 - 10*k, 'filled');
    alpha(0.25)
end

legend('IC 1 Solution', 'IC 2 Solution', 'IC 3 Solution', 'location','southoutside')
ylabel("Total mRNA(M)")
xlabel("Time(hours)")
title("Total mRNA vs Time")
hold off
%%




%%
function P_Conc=getC(t,P,y, hours)
% Protein Function
% Outputs: Concentrations of 5 proteins of interest
% Inputs: Time, concentration, PER maximum degradation rate

v_m=0.65;   % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d=y;      % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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

v_s = 0.3785*cos(2 * 3.1416 / hours * t)+0.3785;
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
% Forcing function for vs


end