% CHEG 831 Project

clc
clear all

% Initial concentrations of the proteins evaluated in equations 1a through 1e
init_1=[0.1 0.25 0.25 0.25 0.25;...
        1.9 0.8 0.8 0.8 0.8;...
        0.7 0.5 0.5 0.5 0.5];

%%
% Figure 1 (Figure 2 of paper)

% Format for ode solver is oe45(function, time_span, initial_cond)
[t,P_Conc] = ode45(@(t,P)getC(t,P,0.65),[0,1000],[0.6;0.5;1.8;0.65;1.2]);
[rows,columns]=size(P_Conc);
for j=1:rows
    P_t(j,1)=(sum(P_Conc(j,2:5)));  % Total PER protein, eq. 2 in paper
end

P_Conc(:,6)=P_t;
for i=1:6
    figure(1)
    plot(t,P_Conc(:,i)) % Plotting all protein concentrations vs time
    hold on
end

ylabel("Per forms of M")
xlabel("Time(h^-1)")
legend("M","Po","P1","P2","PN","Pt")
title("Sustained oscillations of the model equations")

% %%
% % Figure 2 (Figure 3 of paper)
% 
% for k=1:size(init_1)[2] 
% % size(init_1)[2] is number of rows in init_1 (dimension 2)
% P_Conc_2=[] % New vector of protein concentrations
% P_M=[]      % PER mRNA
% P_t2=[]     % total PER protein 
% [t_2,P_Conc_2] = ode45(@(t,P)getC(t,P,0.95),[0,1000],init_1(k,:));
% [rows_2,columns]=size(P_Conc_2);
% size(P_Conc_2)
%     for j=1:rows_2
%         P_t2(j,k)=(sum(P_Conc_2(j,2:5)));   % Total PER protein
%         P_M(j,k)=P_Conc_2(j,1);             % PER mRNA
%     end
% figure(2)
% hold on
% plot(P_t2(:,k),P_M(:,k))
% ylabel("per mRNA(M)")
% xlabel("Total PER protein Pt")
% title("Sustained Oscillations between Total Protein and per mRNA")
% end

%%
% Figure3 (Figure 4 of paper)

count = 1;
for i = 0.5:0.05:2.75
    [t3,P_Conc_3] = ode45(@(t,P)getC(t,P,i),[0,1000],[0.1,0.25,0.25,0.25,0.25]);
    [rows_3,columns_3]=size(P_Conc_3);
    P_t3=zeros(rows_3,1);
    for j=1:rows_3
        P_t3(j,1)=sum((P_Conc_3(j,2:5)));        % Total amount of PER protein
    end
    P_Conc_3(:,6)=P_t3;     
    Peak_f = P_Conc_3(:,6);
    % Find peaks (y values) and locations (x)
    [peaks,locs] = findpeaks(Peak_f);           
    if count>1
        period(:,count) = min(diff(t3(locs)));  % Period of PER oscillations
    end
    count = count+1;
    
end

% V_D=linspace(0.5,2.75,46);  % PER maximum degradation rate
% figure(3)
% plot(V_D(2:end),period(1,2:end))
% ylabel("Period(h^-1)")
% xlabel("V_D")
% title("Dependence of PER oscillations on the maximum rate of PER degradation")



%%
a1stable = [];
z1stable = [];
a1unstab = [];
z1unstab = [];

a2stable = [];
z2stable = [];
a2unstab = [];
z2unstab = [];

a3stable = [];
z3stable = [];
a3unstab = [];
z3unstab = [];

v_s = [0:0.01:5];

% for i = 1:length(v_s)
%     [t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(i),[0,72],[0.6;0.5;1.8;0.65;1.2]);
%     J1 = Jacobian(PConc,v_s(i));
%     eigJ1 = eig(J1);
%     if (lambdaJ11 < 0) && (lambdaJ12 < 0)
%         fprintf('SS solution at z=(0,0) is stable for alpha = %d\n', alpha(i))
%         a1stable = [a1stable, alpha(i)];
%         z1stable = [z1stable, 0];
%     elseif (lambdaJ11 > 0) && (lambdaJ12 > 0)
%         fprintf('SS solution at z=(0,0) is completely unstable for alpha = %d\n', alpha(i))
%         a1unstab = [a1unstab, alpha(i)];
%         z1unstab = [z1unstab, 0];
%     else
%         fprintf('SS solution at z=(0,0) has a saddle point for alpha = %d\n', alpha(i))
%         a1unstab = [a1unstab, alpha(i)];
%         z1unstab = [z1unstab, 0];            
%     end    
% lambdaJ21 = 0.5*(tr2+sqrt(tr2^2-4*detJ2));
% lambdaJ22 = 0.5*(tr2-sqrt(tr2^2-4*detJ2));
%     if (lambdaJ21 < 0) && (lambdaJ22 < 0)
%         fprintf('SS solution at z=(1,0) is stable for alpha = %d\n', alpha(i))
%         a2stable = [a2stable, alpha(i)];
%         z2stable = [z2stable, 1];
%     elseif (lambdaJ21 > 0) && (lambdaJ22 > 0)
%         fprintf('SS solution at z=(1,0) is completely unstable for alpha = %d\n',alpha(i))
%         a2unstab = [a2unstab, alpha(i)];
%         z2unstab = [z2unstab, 1];
%     else
%         fprintf('SS solution at z=(1,0) has a saddle point for alpha = %d\n',alpha(i))
%         a2unstab = [a2unstab, alpha(i)];
%         z2unstab = [z2unstab, 1];
%     end
% lambdaJ31 = 0.5*(tr3+sqrt(tr3^2-4*detJ3));
% lambdaJ32 = 0.5*(tr3-sqrt(tr3^2-4*detJ3));
%     if (lambdaJ31 < 0) && (lambdaJ31 < 0)
%         fprintf('SS solution at z=(1/alpha , 1/beta-1/(alpha*beta)) is stable for alpha = %d\n',alpha(i))
%         a3stable = [a3stable, alpha(i)];
%         z3stable = [z3stable, 1/alpha(i)];
%     elseif (lambdaJ31 > 0) && (lambdaJ32 > 0)
%         fprintf('SS solution at z=(1/alpha , 1/beta-1/(alpha*beta)) is completely unstable for alpha = %d\n',alpha(i))
%         a3unstab = [a3unstab, alpha(i)];
%         z3unstab = [z3unstab, 1/alpha(i)];
%     else
%         fprintf('SS solution at z=(1/alpha , 1/beta-1/(alpha*beta)) has a saddle point for alpha = %d\n',alpha(i))
%         a3unstab = [a3unstab, alpha(i)];
%         z3unstab = [z3unstab, 1/alpha(i)];
%     end
% 
% end
% 

%%
function P_Conc=getC(t,P,y)
% Protein Function
% Outputs: Concentrations of 5 proteins of interest
% Inputs: Time, concentration, PER maximum degradation rate

v_s=y   ;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;   % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d=0.95;   % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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
end

