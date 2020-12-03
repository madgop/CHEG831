clear all
%v_s=y;   % um/hr; accumulation rate of per mRNA in cytosol
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
v_s = [0.01:0.01:5];
count = 1;
count1 = 1;
K_I=1;      % um; threshold constant for repression
n=4;        % [-]; degree of cooperativity
dfdv_s = @(P_conc,v_s) [-n*K_I^n*P_conc(5)^(n-1)/(P_conc(5)^n+K_I^n)^2;0;0;0;0];

[t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(1)),[0,1000],[0.6;0.5;1.8;0.65;1.2]);
P_array(:,1) = P_Conc(end,:);

for j=2:3
    Jacob = Jacobian(P_array(:,j-1),v_s(j-1));
    eigJacob(:,j) = eig(Jacob);
    detJacob = abs(det(Jacob));
    if detJacob < 1e-5
        tp(count) = j;
        count = count+1;
        break
    end
    if abs(real(eigJacob(1))) <1e-5 | abs(real(eigJacob(2))) < 1e-5 | abs(real(eigJacob(3)))<1e-5 | abs(real(eigJacob(4)))<1e-5 | abs(real(eigJacob(5)))<1e-5
        hopfBif(count1) = j;
        count1 = count1 + 1;
    end
        delF = dfdv_s(P_array(:,j-1),v_s(j-1));
        dzdp = -inv(Jacob)*delF;
        P_array(:,j) = P_array(:,j-1)+dzdp*(v_s(j)-v_s(j-1));
        
        f = @(t,P,y) [y*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
            k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
            V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
            V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
            k_1*P(4)-k_2*P(5)];
        
        jac = @(PConc,y) [-v_m*K_m/((K_m.*PConc(1)).^2) 0 0 0 -(y*K_I.^n)/((K_I.^n+PConc(5).^n).^2);
            k_s -(V_1.*K_14)/((K_14+PConc(2)).^2) V_2.*K_14./((K_14+PConc(3)).^2) 0 0;
            0 (V_1.*K_14)/((K_14+PConc(2)).^2)-(V_2.*K_14)/((K_14+PConc(3)).^2) - (V_3.*K_14)./((K_14+PConc(3)).^2) (V_4.*K_14)./((K_14+PConc(4)).^2) 0;
            0 0 V_3.*K_14./((K_14+PConc(3)).^2) (-V_4.*K_14./((K_14+PConc(4)).^2)-k_1 - v_d.*K_d./((K_d+PConc(4)).^2)) k_2;
            0 0 0 k_1 -k_2];
                
        [xsol iter] = Newton(@(P,y) f(P,v_s(j)), @(PConc,y) jac(PConc,v_s(j)), P_array(:,j), 30, 1e-6);
        %Newton(fun,jacobian,x_init,N_max,tolerance)
        P_array(:,j) = xsol;
        count = count + 1;
end


function P_Conc=getC(t,P,y)
% Protein Function
% Outputs: Concentrations of 5 proteins of interest
% Inputs: Time, concentration, PER maximum degradation rate

v_s=y   ;   % um/hr; accumulation rate of per mRNA in cytosol
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
end

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

v_s=y;   % um/hr; accumulation rate of per mRNA in cytosol
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
            dP(i,1) = -v_m*K_m/((K_m*PConc(1))^2);
            dP(i,5) = -(v_s*K_I^n)/((K_I^n+PConc(5)^n)^2);
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

function [x_sol iterations] = Newton(fun,jacobian,x_init,N_max,tolerance)
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