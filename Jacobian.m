function dP = Jacobian(PConc,y)
% Function to calculate Jacobian of system
% Inputs: vector of protein concentrations
% (necessarily at a specific time point)
% and y (v_d parameter)

% PConc vector is defined as follows, from original code: 
% PConc(1): Cytosolic concentration, M
% PConc(2): Unphosphorylated PER, P0
% PConc(3): Monophosphorylated PER, P1
% PConc(4): Biphosphorylated PER, P2
% PConc(5): Nuclear PER, PN

v_s=0.76;   % um/hr; accumulation rate of per mRNA in cytosol
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

dP = zeros(5,5);
    for i = 1:length(C)
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