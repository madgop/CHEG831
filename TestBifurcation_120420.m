clear all
%% Constants and Function Definition

%v_s=y;     % um/hr; accumulation rate of per mRNA in cytosol
            % This is set as a variable later in code
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
n = 4;      % [-]; degree of cooperativity

% f is an anonymous function for the concentrations as defined in 
% equations 1a to 1e of Goldbeter
% Form is f = [f1; f2; f3; f4; f5] = [dM/dt; dP0/dt; dP1/dt; dP2/dt; dPN/dt]

f = @(P,y) [y*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
    k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
    V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
    V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
    k_1*P(4)-k_2*P(5)];

% jac is an anonymous Jacobian for the concentrations of proteins defined in f
% Form is jac = [df1/dM df1/dP0 df1/dP1 df1/dP2 df1/dPN;
%                df2/dM df2/dP0 df2/dP1 df2/dP2 df2/dPN;
%                df3/dM df3/dP0 df3/dP1 df3/dP3 df3/dPN;
%                df4/dM df4/dP0 df4/dP1 df4/dP3 df4/dPN;
%                df5/dM df5/dP0 df5/dP1 df5/dP3 df5/dPN]

jac = @(PConc,y) [-v_m.*K_m./((K_m+PConc(1)).^2) 0 0 0 -(y.*K_I^n)*n*(PConc(5)).^(n-1)./((K_I.^n+PConc(5).^n).^2);
    k_s -(V_1*K_14)./((K_14+PConc(2)).^2) V_2*K_14./((K_14+PConc(3)).^2) 0 0;
    0 (V_1.*K_14)/((K_14+PConc(2)).^2) (-V_2.*K_14)/((K_14+PConc(3)).^2) - (V_3.*K_14)./((K_14+PConc(3)).^2) (V_4.*K_14)./((K_14+PConc(4)).^2) 0;
    0 0 V_3.*K_14./((K_14+PConc(3)).^2) (-V_4.*K_14./((K_14+PConc(4)).^2)-k_1 - v_d.*K_d./((K_d+PConc(4)).^2)) k_2;
    0 0 0 k_1 -k_2];

% dfdv_s is the df/dp vector that we solve for in parametric continuation
% Equation for parametric continuation is: Jacobian*dz/dp = -df/dp 
% Form is dfdv_s = [df1/dvs; df2/dvs; df3/dvs; df4/dvs; df5/dvs] 

dfdv_s = @(P_conc,v_s) [K_I^n/(P_conc(5)^n+K_I^n);0;0;0;0];

%% Bifurcation Algorithm

v_s = [0.01:0.0001:3];  % um/hr; accumulation rate of per mRNA in cytosol
k = 1;                  % counter
i = 1;                  % counter indexing value of v_s parameter

% Syntax is ode45(function, timespan, initialconditions)
[t,P_Conc] = ode45(@(t,P)getC(t,P,v_s(i)),[0,1000],[0.6;0.5;1.8;0.65;1.2]);
P_array(:,1) = P_Conc(end,:);       % save final protein concentrations time point
                                    % in ode45 solver in array P_array
for j=2:length(v_s)
    Jacob = Jacobian(P_array(:,j-1),v_s(j-1));  % compute Jacobian
    eigJacob(:,j) = eig(Jacob);                 % solves for eigenvalues of the Jacobian
    detJacob = abs(det(Jacob));                 % solves for determinant of the Jacobian
    if detJacob < 1e-4
        tp(k) = j;
        k = k+1;
        break
    end
    
    % If the absolute value of the real component of any of the eigenvalues of the Jacobian is less than 1e-4
    if abs(real(eigJacob(1,j)))<1e-4 || abs(real(eigJacob(2,j)))<1e-4 || abs(real(eigJacob(3,j)))<1e-4 ||...
            abs(real(eigJacob(4,j)))<1e-4 || abs(real(eigJacob(5,j)))<1e-4
        hopfBif(i) = j;                             % Store v_s array indices producing any Jacobian eigenvalues
                                                    % (real part magnitude) less than 1e-4 in hopfBif array
                                                    % Would (j-1) be more useful as an index to hopfBif?
        i = i + 1;
    end
    delF = dfdv_s(P_array(:,j-1),v_s(j-1));         % solving for delF in continuation equation
                                                    % Jacob*dzdp = -delF
    dzdp = -inv(Jacob)*delF;                        % solving Jacob*dzdp = -delF
                                                    % consider left division instead? MATLAB recommends this over inv
                                                    % for solving systems of linear equations
                                                    % dzdp = Jacob\(-delF) ? 
    
    % Update protein concentration matrix (new column) for new values of v_s 
    % using first-order parametric continuation                                               
    P_array(:,j) = P_array(:,j-1)+dzdp*(v_s(j)-v_s(j-1)); 
    Jacob1 = det(Jacobian(P_array(:,j),v_s(j)));    % compute new Jacobian at new P_array and v_s values
    % Use Newton's method to solve for protein concentrations (P_array) at new v_s
    [xsol, iter] = Newton(@(P,y) f(P,v_s(j)), @(PConc,y) jac(PConc,v_s(j)), P_array(:,j), 10000, 1e-5);
    P_array(:,j) = xsol;
    k = k + 1;
end



%% Function that returns protein concentrations
function P_Conc=getC(t,P,y)
% Protein Function
% Outputs: Concentrations of 5 proteins of interest
% Inputs: Time, concentration, accumulation rate of per mRNA in cytosol(v_s)

v_s=y   ;   % um/hr; accumulation rate of per mRNA in cytosol
v_m=0.65;   % um/hr; max degradation rate of per mRNA in cytosol
K_m=0.5;    % um; Michaelis constant for cytosolic per mRNA
k_s=0.38;   % hr^-1; first-order rate constant for PER synthesis
v_d= 0.95;  % um/hr; maximum degradation rate of biphosphorylated PER (P2)
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

% Equations 1a through 1e of Goldbeter follow; initialize concentration vector here
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

%% Function that returns Jacobian Matrix
function dP = Jacobian(PConc,y)
% Function to calculate Jacobian of system

% Inputs: vector of protein concentrations
% (necessarily at a specific time point)
% and y (v_s parameter)

% Outputs: dP, the matrix of PConc vectors at various time points
% PConc vector is defined as follows, from original code: 
% PConc(1): Cytosolic concentration, M
% PConc(2): Unphosphorylated PER, P0
% PConc(3): Monophosphorylated PER, P1
% PConc(4): Biphosphorylated PER, P2
% PConc(5): Nuclear PER, PN

PConc = PConc(:)';  % Convert protein concentration array to row vector

v_s=y;      % um/hr; accumulation rate of per mRNA in cytosol
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
V_1=3.2;    % um/hr; V_1 through V_4 are maximum rates for kinases and phosphatases
V_2=1.58;   % um/hr
V_3=5;      % um/hr
V_4=2.5;    % um/hr
n=4;        % [-]; degree of cooperativity


% Following section computes jacobian matrix
% Defining [f1; f2; f3; f4; f5] as [dM/dt; dP0/dt; dP1/dt; dP2/dt; dPN/dt]
% Jacobian is [df1/dM df1/dP0 df1/dP1 df1/dP2 df1/dPN;
%              df2/dM df2/dP0 df2/dP1 df2/dP2 df2/dPN;
%              df3/dM df3/dP0 df3/dP1 df3/dP3 df3/dPN;
%              df4/dM df4/dP0 df4/dP1 df4/dP3 df4/dPN;
%              df5/dM df5/dP0 df5/dP1 df5/dP3 df5/dPN]

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
            dP(i,4) = -V_4*K_14/((K_14+PConc(4))^2) - k_1 - v_d*K_d/((K_d+PConc(4))^2);
            dP(i,5) = k_2;
        else
            dP(i,4) = k_1;
            dP(i,5) = -k_2;
        end
    end
end

%% Function that returns value of newton's method
function [x_sol iterations] = Newton(fun,jacobian,x_init,N_max,tolerance)
% Newton implements Newton's iterative method for the solution of fun(x)=0

% INPUT:
% fun: user supplied function of x, the function that the root(s) of which are sought
% fun has form [dM/dt; dP0/dt; dP1/dt; dP2/dt; dPN/dt]

% jacobian: user supplied function of x, the derivative of fun with respect to x
% x: is the independent variable
% x_init: is the initial guess for the independent variable
% x_init has form [M_0; P0_0; P1_0; P2_0; PN_0];
% N_max: is the maximum number of iterations (10 recommended)
% tolerance: is the absolute tolerance for the solution (10e-6 used previously)

% OUTPUT:
% x_sol: is the solution (if converged) last guess (if not converged) 
% x_sol has form [M; P0; P1; P2; PN];
% iterations: number of iterations (k) required for convergence within given tolerance;
% is equal to -(N_max+1) if no solution is found within N_max iterations;
% it takes a negative value -k if the procedure diverged after k iterations

iterations = 0;
x_sol = x_init;                             % initial solution guess
diff = jacobian(x_sol)\fun(x_sol);          % solves Jacobian(x_sol)*diff = fun(x_sol)
diff_mag_old = norm(diff);                  % computes norm of dx/dp array in J*dx/dp = -df/dp
x_sol = x_sol - diff;                       % update solution

    for k=1:N_max
        diff = jacobian(x_sol)\fun(x_sol);  % solves Jacobian(x_sol)*diff = fun(x_sol)
        diff_mag = norm(diff);
            if(diff_mag < tolerance)
                iterations = k;
                x_sol = x_sol - diff;
                break                       % In nested loops, break only exits 
                                            % from the loop in which it occurs
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