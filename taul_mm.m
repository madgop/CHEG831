% TAUL_MM.M
%
% Simple implementation of Tau Leaping Algorithm to simulate 
% the Michaelis-Menten system.
%
% Based on the scripts of D. Higham in SIAM Reviews, 2008.
% https://epubs.siam.org/doi/10.1137/060666457
% https://epubs.siam.org/doi/pdf/10.1137/060666457
%
% Parameters from Chapter 7 of 
%      Stochastic Modelling for Systems Biology,
%      by Darren J. Wilkinson, Chapman & Hall/CRC, 2006.
%

% Model Description:
%
% Chemical System is modeled with:
% substrate S1, enzyme S2, complex S3, and product S4
%
% System of equations is
%
%  S1 + S2 -(c1)-> S3
%  S3 -(c2)-> S1 + S2
%  S3 -(c3)-> S4 + S2

% stoichiometric matrix 
% each row is a species; each column is a reaction
V = [-1 1 0; -1 1 1; 1 -1 -1; 0 0 1];

%%%%%%%%%% Parameters and Initial Conditions %%%%%%%%%
nA = 6.023e23;             % Avagadro's number
vol = 1e-15;               % volume of system

Y = zeros(4,1);           % vector of molecules of each species
Y(1) = round(5e-7*nA*vol) % molecules of substrate, 'round' rounds to nearest integer
Y(2) = round(2e-7*nA*vol) % molecules of enzyme 

c(1) = 1e6/(nA*vol);      % vector of rate constants for each reaction
c(2) = 1e-4; 
c(3) = 0.1; 

tfinal = 50;       % simulation time
L = 10000;         % number of steps
tau = tfinal/L;    % step size

X =zeros(4,L);
X(:,1)=Y;

for k = 1:L                     % Loop over number of steps
     a(1) = c(1)*Y(1)*Y(2);     % Matrix of first-order forward reaction rate laws
     a(2) = c(2)*Y(3);
     a(3) = c(3)*Y(3);
 
     % Matrix of first-order propensity functions for each reaction
     % The if statement is needed because poissinv is not defined for lamda=0.
 
     if a(1)==0 d(1)=0; else d(1) = poissinv(rand,tau*a(1)); end
     if a(2)==0 d(2)=0; else d(2) = poissinv(rand,tau*a(2)); end
     if a(3)==0 d(3)=0; else d(3) = poissinv(rand,tau*a(3)); end

     % Equation 11 in linked PDF; 
     % Yi(t) is the state vector of the reaction, 
     % describing how species i (from 1 to 4 here) evolves with time
     Y = Y + d(1)*V(:,1) + d(2)*V(:,2) + d(3)*V(:,3);
     X(:,k+1)=Y;     % Populate state vector results for each time step here
end

time = (0:tau:tfinal);  %(initial time: time step : final time)

% plot substrate S1 and product S4 concentration versus time
plot(time,X(1,:),time,X(4,:)) 
title('Molecules of Substrate and Product from Tau Leaping Simulation of Michaelis-Menten Kinetics')
xlabel('Time, min')
ylabel('Molecules')
