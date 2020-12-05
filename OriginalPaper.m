clc
clear all
init_1=[0.1 0.25 0.25 0.25 0.25;1.9 0.8 0.8 0.8 0.8;0.7 0.5 0.5 0.5 0.5];   % Eric Git test
%%Figure1
[t,P_Conc] = ode45(@(t,P)getC(t,P,0.95),[0,72],[0.6;0.5;1.8;0.65;1.2]);

[rows,columns]=size(P_Conc);
for j=1:rows
P_t(j,1)=(sum(P_Conc(j,2:5)));
end
P_Conc(:,6)=P_t;
for i=1:6
    figure(1)
    plot(t,P_Conc(:,i))
    hold on
end
ylabel("Protein Concentration")
xlabel("Time")
legend("M","Po","P1","P2","PN","Pt")
%%Figure2
for k=1:size(init_1)[2]
P_Conc_2=[]
P_M=[]
P_t2=[]
[t_2,P_Conc_2] = ode45(@(t,P)getC(t,P,0.95),[0,1000],init_1(k,:));
[rows_2,columns]=size(P_Conc_2);
size(P_Conc_2)
    for j=1:rows_2
        P_t2(j,k)=(sum(P_Conc_2(j,2:5)));
        P_M(j,k)=P_Conc_2(j,1);
    end
figure(2)
hold on
plot(P_M(:,k),P_t2(:,k))
ylabel("Pt")
xlabel("M")
end

%%Figure3
count = 1;
for i = 0.5:0.05:2.75
    [t3,P_Conc_3] = ode45(@(t,P)getC(t,P,i),[0,1000],[0.1,0.25,0.25,0.25,0.25]);
    [rows_3,columns_3]=size(P_Conc_3);
    P_t3=zeros(rows_3,1);
    for j=1:rows_3
    P_t3(j,1)=sum((P_Conc_3(j,2:5)));
    end
    P_Conc_3(:,6)=P_t3;
    Peak_f = P_Conc_3(:,6);
    [peaks,locs] = findpeaks(Peak_f);
    if count>1
    period(:,count) = min(diff(t3(locs)));
    end
    count = count+1;
    
end
V_D=linspace(0.5,2.75,46);
figure(3)
plot(V_D,period(1,:))
ylabel("Period(h^-1)")
xlabel("V_D")
%%Protein Function
function P_Conc=getC(t,P,y)
v_s=0.76;
v_m=0.65;
K_m=0.5;
k_s=0.38;
v_d=y;
k_1=1.9;
k_2=1.3;
K_I=1;
K_d=0.2;
K_14=2;
V_1=3.2;
V_2=1.58;
V_3=5;
V_4=2.5;
n=4;
P_Conc=zeros(5,1);
P_Conc(1)=v_s*(K_I.^n/(K_I.^n+P(5).^n))-v_m*(P(1)/(K_m+P(1)));
P_Conc(2)=k_s*P(1)-V_1*P(2)/(K_14+P(2))+V_2*(P(3)/(K_14+P(3)));
P_Conc(3)=V_1*(P(2)/(K_14+P(2)))-V_2*(P(3)/(K_14+P(3)))-V_3*(P(3)/(K_14+P(3)))+V_4*(P(4)/(K_14+P(4)));
P_Conc(4)=V_3*(P(3)/(K_14+P(3)))-V_4*(P(4)/(K_14+P(4)))-k_1*P(4)+k_2*P(5)-v_d*(P(4)/(K_d+P(4)));
P_Conc(5)=k_1*P(4)-k_2*P(5);
end

