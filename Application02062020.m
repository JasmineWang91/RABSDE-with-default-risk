clear all
clear
clc

N = 400;
T = 1;
p = 10000; % panalization parameter
stepsize = T/N;
%anticipate = rand(1,1); % anticipated part 
anticipate = 0.3;
tau = 0.2;
T_an = T+anticipate; % anticipated part (total time)
n_an = ceil(anticipate/stepsize); % anticipated part (extra steps)
N_an = N+n_an; % anticipated part (total steps)
t = (0:stepsize:T)'; % discrete time sequence t=[0,delta,2 delta,..., T]
%t_an = [t,((T+stepsize):stepsize:T_an)]; % discrete time sequence st=[0,delta,2 delta,..., T+anticipate]
t_an = (0:stepsize:(stepsize*N_an))';

%------------ explicit penalization scheme  ---------


%------------ explicit reflected scheme  ---------



%-------------------- 5.1 Discrete Time Framework ---------------------

%--------- 5.1.1 Random Walk Approximation of Brownian Motion -------

%------ (1) create Bernoulli sequence {1,-1}, P(1)=P(-1)=0.5 -------
e=rand(1,N_an);

for i=1:N_an
    if e(i)<0.5
        e(i)=-1;
    else
        e(i)=1;
    end
end

%------ (2) create Brownian Motion ------
E = cumsum(e);

B(1) = 0;  % first element = 0

for i = 2:(N_an+1)
    B(i) = sqrt(stepsize)*E(i-1);  % Brownian Motion
    b(i-1) = B(i)-B(i-1);  % B(i)-B(i-1)
end

%-------------- 5.1.2 Approximation of Defaultable Model ------------
%tau = rand(1,1);  % default time

h = zeros(N+1,1); 
h(1) = 0;
for i = 2:(N+1)
    if t(i) < tau
        h(i) = 0;
    else
        h(i) = 1;
    end    
end

%------------------ default martingale ------------------
%sh = cumsum(1-h);

M(1) = 0;  % first element = 0   
gamma(1) = 0;

for i = 2:(N+1)
    if h(i-1) == 0 && h(i) == 0
    gamma(i) = 1/(T-t(i));
    m(i-1) = h(i)-h(i-1)-stepsize*gamma(i)*(1-h(i));  % M(i)-M(i-1)
    
    elseif h(i-1) == 0 && h(i) == 1
        gamma(i) = 1/(T-t(i));
        m(i-1) = h(i)-h(i-1)-stepsize*gamma(i)*(1-h(i));  % M(i)-M(i-1)
        
    elseif h(i-1) == 1
        gamma(i) = 0;
        m(i-1) = 0;
    end
end

M = [0, cumsum(m)];  % default martingale

%--------- 5.1.4 Approximation of Terminal Value and Generator -----------
C(1) = 1.5;
S(1) = 1.5;
S_1(1) = 1.5;
r = 1;
for i = 1: N
    C(i+1) = C(i) + 0.3*C(i)*stepsize;
    S(i+1) = S(i) + S(i)*(1.1*stepsize + 0.5*b(i) + 0.2*m(i));
    S_1(i+1) = S_1(i) + S_1(i)*(1.1*stepsize + 0.5*b(i));
end

for i = 1:(N+1)
    L(i) = exp(-r*t(i))*max(S(i)-1,0); % lower obstacle
    V(i) = exp(-r*t(i))*2*max(S(i)-1,0); % upper obstacle
    
    L_1(i) = exp(-r*t(i))*max(S_1(i)-1,0); % lower obstacle
    V_1(i) = exp(-r*t(i))*2*max(S_1(i)-1,0); % upper obstacle
  
    %L(i) = -abs(B(i))+1/2*M(i)+T; % lower obstacle
    %V(i) = abs(B(i)+1)+1/2*M(i)+2*T-1/2; % upper obstacle
end   

xi = exp(-r*t(N+1))*1.2*max(S(N+1)-1,0);
xi_1 = exp(-r*t(N+1))*1.2*max(S_1(N+1)-1,0);
%---------- numerical scheme: explicit penalization scheme ---------------
y_pena(N+1)= xi; % Terminal Value
y_ref(N+1)= xi;
y_ref_1(N+1)= xi_1;


%-------------------------------------------------------------------------
for i = N:-1:1         
%---------- anticipated part definition   ------------------------------         
    
%---------- explicit part ---------------
    if  h(i) == 1 && h(i+1) == 1
        Ey_pena(i) = y_pena(i+1);
        Ey_ref(i) = y_ref(i+1);
        u_pena(i) = 0; 
        u_ref(i) = 0;        
        
    elseif  h(i) == 0 && h(i+1) == 1
        Ey_pena(i) = (stepsize/((T-t(i))))*y_pena(i+1);
        Ey_ref(i) = (stepsize/((T-t(i))))*y_ref(i+1);
        u_pena(i) = (stepsize*y_pena(i+1))/((stepsize+(stepsize*gamma(i+1))^2*(T-t(i+1))));
        u_ref(i) = (stepsize*y_ref(i+1))/((stepsize+(stepsize*gamma(i+1))^2*(T-t(i+1))));
        
    elseif h(i) == 0 && h(i+1) == 0
        Ey_pena(i) = ((T-t(i+1))/((T-t(i))))*y_pena(i+1);  
        Ey_ref(i) = ((T-t(i+1))/((T-t(i))))*y_ref(i+1);      
        u_pena(i) = -((T-t(i+1))*stepsize*gamma(i+1)*y_pena(i+1))/((stepsize+(stepsize*gamma(i+1))^2*(T-t(i+1))));    
        u_ref(i) = -((T-t(i+1))*stepsize*gamma(i+1)*y_ref(i+1))/((stepsize+(stepsize*gamma(i+1))^2*(T-t(i+1))));          
    end        
   
    
%---------- panaliztion solution ---------------
  z_pena(i) = e(i+1)*(y_pena(i+1)/(2*sqrt(stepsize)));
  
  k1_pena(i) = -min((Ey_pena(i)-L(i)),0)*((p*stepsize)/(1+p*stepsize)); 
  k2_pena(i) = max((Ey_pena(i)-V(i)),0)*((p*stepsize)/(1+p*stepsize)); 
  
  y_pena(i) = Ey_pena(i)+k1_pena(i)-k2_pena(i);  
  
%---------- reflected solution ---------------  
  z_ref(i) = e(i+1)*(y_ref(i+1)/(2*sqrt(stepsize)));
  
  k1_ref(i) = -min((Ey_ref(i)-L(i)),0); 
  k2_ref(i) = max((Ey_ref(i)-V(i)),0); 
  
  y_ref(i) = Ey_ref(i)+k1_ref(i)-k2_ref(i); 
  
%-------------------- without default -------------------------------
  z_ref_1(i) = e(i+1)*(y_ref_1(i+1)/(2*sqrt(stepsize)));
  Ey_ref_1(i) = y_ref_1(i+1);
  
  k1_ref_1(i) = -min((Ey_ref_1(i)-L_1(i)),0); %------ without default
  k2_ref_1(i) = max((Ey_ref_1(i)-V_1(i)),0); %------ without default
  
  y_ref_1(i) = Ey_ref_1(i)+k1_ref_1(i)-k2_ref_1(i);  
   
end


K1_pena = cumsum(k1_pena); 
K2_pena = cumsum(k2_pena);
K1_ref = cumsum(k1_ref); 
K2_ref = cumsum(k2_ref);


%--------------- outcome ---------------- 
xi(1)
L(1)
V(1)
L(N+1)
V(N+1)
y_pena(1)
y_ref(1)
y_ref_1(1)


%--------------- creates graphs ---------------- 
figure(1)  
plot(t,B(1:N+1)) % Brownian motion
title('Trajectories of the Brownian motion');
xlabel('t')
ylabel('B')

figure(2)  
plot(t,M) % Default martinagle
legend('Default martinagle');
title('Trajectories of the default martinagle');
xlabel('t')
ylabel('M')

%subplot(2,2,3); % terminal value xi
%plot(t_an(N:N_an),xi)
%title('Trajectories of the anticipated process');
%xlabel('t')
%ylabel('xi')

figure(3) 
plot(t,L)
hold on 
plot(t,V)
hold on 
plot(t,y_pena(1:N+1))
legend('Lower obstacle','Upper obstacle','Penalization solution');
title('Trajectories of the explicit penalization solution');
xlabel('t')
ylabel('Y')

figure(4) 
plot(t,L)
hold on 
plot(t,V)
hold on 
plot(t,y_ref(1:N+1))
legend('Lower obstacle','Upper obstacle','Reflected solution');
title('Trajectories of the explicit reflected solution');
xlabel('t')
ylabel('')

figure(5) 
plot(t,L_1)
hold on 
plot(t,V_1)
hold on 
plot(t,y_ref_1(1:N+1))
legend('Lower obstacle','Upper obstacle','Reflected solution');
title('Trajectories of the explicit reflected solution (without default)');
xlabel('t')
ylabel('')


figure(6) 
plot(t(1:N),K1_pena)
hold on
plot(t(1:N),K2_pena)
hold on
plot(t(1:N),y_pena(1:N))
legend('k+','k-','Penalization solution');
title('Trajectories of increasing processes of the explicit penalization scheme');
xlabel('t')
ylabel('K')

figure(7) 
plot(t(1:N),K1_ref)
hold on
plot(t(1:N),K2_ref)
hold on
plot(t(1:N),y_ref(1:N))
legend('k+','k-','reflected solution');
title('Trajectories of increasing processes of the explicit reflected scheme');
xlabel('t')
ylabel('K')

figure(8) 
plot(t,y_pena(1:N+1))
hold on
plot(t,y_ref(1:N+1))
legend('Penalization solution','Reflected solution');
title('Trajectories of solutions Y of the explicit penalization and explicit reflected scheme');
xlabel('t')
ylabel('Y')

figure(9) 
plot(t(1:N),z_pena)
hold on
plot(t(1:N),z_ref)
title('Trajectories of Z of the explicit reflected scheme');
xlabel('t')
ylabel('Z')

figure(10) 
plot(t(1:N),u_pena)
hold on
plot(t(1:N),u_ref)
title('Trajectories of U of the explicit reflected scheme');
xlabel('t')
ylabel('U')

figure(11) 
plot(t(1:N),S(1:N))
hold on
plot(t(1:N),C(1:N))
title('Trajectories of S and C');
xlabel('t')
ylabel('S')


