%Simulations of Long-Period FBG for ION
%Adam Ignaciuk

%%%%%%%%%%%%%%%%%%%
%Dane z artykułu na później
a1 = 2.625 * 10e-6;
a2 = 62.4 * 10e-6;
a3 = a2+300*10e-9;
n_0 = 1.4681;
n_1 = 1.458;
n_2 = 1.45;
n_3 = 1.5;
n_4 = 1;
sigma_test = 4*10e-4;
Period_test = 450*10e-6;
L_test = 1.2*10e-2;
%%%%%%%%%%%%%%%%%%% 
%5CB liquid crystal - podają łączny czyli chyba n_eff, później znajdę dane
%z n_o i n_e
n_5CB_glass = 1.4990;
%%%%%%%%%%%%%%%%%%%
%Dane dla femtosekundowej siatki, o okresie nanometrowym
n_eff = 1.447;
n_modul = 0.0017; %modulacja RI
Period = 530*1e-9; %okres 
L = 0.001; %długość siatki 
Lambda_d=2*n_eff*Period; %Centralna Długość Fali = 1533,82nm

Lambda=(1533.80:0.0005:1533.84)*1e-9; %Fale
dz = 0.00005; %Zmiana odłegości do macierzy
z = (0 : dz : L); %Odległość na całej siatce 

delta=2*pi*n_eff*((1./Lambda)-(1./Lambda_d)); % w rzeczywistości, tą wielkość trzeba jeszcze przemnożyć przez detuning 
% detuning = tutaj mam mały problem 

sigma=(2*pi./Lambda)*n_modul;
kappa=(pi./Lambda)*n_modul;
delta_x= delta.*sigma; %wartość oznaczana jako delta z daszkiem 
gamma = sqrt(kappa.^2+delta_x.^2); %dla siatek działających tylko w transmisji 

T=[];
F_i=[];

for i=1:length(Lambda)
    F_i=[1 0;0 1];
    for j=1:length(z) %poniżej macierz dla siatek d. w transmisji
         
       F_i = F_i*[cos(dz*gamma(i))+1i*(delta_x(i)/gamma(i))*sin(gamma(i)*dz) 1i*(kappa(i)/gamma(i))*sin(gamma(i)*dz);1i*(kappa(i)/gamma(i))*sin(gamma(i)*dz) cos(gamma(i)*dz)-1i*(delta_x(i)/gamma(i))*sin(gamma(i)*dz)];
    end
    x = abs(1/F_i(1,1))^2; %to jest niby Transmisja 
    T=[T,x];
end

plot(Lambda,T);


