clc;
clear; 
clear all;
close all;
disp('Programmer: Seyed Saeed Mirbagheri ( Student number: 400126116 ). ')
method=input('Enter your desired solution method: 1=FTCS    2=CTCS   3=BTCS  4=Crank Niclson : ');
Range=input('Choose one of two intervals:    1=Error analysis for space      2=Error analysis for time  : ');
%%%Domain segmentation

L=1;
Fo_end=1;
if Range==1
nwo=0.45;   
n=[5  9 17 33 65];
delta_e=L./(n-1);
delta_Fo=nwo.*delta_e.^2;
Length_n=length(n);
Length_delta_Fo=length(delta_Fo);

else 
nwo=[0.45,0.225,0.1125,0.05625,0.028125];
n=33;
delta_e=L/(n-1);
delta_Fo=nwo.*(delta_e^2);
Length_n=length(delta_e);
Length_delta_Fo=length(delta_Fo);

end



%%%%%Define problem constants
Bi=0.1;
Ar=20;
Meu=(Ar*Bi)^2;
Teta_infinite=0;
Teta_zero=1;
I=1;

%%%%%%% Formation of coefficient matrix (AA) and information vector (BB)
if method==1
[Teta,N,NN]=FTCS(Bi,Ar,L  ,Fo_end , Teta_zero , Teta_infinite , Length_n ,Length_delta_Fo , delta_e , delta_Fo , Range ) ;
elseif method==2
[Teta,N,NN]=CTCS(Bi,Ar,L  ,Fo_end , Teta_zero , Teta_infinite , Length_n ,Length_delta_Fo , delta_e , delta_Fo , Range ) ;
elseif method==3
[Teta,N,NN]=BTCS(Bi,Ar, L , Fo_end , Teta_zero , Teta_infinite , Length_n ,Length_delta_Fo , delta_e , delta_Fo, Range ) ;
elseif method==4  
[Teta,N,NN]=Crank_Nicholson(Bi,Ar, L , Fo_end , Teta_zero , Teta_infinite , Length_n,Length_delta_Fo , delta_e , delta_Fo,Range ) ;
end

if Range==1
%%%%%Call a function to divide the domain into different n
e=zeros(N(Length_n),N(Length_n));
for i=1:Length_n
e(i,1:N(i))=Vector_x(delta_e(i),L);
end

[error,nn]=Error(N,NN,Length_n,Length_delta_Fo,e,Teta,Range );
%%%%Draw in \theta the function e
figure(1)
plot(e(1,1:N(1)),Teta(1:N(1),NN(1),1),'bo','linewidth',2)
if Range==1
title('Draw the \theta variable according to position e for different intervals (Fo=1 , \nu=0.45)');
else
title('Draw the \theta variable according to position e for different intervals (Fo=1 , n=5)');
end
xlabel('e');
ylabel('\theta');
grid on
hold on
plot(e(2,1:N(2)),Teta(1:N(2),NN(2),2),'rs','linewidth',2)
plot(e(3,1:N(3)),Teta(1:N(3),NN(3),3),'k*','linewidth',2)
plot(e(4,1:N(4)),Teta(1:N(4),NN(4),4),'gd','linewidth',2)
plot(e(5,1:N(5)),Teta(1:N(5),NN(5),5),'mx','linewidth',2)
legend('n=5','n=9','n=17','n=33','n=65')

%%%%%Calling the point return function (0.25, 0.5, 0.75 and 1)


%%%%%Draw the relative error value in terms of H
figure(2)
loglog(delta_e(1:end-1),abs(error(:,1)),'bo','linewidth',2);
title('Solution convergence (Fo=1 , \nu=0.45))');
ylabel('successive error');
xlabel('h');
hold on
grid on
loglog(delta_e(1:end-1),abs(error(:,2)),'rs','linewidth',2);
loglog(delta_e(1:end-1),abs(error(:,3)),'k*','linewidth',2);
loglog(delta_e(1:end-1),abs(error(:,4)),'gd','linewidth',2);
legend('e=0.25','n=0.5','n=0.75','n=1')

%%%%%%%%Error slope calculations
 R_e_delta_e=Error_Slope( error,delta_e );
%%%%Draw the error slope in terms of h
figure(3)
plot(delta_e(1:end-2),abs( R_e_delta_e(:,1)),'b','linewidth',2);
title('Error slope (Fo=1 , \nu=0.45)');
ylabel('Slope');
xlabel('h');
hold on
grid on
plot(delta_e(1:end-2),abs( R_e_delta_e(:,2)),'r','linewidth',2);
plot(delta_e(1:end-2),abs( R_e_delta_e(:,3)),'k','linewidth',2);
plot(delta_e(1:end-2),abs( R_e_delta_e(:,4)),'g','linewidth',2);
legend('e=0.25','n=0.5','n=0.75','n=1')



else   
 %%%%%Call a function to divide the domain into different n
e=zeros(N(Length_n),1);
for i=1:Length_n
e(i,1:N(i))=Vector_x(delta_e(i),L);
end



[error,nn]=Error(N,NN,Length_n,Length_delta_Fo,e,Teta,Range );
%%%%Draw in \theta the function e
figure(1)
plot(e(1,1:N(2)),Teta(1:N(1),NN(1),1),'bo','linewidth',2)
title('Draw the \theta variable according to position e for different intervals (Fo=1 , n=5)');
xlabel('e');
ylabel('\theta');
grid on
hold on
plot(e(1,1:N(2)),Teta(1:N(2),NN(2),2),'rs','linewidth',2)
plot(e(1,1:N(3)),Teta(1:N(3),NN(3),3),'k*','linewidth',2)
plot(e(1,1:N(4)),Teta(1:N(4),NN(4),4),'gd','linewidth',2)
plot(e(1,1:N(5)),Teta(1:N(5),NN(5),5),'mx','linewidth',2)
legend('\nu=0.45' ,'\nu=0.225' ,'\nu=0.1125' ,'\nu=0.05625' ,'\nu=0.028125')


%%%%%Calling the point return function (0.25, 0.5, 0.75 and 1)


%%%%%Draw the relative error value in terms of H
figure(2)
loglog(delta_Fo(1:end-1),abs(error(:,1)),'b','linewidth',2);
title('Solution convergence (Fo=1 , n=5)');
ylabel('successive error');
xlabel('delta_Fo');
hold on
grid on
loglog(delta_Fo(1:end-1),abs(error(:,2)),'r','linewidth',2);
loglog(delta_Fo(1:end-1),abs(error(:,3)),'k','linewidth',2);
loglog(delta_Fo(1:end-1),abs(error(:,4)),'g','linewidth',2);
legend('\nu=0.45' ,'\nu=0.225' ,'\nu=0.1125' ,'\nu=0.05625' )    
end

%%%%Calculate the actual error
