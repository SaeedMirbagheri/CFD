clc;
clear; 
clear all;
close all;
disp('Programmer: Seyed Saeed Mirbagheri ( Student number: 400126116 ). ')
Range=input('Choose one of two intervals:   1=[5 9 17 33 65]     2=[65 129 257 513 1025]   3=[513 1025 2049 4097 8193]: ');
%%%Domain segmentation
if Range==1
n=[5 9 17 33 65];
elseif Range==3
n=[257 513 1025 2049 4097]; 
else   
n=[65 129 257 513 1025];
end
Length_n=length(n);
L=1;
delta_e=L./(n-1);
%%%%%Define problem constants
Bi=0.1;
Ar=20;
Meu=(Ar*Bi)^2;
I=1;

AA=zeros(1,1);
N=zeros(1,1);
BB=zeros(1,1);
%%%%%%% Formation of coefficient matrix (AA) and information vector (BB)

for I=1:Length_n
e=Vector_x(delta_e(I),L);

N(I)=length(e);
AA(1,1,I)=1;
BB(1,I)=1;

for i=2:N(I)-1
    for j=1:N(I)
        if i==j
            AA(i,j,I)=-(2+(Meu*delta_e(I))^2);
            
        elseif i==j+1
            AA(i,j,I)=1;
        
        elseif i==j-1
            AA(i,j,I)=1;
        end
    end
    
    
      BB(i,I)=0;
end
AA(N(I),N(I)-1,I)=2;

AA(N(I),N(I),I)=-(2+(Meu*delta_e(I))^2+2*delta_e(I)*Bi);


BB(N(I),I)=0;

end
A1=sparse(AA(1:N(1),1:N(1),1));
A2=sparse(AA(1:N(2),1:N(2),2));
A3=sparse(AA(1:N(3),1:N(3),3));
A4=sparse(AA(1:N(4),1:N(4),4));
A5=sparse(AA(1:N(5),1:N(5),5));

B1=BB(1:N(1),1);
B2=BB(1:N(2),2);
B3=BB(1:N(3),3);
B4=BB(1:N(4),4);
B5=BB(1:N(5),5);

Teta1=A1\B1;
Teta2=A2\B2;
Teta3=A3\B3;
Teta4=A4\B4;
Teta5=A5\B5;

%%%%%Call a function to divide the domain into different n
e1=Vector_x(delta_e(1),L);
e2=Vector_x(delta_e(2),L);
e3=Vector_x(delta_e(3),L);
e4=Vector_x(delta_e(4),L);
e5=Vector_x(delta_e(5),L);

%%%%Draw in ? the function e
figure(1)
scatter(e1,Teta1,'bo','linewidth',2)
title('Draw the \theta variable according to position e for different intervals');
xlabel('e');
ylabel('\theta');
hold on
grid on
scatter(e2,Teta2,'rs','linewidth',2)
scatter(e3,Teta3,'k*','linewidth',2)
scatter(e4,Teta4,'gd','linewidth',2)
scatter(e5,Teta5,'mx','linewidth',2)
legend('n=65','n=101','n=201','n=501','n=1001')

%%%%%Calling the point return function (0.25, 0.5, 0.75 and 1)
n1=Number_Teta(N(1),e1);
n2=Number_Teta(N(2),e2);
n3=Number_Teta(N(3),e3);
n4=Number_Teta(N(4),e4);
n5=Number_Teta(N(5),e5);

%%%%Calculate the relative error value
error=zeros(1,1);
for i=1
    for j=1:4
        error(i,j)=Teta2(n2(j))-Teta1(n1(j));
    end
end
for i=2
    for j=1:4
        error(i,j)=Teta3(n3(j))-Teta2(n2(j));
    end
end
for i=3
    for j=1:4
        error(i,j)=Teta4(n4(j))-Teta3(n3(j));
    end
end
for i=4
    for j=1:4
        error(i,j)=Teta5(n5(j))-Teta4(n4(j));
    end
end

%%%%%Draw the relative error value in terms of H
figure(2)
loglog(delta_e(1:end-1),abs(error(:,1)),'bo','linewidth',2);
title('Solution convergence');
ylabel('successive error');
xlabel('h');
hold on
grid on
loglog(delta_e(1:end-1),abs(error(:,2)),'rs','linewidth',2);
loglog(delta_e(1:end-1),abs(error(:,3)),'k*','linewidth',2);
loglog(delta_e(1:end-1),abs(error(:,4)),'gd','linewidth',2);
legend('e=0.25','n=0.5','n=0.75','n=1')

%%%%%%%%Error slope calculations
R_e=zeros(1,1);
R_delta_e=zeros(1,1);
R_e_delta_e=zeros(1,1);
for i=1:3
    for j=1:4
        R_e(i,j)=log(error(i+1,j)/error(i,j));

    end
        R_delta_e(i)=log(delta_e(i+1)/delta_e(i));

        
    for j=1:4
        R_e_delta_e(i,j)=R_e(i,j)/R_delta_e(i);
    end
end

%%%%Draw the error slope in terms of h
figure(3)
plot(delta_e(1:end-2),abs( R_e_delta_e(:,1)),'b','linewidth',2);
title('Error slope');
ylabel('Slope');
xlabel('h');
hold on
grid on
plot(delta_e(1:end-2),abs( R_e_delta_e(:,2)),'r','linewidth',2);
plot(delta_e(1:end-2),abs( R_e_delta_e(:,3)),'k','linewidth',2);
plot(delta_e(1:end-2),abs( R_e_delta_e(:,4)),'g','linewidth',2);
legend('e=0.25','n=0.5','n=0.75','n=1')

%%%%Calculate the actual error
figure(4)
Teta_Real=(cosh(Meu*(1-e))+(Bi/Meu)*sinh(Meu*(1-e)))/(cosh(Meu)+(Bi/Meu)*sinh(Meu)); 
Error_Real=Teta_Real-Teta5';
plot(e,Error_Real,'r','linewidth',2)
title('Real error');
ylabel('Error_Real');
xlabel('h');
grid on
