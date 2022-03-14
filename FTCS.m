function [ Teta,N,NN ] = FTCS(Bi,Ar,L  ,Fo_end , Teta_zero , Teta_infinite , Length_n ,Length_delta_Fo , delta_e , delta_Fo , Range )
N=zeros(1,1);
NN=zeros(1,1);
Teta=zeros(1,1);
if Range==1
for I=1:Length_n
e=Vector_x(delta_e(I),L);
Fo=Vector_x(delta_Fo(I),Fo_end);
N(I)=length(e);
NN(I)=length(Fo);
Teta(2:N(I),1,I)=Teta_infinite;
Teta(1,1,I)=Teta_zero;

for j=1:NN(I)-1
       Teta(1,j+1,I)=Teta_zero;
       for i=2:N(I)-1
         
            Teta(i,j+1,I)=(delta_Fo(I)/((delta_e(I))^2))*Teta(i+1,j,I)+(delta_Fo(I)/(delta_e(I)^2))*Teta(i-1,j,I)-((2*delta_Fo(I)/(delta_e(I)^2))+delta_Fo(I)*Bi*Ar-1)*Teta(i,j,I);
       end
            Teta(N(I),j+1,I)=(2*delta_Fo(I)/(delta_e(I)^2))*Teta(N(I)-1,j,I)-((2*delta_Fo(I)/(delta_e(I)^2))+((2*Bi*delta_Fo(I))/(delta_e(I)))+delta_Fo(I)*Bi*Ar-1)*Teta(N(I),j,I);  
end
end
else
for I=1:Length_delta_Fo
    e=Vector_x(delta_e(1),L);
    Fo=Vector_x(delta_Fo(I),Fo_end);
    N(I)=length(e);
    NN(I)=length(Fo);
    Teta(2:N(I),1,I)=Teta_infinite;
    Teta(1,1,I)=Teta_zero;

for j=1:NN(I)-1
       Teta(1,j+1,I)=Teta_zero;
       for i=2:N(I)-1
         
            Teta(i,j+1,I)=(delta_Fo(I)/((delta_e(1))^2))*Teta(i+1,j,I)+(delta_Fo(I)/(delta_e(1)^2))*Teta(i-1,j,I)-((2*delta_Fo(I)/(delta_e(1)^2))+delta_Fo(I)*Bi*Ar-1)*Teta(i,j,I);
       end
            Teta(N(I),j+1,I)=(2*delta_Fo(I)/(delta_e(1)^2))*Teta(N(I)-1,j,I)-((2*delta_Fo(I)/(delta_e(1)^2))+((2*Bi*delta_Fo(I))/(delta_e(1)))+delta_Fo(I)*Bi*Ar-1)*Teta(N(I),j,I);  
end
end   
end
end

