function [  Teta,N,NN  ] = BTCS(Bi,Ar, L , Fo_end , Teta_zero , Teta_infinite , Length_n,Length_delta_Fo , delta_e , delta_Fo, Range )
N=zeros(1,1);
NN=zeros(1,1);
AA=zeros(1,1);
Teta=zeros(1,1);
BB=zeros(1,1);

if Range==1    
for I=1:Length_n
    
e=Vector_x(delta_e(I),L);
Fo=Vector_x(delta_Fo(I),Fo_end);
N(I)=length(e);
NN(I)=length(Fo);
AA(1,1,I)=1;

Teta(1,1,I)=Teta_zero;
Teta(2:N(I),1,I)=Teta_infinite;
    
for J=2:NN(I)
    Teta(1,J,I)=Teta_zero;
    AA(1,1,I)=1+(2*delta_Fo(I))/(delta_e(I)^2)+Bi*Ar*delta_Fo(I);
    BB(1,I)=Teta(2,J-1,I)+((delta_Fo(I))/(delta_e(I)^2))*Teta(1,J,I);
    AA(1,2,I)=-(delta_Fo(I))/(delta_e(I)^2);

for i=2:N(I)-2
    for j=1:N(I)-1
        
        if i==j
            AA(i,j,I)=1+(2*delta_Fo(I))/(delta_e(I)^2)+Bi*Ar*delta_Fo(I);
            
        elseif j==i+1
            AA(i,j,I)=-(delta_Fo(I))/(delta_e(I)^2);
        
        elseif j==i-1
            AA(i,j,I)=-(delta_Fo(I))/(delta_e(I)^2);
        end

    end
    
      BB(i,I)=Teta(i+1,J-1,I);
end
AA(N(I)-1,N(I)-2,I)=-(2*delta_Fo(I))/(delta_e(I)^2);
AA(N(I)-1,N(I)-1,I)=1+(2*delta_Fo(I))/(delta_e(I)^2)+Bi*Ar*delta_Fo(I)+(2*delta_Fo(I)*Bi)/(delta_e(I));
BB(N(I)-1,I)=Teta(N(I),J-1,I);
teta=AA(:,:,I)\BB(:,I);
Teta(2:N(I),J,I)=teta;
end
end

else
for I=1:Length_delta_Fo         
e=Vector_x(delta_e(1),L);
Fo=Vector_x(delta_Fo(I),Fo_end);
N(I)=length(e);
NN(I)=length(Fo);
AA(1,1,I)=1;

Teta(1,1,I)=Teta_zero;
Teta(2:N(I),1,I)=Teta_infinite;
    
for J=2:NN(I)
    Teta(1,J,I)=Teta_zero;
    AA(1,1,I)=1+(2*delta_Fo(I))/(delta_e(1)^2)+Bi*Ar*delta_Fo(I);
    BB(1,I)=Teta(2,J-1,I)+((delta_Fo(I))/(delta_e(1)^2))*Teta(1,J,I);
    AA(1,2,I)=-(delta_Fo(I))/(delta_e(1)^2);

for i=2:N(I)-2
    for j=1:N(I)-1
        
        if i==j
            AA(i,j,I)=1+(2*delta_Fo(I))/(delta_e(1)^2)+Bi*Ar*delta_Fo(I);
            
        elseif j==i+1
            AA(i,j,I)=-(delta_Fo(I))/(delta_e(1)^2);
        
        elseif j==i-1
            AA(i,j,I)=-(delta_Fo(I))/(delta_e(1)^2);
        end

    end
    
      BB(i,I)=Teta(i+1,J-1,I);
end
AA(N(I)-1,N(I)-2,I)=-(2*delta_Fo(I))/(delta_e(1)^2);
AA(N(I)-1,N(I)-1,I)=1+(2*delta_Fo(I))/(delta_e(1)^2)+Bi*Ar*delta_Fo(I)+(2*delta_Fo(I)*Bi)/(delta_e(1));
BB(N(I)-1,I)=Teta(N(I),J-1,I);
teta=AA(:,:,I)\BB(:,I);
Teta(2:N(I),J,I)=teta;
end
end
end

