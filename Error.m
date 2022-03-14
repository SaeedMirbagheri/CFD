function [ error,n ] = Error(N,NN,Length_n,Length_delta_Fo,e,Teta,Range )
n=zeros(Length_n,4);
for i=1:Length_n
n(i,1:4)=Number_Teta(N(i),e(i,1:N(i)));
end
error=zeros(1,1);
if Range==1
%%%%Calculate the relative error value
for i=1:Length_n-1
    for j=1:4
        error(i,j)=Teta(n(i+1,j),NN(i+1),i+1)-Teta(n(i,j),NN(i),i);
    end
end
else
for i=1:Length_delta_Fo-1
    for j=1:4
        error(i,j)=Teta(n(1,j),NN(i+1),i+1)-Teta(n(1,j),NN(i),i);
    end 
    
end        
end

