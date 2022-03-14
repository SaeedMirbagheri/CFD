function [ I ] = Number_Teta( N,e )
%%%%%Calling the point return function (0.25, 0.5, 0.75 and 1)

I=zeros(1,4);

for i=1:N
    ee=e(i);
     if ee==0.25
       I(1)=i;
     elseif ee==0.5
        I(2)=i;   
     elseif ee==0.75
        I(3)=i;
     elseif ee==1
         I(4)=i;
     end
end

