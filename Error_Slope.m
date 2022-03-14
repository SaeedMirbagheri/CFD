function [ R_e_delta_e ] = Error_Slope( error,delta_e )
R_e=zeros(1,1);
R_delta_e=zeros(1,1);
R_e_delta_e=zeros(1,1);
for i=1:3
    for j=1:4
        R_e(i,j)=log10(error(i+1,j)/error(i,j));

    end
        R_delta_e(i)=log10(delta_e(i+1)/delta_e(i));

        
    for j=1:4
        R_e_delta_e(i,j)=abs(R_e(i,j)/R_delta_e(i));
    end
end

end

