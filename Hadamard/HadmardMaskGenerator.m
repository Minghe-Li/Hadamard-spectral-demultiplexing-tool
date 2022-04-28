%Sylvester-type Hadamard matrix generator for QCL spectral analysis
%Written by Garth J. Simpson, 2020

clearvars;

Dims = 32; %size of the desired square Hadamard mask.

H_dims = 48;
temp = hadamard(H_dims );%note - size must be greater than Dims.
for r = 1:H_dims
    for c = 1:H_dims
        if temp(r,c)<0
            temp(r,c)=0;
        end
    end
end

H_Mat = temp(2:Dims+1,2:Dims+1); %removes first row/column of all 1's. 

H_Inv = 8*inv(H_Mat);