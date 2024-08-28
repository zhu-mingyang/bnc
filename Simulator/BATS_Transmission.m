function [A,G,H,Y] = BATS_Transmission(B,M,q,Psi,k,e,N_INAC,H_reg,T,n,m,BATS_Type,BATS_profile)
A = cell(1,n+m);
G = cell(1,n+m);
H = cell(1,n+m);
Y = cell(1,n+m);
if (strcmp(BATS_Type,'random'))
    for j = 1:n                                                             % transmit n batches
        [A{j},G{j},H{j},Y{j}] = Transmit_a_Batch(B,M,q,Psi,k,e,N_INAC);
    end
else
    for j = 1:n
        [A{j},G{j},H{j},Y{j}] = Transmit_a_Batch_with_Profile(B,M,q,k,e,N_INAC,BATS_profile{j});
    end
end
for j = n+1:n+m
    [A{j},G{j},H{j},Y{j}] = Construct_Precoded_Batches(H_reg,j-n,q,T,N_INAC);
end
end