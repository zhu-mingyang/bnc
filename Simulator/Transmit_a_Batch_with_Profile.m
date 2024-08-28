function [A,G,H,Y] = Transmit_a_Batch_with_Profile(B,M,q,k,e,N_INAC,Profile)
Y = gf(zeros(1+N_INAC,M),log2(q));
m = log2(q);
A = Profile;                                                         % the index of contributors to this batch. randperm(n,k) generates k distinct integers in 1~n
dg = size(A,2);
G = randi([0 q-1],dg,M);
G = gf(G,m);
X = B(:,A)*G;
H = Erasure_Mat(M,e(1));                                                       % H_1 = E_1
H = gf(H,m);
for i = 2:k
    E = Erasure_Mat(M,e(i));
    Phi = randi([0 q-1],M,M);
    Phi = gf(Phi,m);
    H = H*Phi*E;                                                            % H_k = H_k-1*Phi*E
end
Y(1,:) = X*H;
end