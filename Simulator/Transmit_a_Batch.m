function [A,G,H,Y] = Transmit_a_Batch(B,M,q,Psi,k,e,N_INAC)
Y = gf(zeros(1+N_INAC,M),log2(q));
[T,K] = size(B);
m = log2(q);
dg = randsrc(1,1,[1:K; Psi]);
A = randperm(K,dg);                                                         % the index of contributors to this batch. randperm(n,k) generates k distinct integers in 1~n
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