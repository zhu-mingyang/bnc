function [A,G,H,Y] = Construct_Precoded_Batches(PCM,row_index,q,T,N_INAC)
Y = gf(zeros(1+N_INAC,1),log2(q));
dg = sum(PCM(row_index,:)~=0);
A = find(PCM(row_index,:)~=0);
G = PCM(row_index,A)';
H = gf(1,log2(q));
Y(1,1) = gf(zeros(T,1),log2(q));
end