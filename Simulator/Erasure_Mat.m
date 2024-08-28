function [E] = Erasure_Mat(M,e)
E = randsrc(M,M,[0,1; e,1-e]);
i = logical(1-eye(M));
E(i) = 0;
end