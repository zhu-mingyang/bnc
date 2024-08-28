function [B_lift,B_lift_QC] = Lifting_Protograph(B,Z1,Z2)

if (Z1 < max(max(B)))
    error('Z1 is too small to eliminate all parallel edges.')
end

[M,N] = size(B);
B1 = zeros(M*Z1,N*Z1);
for i = 1:M
    for j = 1:N
        row = zeros(1,Z1);
        row(randperm(Z1,B(i,j))) = 1;
        for k = 1:Z1
            B1((i-1)*Z1+k,(j-1)*Z1+1:j*Z1) = circshift(row,k-1);
        end
    end
end

B2 = B1;
for i = 1:M*Z1
    for j = 1:N*Z1
        if (B2(i,j) == 0)
            B2(i,j) = -1;
        else
            B2(i,j) = randi([0,Z2-1],1,1);
        end
    end
end

B_lift_QC = B2;
B_lift = zeros(M*Z1*Z2,N*Z1*Z2);

for i = 1:M*Z1
    for j = 1:N*Z1
        if (B_lift_QC(i,j) >= 0)
            B_lift((i-1)*Z2+1:i*Z2,(j-1)*Z2+1:j*Z2) = circshift(eye(Z2),B_lift_QC(i,j),2);
        end
    end
end

end