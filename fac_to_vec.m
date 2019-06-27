function x = fac_to_vec(A)
N = length(A);

sz = zeros(N-2,1);
for n = 1:N-2
    sz(n) = size(A{n},1);
end
R = size(A{1},2);
P = sum(sz)*R+R+sz(2);

x = zeros(P,1);
for n = 1:N-2
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    x(idx1:idx2) = reshape(A{n},sz(n)*R,1);
end
idx1 = idx2+1;
idx2 = idx2+R;
x(idx1:idx2) = A{n+1};
idx1 = idx2+1;
idx2 = idx2+sz(2);
x(idx1:idx2) = A{n+2};
