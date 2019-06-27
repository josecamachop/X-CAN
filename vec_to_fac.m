function A = cp_vec_to_fac(x,Z)

P = length(x);
N = ndims(Z);
sz = size(Z);

R = (P - sz(2)) / (sum(sz)+1);

A = cell(N,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end

A{3} = x(idx2+1:idx2+R);
A{4} = x(idx2+R+1:end);