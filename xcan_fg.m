function [f,g] = xcan_fg(Z, A, pen1, lambda1, XtX, lambda2, XXt, withm)
% Data structure:
% Z: data set
% lambda: regularization coeficient
% XX: structural correlation, can be the correlation matrix, MEDA or similar. 

% Data setup 
[m,n] = size(Z); 
U = A{1};
V = A{2};
D = A{3};
if withm, v0 = A{4}; end;
k = size(A{1},2);
% Function value (residual) 
if withm,
    E = (Z - ones(size(Z,1),1)*v0' - U * diag(D) * V'); 
else
    E = (Z - U * diag(D) * V'); 
end

f0 = 0;
for a=1:k,
    f0 = f0 + (V(:,a)'*V(:,a)-1)^2 + (U(:,a)'*U(:,a)-1)^2;
end

f1 = 0;
for a=1:k,
    f1 = f1 + sum(sum((V(:,a)*V(:,a)'./XtX).^2));
end

% f1 = 0;
% for j1=1:n,
%     for j2=1:n,
%         for a=1:k,
%             f1 = f1 + (V(j1,a)*V(j2,a)/XtX(j1,j2))^2;
%         end
%     end
% end

f2 = 0;
for a=1:k,
    f2 = f2 + sum(sum((U(:,a)*U(:,a)'./XXt).^2));
end

% f2 = 0;
% for i1=1:m,
%     for i2=1:m,
%         for a=1:k,
%             f2 = f2 + (U(i1,a)*U(i2,a)/XXt(i1,i2))^2;
%         end
%     end
% end

f = norm(E, 'fro')^2 + pen1*f0 + lambda1*f1  + lambda2*f2;

% First derivatives computed in matrix form 

for a=1:k,
   g0(:,a) =  (V(:,a)'*V(:,a)-1)*V(:,a); 
   g0b(:,a) =  (U(:,a)'*U(:,a)-1)*U(:,a); 
end

% for j=1:n,
%    for a=1:k,
%        g0(j,a) =  (V(:,a)'*V(:,a)-1)*V(j,a);   
%    end
% end

% for i=1:m,
%    for a=1:k,
%        g0b(i,a) =  (U(:,a)'*U(:,a)-1)*U(i,a);   
%    end
% end

for j=1:n,
   for a=1:k,
       g1(j,a) =  V(j,a)*(V(:,a)./XtX(:,j))'*(V(:,a)./XtX(:,j));   
   end
end

for i=1:m,
   for a=1:k,
       g2(i,a) =  U(i,a)*(U(:,a)./XXt(:,i))'*(U(:,a)./XXt(:,i));   
   end
end

g = zeros((m+n)*k+k+n,1); 

g(1:m*k) = reshape(-2 * E * V * diag(D)' + 4 * pen1 * g0b + 4 * lambda2 * g2, m*k, 1);

g(m*k+1:(m+n)*k) = reshape(-2 * E' * U * diag(D) + 4 * pen1 * g0 + 4 * lambda1 * g1, n*k, 1);

g((m+n)*k+1:(m+n)*k+k) = reshape(-2 * diag(U' * E * V), k, 1);

if withm,
    g((m+n)*k+k+1:end) = reshape(-2 * E' * ones(size(Z,1),1), size(E,2), 1);
end
