 n = 300;
% p = 5*(log(n)/n);
% q = 1*(log(n)/n);

%[RP,overlap, comm1,comm2] = permutate(p,q,n)

% overlap_vector = zeros(20,1);
% score = zeros(65,50);
% for a = 5:70
%     p = a/n;
%     for b = 1:50
%         q = b/n;
%         for i = 1:20 
%             overlap = permutate(p,q,n);
%             overlap_vector(i) = overlap;
%         end
%         score(a-4,b) = mean(overlap_vector);
%     end
% end
           
load /Users/Adam/Downloads/zachary.mat;
%Finding the two communities with omegat_tilda
[C1,C2, Omega_Tilda, Neg_Omega_Tilda] = partition(34,A);

%Finding the true value of omega
omega = zeros(34,1);
for i = 1:34
    if i <= 34/2
        omega(i) = 1;
    else 
        omega(i) = -1;
    end
end 
%Computing Delta
for i = 1:34
    if Omega_Tilda(i) == omega(i)
        delta(i) = 1; 
    else 
        delta(i) = 0;
    end
end
%Computing Negative Delta
for i = 1:34
    if Neg_Omega_Tilda(i) == omega(i)
        Neg_delta(i) = 1; 
    else 
        Neg_delta(i) = 0;
    end
end
%Computing the over lap%
rawoverlap = max(sum(delta),sum(Neg_delta)); 
overlap = (2/34)*rawoverlap-1

function [overlap, comm1, comm2] = permutate(p,q,n)
% Genertation of a random permutation matrix T
I = eye(n);
ix = randperm (n);
T = I(ix,:);


% Generation of the adjacency matrix A of a stochastic block model defined
%over n nodes
n2 = n/2;
P = random('bino', 1, p, n2, n2); % upper left block
dP2 = random('bino', 1, p, n2, 1); % diagonal of the lower right
%block
Q = random('bino', 1, q, n2, n2); % upper right block
% carve the two triangular and diagonal matrices that we need
U = triu(P, 1);
L = tril(P,-1);
dP = diag(P);
A0 = U + U' + diag(dP);
A1 = Q;
A2 = Q';
A3 = L + L' + diag(dP2);
A =[A0 A1;A2 A3];

RP = T*A*T.';

omega = zeros(n,1);
for i = 1:n
    if i <= n/2
        omega(i) = 1;
    else 
        omega(i) = -1;
    end
end 

[V,l] = eigs(RP,2);
v2 = V(3:end);
comm1 = zeros(n,1);
comm2 = zeros(n,1);
omega_hat = zeros(n,1);
for i = 1:n
    
    if v2(i) > 0
        comm1(i) = i;
        omega_hat(i) = 1;
        neg_omega_hat(i) = -1;
    else 
        comm2(i) = i;
        omega_hat(i) = -1;
        neg_omega_hat(i) = 1;
    end  
end

%Computing Delta
for i = 1:n
    if omega_hat(i) == omega(i)
        delta(i) = 1; 
    else 
        delta(i) = 0;
    end
end
%Computing Negative Delta
for i = 1:n
    if neg_omega_hat(i) == omega(i)
        Neg_delta(i) = 1; 
    else 
        Neg_delta(i) = 0;
    end
end
%Computing the over lap%
rawoverlap = max(sum(delta),sum(Neg_delta)); 
overlap = (2/n)*rawoverlap-1;
end 

function [comm1,comm2, omega_hat, neg_omega_hat] = partition(n,A)
[V,lamb] = eigs(A,2);
v2 = V(3:end);
comm1 = zeros(n,1);
comm2 = zeros(n,1);
for i = 1:n
    
    if v2(i) > 0
        comm1(i) = i;
        omega_hat(i) = 1;
        neg_omega_hat(i) = -1;
    else 
        comm2(i) = i;
        omega_hat(i) = -1;
        neg_omega_hat(i) = 1;
    end  
end
end



