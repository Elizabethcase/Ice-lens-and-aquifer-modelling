function [fadvp,fadvm]=AdvectiveFlux(N,velocity,field)
% velocity is (N+1)x1
% field is Nx1

fadv = zeros(N+1,1);
for i = 1:(N+1)
    if and(velocity(i)>0,i==1)
        % fadv(i) = 0; % will be fixed by the boundary condition
        fadv(i) = velocity(i)*field(i); % will usually be fixed by the boundary condition
    elseif and(velocity(i)>0,i~=1)
        fadv(i) = velocity(i)*field(i-1);
    elseif and(velocity(i)<=0,i==(N+1))
        % fadv(i) = 0; % will be fixed by the boundary condition
        fadv(i) = velocity(i)*field(i-1); % natural inflow
    elseif and(velocity(i)<=0,i~=(N+1)) % i.e velocity(i)<=0 and not at the end point
        fadv(i) = velocity(i)*field(i);
    end
end
fadvp = fadv(2:(N+1));
fadvm = fadv(1:N);

% vp = velocity>0;
% vm = velocity<0;
% 
% V = [zeros(N-1,1)', vp', vm(1:(N-1))'];
% I = [2:N, 1:N, 1:(N-1)];
% J = [1:(N-1), 1:N, 2:N];
% Dp = sparse(I,J,V,N,N);
% 
% V = [vp(1:(N-1))', [0; vm(1:(N-1))]', zeros(N-1,1)'];
% I = [2:N, 1:N, 1:(N-1)];
% J = [1:(N-1), 1:N, 2:N];
% Dm = sparse(I,J,V,N,N);

% V = [-vp(2:N)', (vp-vm)', vm(1:(N-1))'];
% I = [2:N, 1:N, 1:(N-1)];
% J = [1:(N-1), 1:N, 2:N];
% D = sparse(I,J,V,N,N);
% D2 = Dp-Dm;
end