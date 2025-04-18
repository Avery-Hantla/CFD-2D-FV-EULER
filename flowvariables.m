%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Function to Calculate Flow Variables from Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho,u,v,E,P,c] = flowvariables(Q,gamma)
    rho = Q(1,:,:);
    u = Q(2,:,:)./rho;
    v = Q(3,:,:)./rho;
    E = Q(4,:,:);
    P = (E-0.5.*rho.*(u.^2+v.^2)).*(gamma-1);
    c = sqrt(gamma.*(P./rho));

    % rho = reshape(rho,length(rho(1,:,1)),length(rho(1,1,:)));
    % u = reshape(u,length(u(1,:,1)),length(u(1,1,:)));
    % v = reshape(v,length(v(1,:,1)),length(v(1,1,:)));
    % E = reshape(E,length(E(1,:,1)),length(E(1,1,:)));
    % P = reshape(P,length(P(1,:,1)),length(P(1,1,:)));
    % c = reshape(c,length(c(1,:,1)),length(c(1,1,:)));
end