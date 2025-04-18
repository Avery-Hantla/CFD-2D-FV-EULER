%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SSP RK2 Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q_np1,max_res,QR_jphalf] = SSP_RK2(Qbar,sigma,gamma,islimon,order,BC_type,MiFaces,MjFaces,vol)
    % Find flow variables
    [dQdt,QL_iphalf,QR_iphalf,QL_jphalf,QR_jphalf] = res(Qbar,order,islimon,BC_type,gamma,MiFaces,MjFaces,vol);
    max_res = max(dQdt,[],2:3);

    % Find current iteration time step
    % [QL_iphalf, QR_iphalf] = reconstruction(Qbar,Qbar(:,1,:),Qbar(:,end,:),Qbar(:,2,:),Qbar(:,end-1,:),order,islimon,{BC_type{3},BC_type{4}},gamma,MiFaces,BC_type{5},2);
    % [QL_jphalf, QR_jphalf] = reconstruction(Qbar,Qbar(:,:,1),Qbar(:,:,end),Qbar(:,:,2),Qbar(:,:,end-1),order,islimon,{BC_type{1},BC_type{2}},gamma,MjFaces,BC_type{5},3);
    dt = finddt(sigma,gamma,MiFaces,MjFaces,vol,QL_iphalf,QR_iphalf,QL_jphalf,QR_jphalf);

    % Calculate Qstar
    Q_np1 = Qbar + dt.*dQdt;

    % Find Q for next time step
    Q_np1 = 0.5.*(Qbar+Q_np1+dt.*res(Q_np1,order,islimon,BC_type,gamma,MiFaces,MjFaces,vol));

end

