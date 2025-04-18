%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to Compute Time Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = finddt(sigma,gamma,MiFaces,MjFaces,vol,QL_iphalf,QR_iphalf,QL_jphalf,QR_jphalf)
    % [QL_iphalf, QR_iphalf] = reconstruction(Qbar,Qbar(:,1,:),Qbar(:,end,:),Qbar(:,2,:),Qbar(:,end-1,:),order,islimon,{BC_type{3},BC_type{4}},gamma,MiFaces,BC_type{5},2);
    % [QL_jphalf, QR_jphalf] = reconstruction(Qbar,Qbar(:,:,1),Qbar(:,:,end),Qbar(:,:,2),Qbar(:,:,end-1),order,islimon,{BC_type{1},BC_type{2}},gamma,MjFaces,BC_type{5},3);

    Q_i = (QL_iphalf+QR_iphalf)./2;
    Q_j = (QL_jphalf+QR_jphalf)./2;

    % Find flow variables
    [~,u_i,v_i,~,~,c_i]  = flowvariables(Q_i,gamma); 
    [~,u_j,v_j,~,~,c_j]  = flowvariables(Q_j,gamma); 

    u_i = reshape(u_i,length(u_i(1,:,1)),[]);
    v_i = reshape(v_i,length(v_i(1,:,1)),[]);
    c_i = reshape(c_i,length(c_i(1,:,1)),[]);

    u_j = reshape(u_j,length(u_j(1,:,1)),[]);
    v_j = reshape(v_j,length(v_j(1,:,1)),[]);
    c_j = reshape(c_j,length(c_j(1,:,1)),[]);

    % Find Face Area
    Si = MiFaces(:,:,1);
    Sj = MjFaces(:,:,1);

    % Find nomral vectors
    nxi = MiFaces(:,:,2);
    nyi = MiFaces(:,:,3);

    nxj = MjFaces(:,:,2);
    nyj = MjFaces(:,:,3);

    % Find normal velocities
    vni = u_i.*nxi+v_i.*nyi;
    vnj = u_j.*nxj+v_j.*nyj;

    % Compute time step
    left_face = (abs(vni(1:end-1,:))+c_i(1:end-1,:)).*Si(1:end-1,:);
    right_face = (abs(vni(2:end,:))+c_i(2:end,:)).*Si(2:end,:);
    top_face = (abs(vnj(:,2:end))+c_j(:,2:end)).*Sj(:,2:end);
    bottom_face = (abs(vnj(:,1:end-1))+c_j(:,1:end-1)).*Sj(:,1:end-1);

    dt = (2*sigma*vol')./(left_face+right_face+top_face+bottom_face);
    dt = min(dt,[],'all');
end