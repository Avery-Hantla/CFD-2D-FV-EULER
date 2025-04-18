%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to calculate the residual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res,QL_iphalf,QR_iphalf,QL_jphalf,QR_jphalf] = res(Qbar,order,islimon,BC_type,gamma,MiFaces,MjFaces,vol)
    [QL_iphalf, QR_iphalf] = reconstruction(Qbar,Qbar(:,1,:),Qbar(:,end,:),Qbar(:,2,:),Qbar(:,end-1,:),order,islimon,{BC_type{3},BC_type{4}},gamma,MiFaces,BC_type{5},2);
    F = riemann(gamma,QL_iphalf,QR_iphalf,MiFaces);
    F_iphalf = F(:,2:end,:);
    F_imhalf = F(:,1:end-1,:);

    [QL_jphalf, QR_jphalf] = reconstruction(Qbar,Qbar(:,:,1),Qbar(:,:,end),Qbar(:,:,2),Qbar(:,:,end-1),order,islimon,{BC_type{1},BC_type{2}},gamma,MjFaces,BC_type{5},3);
    F = riemann(gamma,QL_jphalf,QR_jphalf,MjFaces);
    F_jphalf = F(:,:,2:end);
    F_jmhalf = F(:,:,1:end-1);
    
    res = -1./(reshape(vol',1,length(vol(1,:)),[])).*(F_iphalf.*reshape(MiFaces(2:end,:,1),1,length(F_iphalf(1,:,1)),[])...
        -F_imhalf.*reshape(MiFaces(1:end-1,:,1),1,length(F_imhalf(1,:,1)),[])...
        +F_jphalf.*reshape(MjFaces(:,2:end,1),1,length(F_jphalf(1,:,1)),[])...
        -F_jmhalf.*reshape(MjFaces(:,1:end-1,1),1,length(F_jphalf(1,:,1)),[]));

end