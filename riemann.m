%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Function that calulcates the rosuvinov flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = riemann(gamma,QL_phalf,QR_phalf,MFaces) %% Validated
    % Shift Indexs
    % [QL_iphalf, QR_iphalf] = reconstruction(Qbar,dX,order,islimon,QBC);

    % Find value for left side
    [rho,u,v,E,P,cL]  = flowvariables(QL_phalf,gamma); 
    nx = reshape(MFaces(:,:,2),1,length(MFaces(:,1,1)),[]);
    ny = reshape(MFaces(:,:,3),1,length(MFaces(:,1,1)),[]);
    % vnL = [u,v].*[MFaces(:,:,2),MFaces(:,:,3)];
    vnL = u.*nx+v.*ny;
    FL_phalf = cat(1,reshape(rho.*vnL,[1,length(QL_phalf(1,:,1)),length(QL_phalf(1,1,:))]),...
        rho.*u.*vnL+P.*nx,...
        rho.*v.*vnL+P.*ny,...
        reshape(vnL.*(E+P),[1,length(QL_phalf(1,:,1)),length(QL_phalf(1,1,:))]));
    % FL_phalf = [rho.*u;rho.*u.^2+P;rho.*u.*v;u.*(E+P)].*MFaces(:,:,2) + [rho.*v;rho.*v.*u;rho.*v.^2+P;v.*(E+P)].*MFaces(:,:,3);

    % Find values for right side
    [rho,u,v,E,P,cR]  = flowvariables(QR_phalf,gamma); 
    % vnR = [u,v].*[MFaces(:,:,2),MFaces(:,:,3)];
    vnR = u.*nx+v.*ny;
    FR_phalf = cat(1,reshape(rho.*vnR,[1,length(QR_phalf(1,:,1)),length(QR_phalf(1,1,:))]),...
        rho.*u.*vnR+P.*nx,...
        rho.*v.*vnR+P.*ny,...=
        reshape(vnR.*(E+P),[1,length(QR_phalf(1,:,1)),length(QR_phalf(1,1,:))]));
    % FR_phalf = [rho.*vnR;rho.*u.*vnR+P.*MFaces(:,:,2);rho.*v.*vnR+P.*MFaces(:,:,3);vnR.*(E+P)];
    % FR_phalf = [rho.*vnR;rho.*u.*vnR+P.*nx;rho.*v.*vnR+P.*ny;vnR.*(E+P)];

    % Find average u and c
    vn = (abs(vnL)+abs(vnR))./2;
    c = (cL+cR)./2;

    % Solve for F i plus half
    F = ((FL_phalf+FR_phalf)./2) - bsxfun(@times, (QR_phalf-QL_phalf), reshape(((vn+c)./2),[1,length(QR_phalf(1,:,1)),length(QR_phalf(1,1,:))]));  %((vn+c)./2).*(QR_phalf-QL_phalf);
    % F_phalf = F(4,:,2:end);
    % F_mhalf = F(4,:,1:end-1);
end
