%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Function to Find Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function QBC = boundaries(Qbar,BC_type,FreeStream,MFaces,gamma)
    if strcmp(BC_type,'Symmetry') %% Validated 
        [rho,u,v,E,~,~] = flowvariables(Qbar,gamma);
        v = -v;
        QBC = [rho;rho.*u;rho.*v;E];
        QBC = reshape(QBC,length(QBC(:,1)),1,length(QBC(1,:)));
    elseif strcmp(BC_type,'Extrapolation') 
        QBC = Qbar;
    elseif strcmp(BC_type,'Slip_Wall') %% Validated 
        [rho,u,v,E,~,~] = flowvariables(Qbar,gamma);
        nx = reshape(MFaces(:,1,2),[],length(MFaces(:,1,1)));
        ny = reshape(MFaces(:,1,3),[],length(MFaces(:,1,1)));
        vn = u.*nx+v.*ny;
        u = u-2.*vn.*nx;
        v = v-2.*vn.*ny;
        QBC = [rho;rho.*u;rho.*v;E];
    else % Freestream %% Validated
        QBC = FreeStream.*ones(4,length(Qbar(1,:)));
    end
end