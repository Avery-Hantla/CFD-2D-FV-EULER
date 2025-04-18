%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Function to reconstruct cells for 1st and 2nd order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QL_phalf, QR_phalf] = reconstruction(Qbar,QBC1,QBCEND,QBC2,QBCEND_m1,order,islimon,BC_type,gamma,MFaces,FreeStream,axis)
    if order == 2
        
        % Find Slope
        if islimon == true
            S_i = minmod(Qbar,QBC1,QBCEND,QBC2,QBCEND_m1,axis);
        else 
            % Shift Indexs
            Qbar_p1 = resize(Qbar,size(Qbar,axis)-2,Dimension=axis,Side="leading"); 
            Qbar_m1 = resize(Qbar,size(Qbar,axis)-2,Dimension=axis,Side="trailing");

            S_i = (Qbar_p1-Qbar_m1)./(2*dX);
            S_1 = (QBC2-QBC1)./dX;
            S_N = (QBCEND-QBCEND_m1)./dX;
            S_i = cat(axis,S_1,S_i,S_N); 
        end                

        % Reconstructs Q +-1/2
        QL_phalf = Qbar + S_i.*(0.5);%*dX);

        QR_phalf =  Qbar - S_i.*(0.5);%.*dX);
        
        QBC1 = resize(QR_phalf,1,Dimension=axis,Side="trailing");
        QBCEND = resize(QL_phalf,1,Dimension=axis,Side="leading");

        % Put boundary conditions in Q Reconstruction
        QBC1 = boundaries(QBC1,BC_type{1},FreeStream,MFaces,gamma);
        QBCEND = boundaries(QBCEND,BC_type{2},FreeStream,MFaces,gamma);
        QL_phalf = cat(axis,QBC1,QL_phalf);
        QR_phalf = cat(axis,QR_phalf,QBCEND);

    elseif order == 1
        QBC1 = boundaries(QBC1,BC_type{1},FreeStream,MFaces,gamma);
        QBCEND = boundaries(QBCEND,BC_type{2},FreeStream,MFaces,gamma);
        QL_phalf = cat(axis,QBC1,Qbar);
        QR_phalf = cat(axis,Qbar,QBCEND);
    end
end