%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Function to Find Force Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CL,CD] = findForceCoefficients(Qbar_wall,Mfaces,gamma,FreeStream_Density,FreeStream_X,FreeStream_Y,refArea)
    [~,~,~,~,P,~] = flowvariables(Qbar_wall,gamma);
    Force_X = P'.*Mfaces(:,1,1).*-Mfaces(:,1,2);
    Force_Y = P'.*Mfaces(:,1,1).*-Mfaces(:,1,3);
    CL = (-sum(Force_Y)+sum(Force_Y))./(0.5*FreeStream_Density.*sqrt(FreeStream_X.^2+FreeStream_Y.^2).^2.*refArea);
    CD = (sum(Force_X)+sum(Force_X))./(0.5*FreeStream_Density.*sqrt(FreeStream_X.^2+FreeStream_Y.^2).^2.*refArea);
end