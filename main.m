%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           MATLAB 2D Euler Code
%                               Avery Hantla
%                              November, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID = 160;                           % Number of Cells in I direction
JD = 80;                            % Number of Cells in J direction
sigma = 0.75;                       % CFL Number                 
order = 1;                          % Desired Order of Error
islimiteron = true;                 % Use limiter? true/false
n = 5000;                          % Number of Iterations
isplot = true;                     % Plot during sim? true/false
output_step = 50;                  % How many iterations between plots     
restart = false;                    % Restart flag
restart_file = 'second_order_coarse.mat';  % Restart File Name

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain 1
Domain_1_1 = fliplr(circle([0,0],1,120,180));
Domain_1_2 = straight_line([Domain_1_1(1,end),Domain_1_1(2,end)],[1.2320508,1.8660254]);
Domain_1 = [Domain_1_1(:,1:end-1),Domain_1_2];
Domain_1_type = "Slip_Wall";

% Domain 2
Domain_2_1 = fliplr(circle([0,0],2,150,180));
Domain_2_2 = straight_line([Domain_2_1(1,end),Domain_2_1(2,end)],[0,4]);
Domain_2 = [Domain_2_1(:,1:end-1),Domain_2_2];
Domain_2_type = "FreeStream";
FreeStream_Pressure = 1;
FreeStream_Density = 1.4;
FreeStream_X = 5;
FreeStream_Y = 0;
gamma = 1.4;
refArea = 3.3013*2;

% Domain 3
Domain_3 = straight_line([-1,0],[-2,0]);
Domain_3_type = "Symmetry";

% Domain 4
Domain_4 = straight_line([1.2320508,1.8660254],[0,4]);
Domain_4_type = "Extrapolation";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initilize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[grid,MiFaces,MjFaces,vol,cell_centers] = structuredMesh(ID,JD,Domain_1,Domain_2,Domain_3,Domain_4);
FreeStream_Energy = FreeStream_Pressure/(gamma-1)+0.5*(FreeStream_Density*(FreeStream_X^2+FreeStream_Y^2));
FreeStream_Q = [FreeStream_Density;FreeStream_Density*FreeStream_X;FreeStream_Density*FreeStream_Y;FreeStream_Energy];
BC_type = {Domain_1_type,Domain_2_type,Domain_3_type,Domain_4_type,FreeStream_Q};
Wall_Boundary = find(strcmp(cellstr(BC_type(1:4)),"Slip_Wall")==1);

res = zeros(4,n);
CL = zeros(1,n);
CD = zeros(1,n);

if restart == true 
    try saved_solution = load(['sol/',restart_file]);
    catch; disp('ERROR: No saved solution'); return; end

    if saved_solution.ID ~= ID || saved_solution.JD ~= JD 
            disp('ERROR: Saved solution conflicts with inputted grid dimentions'); return 
    end
    n_start = saved_solution.iteration;
    Qbar = saved_solution.Q;
    res(:,1:n_start) = saved_solution.res;
    CD(1:n_start) = saved_solution.CD;
    CL(1:n_start) = saved_solution.CL;
else 
    Qbar = FreeStream_Q.*ones(4,ID,JD); 
    n_start = 1;
end

if isplot == true
    [~,u,v,~,~,c] = flowvariables(Qbar,gamma);
    Mach = sqrt(u.^2+v.^2)./c; Mach(1,1,1) = 0;
    figure; cMap=jet(256);
    [~,h] = contourf(grid(1:end-1,1:end-1,1)',grid(1:end-1,1:end-1,2)',reshape(Mach,ID,JD),1000);
    set(h, 'edgecolor','none');
    colormap(cMap); cb=colorbar; ylabel(cb,'Mach (~)','FontSize',12,'Rotation',270); axis equal
    xlabel('X Coordinate'); ylabel('Y Coordinate');

    M(n/output_step) = struct('cdata',[],'colormap',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Sim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = n_start:n % Time Intergration Loop
    [Qbar,res(:,n),QR_jphalf] = SSP_RK2(Qbar,sigma,gamma,islimiteron,order,BC_type,MiFaces,MjFaces,vol);
    [CL(n),CD(n)] = findForceCoefficients(QR_jphalf(:,:,1),MjFaces(:,1,:),gamma,FreeStream_Density,FreeStream_X,FreeStream_Y,refArea);

    if mod(n,output_step) == 0
    fprintf('n: %d, res: %d, %d, %d, %d \n',n, (res(1,n)), (res(2,n)), (res(3,n)), (res(4,n)))
        if isplot == true
            [~,u,v,~,~,c] = flowvariables(Qbar,gamma);
            Mach = sqrt(u.^2+v.^2)./c;
            figure(1)
            h.ZData = reshape(Mach,ID,JD);

            M(n/output_step) = getframe(gcf);
        end
    end
end

% Save Solution
save_struct = struct('iteration',n,'Q',Qbar,'ID',ID,'JD',JD,'res',res,'CD',CD,'CL',CL,'grid',grid);
save(['sol/',restart_file],"-struct","save_struct");

% Save Video
if isplot == true
    movie(figure,M);
    v = VideoWriter("sol/Mach_History");
    open(v)
    writeVideo(v,M)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% Post Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Residuals
figure; semilogy(1:n,res(1,:),1:n,res(2,:),1:n,res(3,:),1:n,res(4,:)); 
legend('Continuity','X Momentum','Y Momentum','Energy'); xlabel('Iterations'); ylabel('Residual')

% % Plot CL
% figure; plot(1:n,CL); xlabel('Iterations'); ylabel('C_L')
% % Plot CD 
% figure; plot(1:n,CD); xlabel('Iterations'); ylabel('C_D')

[rho,u,v,E,P,c] = flowvariables(Qbar,gamma);
Mach = sqrt(u.^2+v.^2)./c;
Mach_cell = griddata(cell_centers(:,:,1),cell_centers(:,:,2),reshape(Mach,ID,JD),grid(:,:,1),grid(:,:,2),"v4");
figure; cMap=jet(256);
[~,h] = contourf(grid(:,:,1),grid(:,:,2),Mach_cell,1000); hold on
set(h, 'edgecolor','none');
colormap(cMap); cb=colorbar; ylabel(cb,'Mach Number (~)','FontSize',12,'Rotation',270); axis equal
xlabel('X Coordinate'); ylabel('Y Coordinate');

% Plot Mach Contour
[rho,u,v,E,P,c] = flowvariables(Qbar,gamma);
Mach = sqrt(u.^2+v.^2)./c;
figure; cMap=jet(256);
[~,h] = contourf(grid(1:end-1,1:end-1,1)',grid(1:end-1,1:end-1,2)',reshape(Mach,ID,JD),1000); hold on
set(h, 'edgecolor','none');
colormap(cMap); cb=colorbar; ylabel(cb,'Mach Number (~)','FontSize',12,'Rotation',270); axis equal
xlabel('X Coordinate'); ylabel('Y Coordinate');

% figure;
% quiver(grid(1:end-1,1:end-1,1)',grid(1:end-1,1:end-1,2)',reshape(u,ID,JD),reshape(v,ID,JD),0.75,'Color','black')
% streamline(grid(1:end-1,1:end-1,1)',grid(1:end-1,1:end-1,2)',reshape(u,ID,JD),reshape(v,ID,JD))

% Plot Pressure Contour
Cp = (P-FreeStream_Pressure)./(0.5*FreeStream_Density.*sqrt(FreeStream_X.^2+FreeStream_Y.^2).^2);
figure; cMap=jet(256);
[~,h] = contourf(grid(1:end-1,1:end-1,1)',grid(1:end-1,1:end-1,2)',reshape(Cp,ID,JD),1000);
set(h, 'edgecolor','none');
colormap(cMap); cb=colorbar; ylabel(cb,'Cofficient of Pressure, Cp (~)','FontSize',12,'Rotation',270); axis equal
xlabel('X Coordinate'); ylabel('Y Coordinate');
