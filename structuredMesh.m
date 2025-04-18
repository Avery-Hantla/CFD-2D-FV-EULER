%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Function to Create Structured Grid
%
%                              Domain 2
%                             __________
%                            |          |
%                   Domain 3 |          | Domain 4
%                            |          |
%                            |__________|
%              ^               Domain 1
%             eta
%                 xi >
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grid,MIFaces,MJFaces,Vol,cell_centers] = structuredMesh(ID,JD,Domain_1,Domain_2,Domain_3,Domain_4)
    eta_points = JD+1;
    xi_points = ID+1;
    % plot(Domain_1(1,:),Domain_1(2,:)); hold on
    % plot(Domain_2(1,:),Domain_2(2,:));
    % plot(Domain_3(1,:),Domain_3(2,:));
    % plot(Domain_4(1,:),Domain_4(2,:));

    % Discretize Domain 1 and 2
    s_12 = linspace(0,1,xi_points);
    Domain_1_dis(1,:) = pdearcl(Domain_1(1,:),Domain_1,s_12,0,1);
    Domain_1_dis(2,:) = pdearcl(Domain_1(2,:),Domain_1,s_12,0,1);
    
    Domain_2_dis(1,:) = pdearcl(Domain_2(1,:),Domain_2,s_12,0,1);
    Domain_2_dis(2,:) = pdearcl(Domain_2(2,:),Domain_2,s_12,0,1);
    
    % Interpolate for each mesh point assuming Domains 3 and 4 are linear
    s_34 = linspace(0,1,eta_points);
    grid_star = zeros(eta_points,xi_points,2);
    for idx = 1:xi_points
        lin_line = straight_line([Domain_1_dis(1,idx),Domain_1_dis(2,idx)],[Domain_2_dis(1,idx),Domain_2_dis(2,idx)]);
        grid_star(:,idx,1) = pdearcl(lin_line(1,:),lin_line,s_34,0,1); %linspace(Domain_1_dis(1,idx),Domain_2_dis(1,idx),eta_points);
        grid_star(:,idx,2) = pdearcl(lin_line(2,:),lin_line,s_34,0,1); %linspace(Domain_1_dis(2,idx),Domain_2_dis(2,idx),eta_points);
    end
    
    % Linear Interpolation Correction for Domains 3 and 4 
    % Discretize Domain 3 and 4
    Domain_3_dis(1,:) = pdearcl(Domain_3(1,:),Domain_3,s_34,0,1);
    Domain_3_dis(2,:) = pdearcl(Domain_3(2,:),Domain_3,s_34,0,1);
    
    Domain_4_dis(1,:) = pdearcl(Domain_4(1,:),Domain_4,s_34,0,1);
    Domain_4_dis(2,:) = pdearcl(Domain_4(2,:),Domain_4,s_34,0,1);
    
    % Find error from grid 1 from domain 3 and 4
    e3 = Domain_3_dis - [grid_star(:,1,1)';grid_star(:,1,2)'];
    e4 = Domain_4_dis - [grid_star(:,end,1)';grid_star(:,end,2)'];
    
    % Interpolate and correct error
    e_correct = zeros(eta_points,xi_points,2);
    for jdx = 1:eta_points
        e_correct(jdx,:,1) = linspace(e3(1,jdx),e4(1,jdx),xi_points);
        e_correct(jdx,:,2) = linspace(e3(2,jdx),e4(2,jdx),xi_points);
    end
    
    grid = grid_star+e_correct;
    % 
    % scatter(grid(:,:,1),grid(:,:,2),5,'black','filled'); axis equal
    % xlabel('x'); ylabel('y')

    % Find volumes, area of faces, and normal vectors
    grid_jp1 = grid(:,2:end,:);
    grid_ip1 = grid(2:end,:,:);
    
    % Find Area of Faces
    areaJFaces = sqrt((grid_jp1(:,:,1)-grid(:,1:end-1,1)).^2+(grid_jp1(:,:,2)-grid(:,1:end-1,2)).^2);
    areaIFaces = sqrt((grid_ip1(:,:,1)-grid(1:end-1,:,1)).^2+(grid_ip1(:,:,2)-grid(1:end-1,:,2)).^2);

    % Find normal Vector of Faces
    nxIFaces = (grid(2:end,:,2)-grid(1:end-1,:,2))./areaIFaces;
    nyIFaces = (grid(1:end-1,:,1)-grid(2:end,:,1))./areaIFaces;

    nxJFaces = -(grid(:,2:end,2)-grid(:,1:end-1,2))./areaJFaces;
    nyJFaces = -(grid(:,1:end-1,1)-grid(:,2:end,1))./areaJFaces; 

    % quiver(grid(1:end-1,:,1),grid(1:end-1,:,2),nxIFaces,nyIFaces,0.5,'Color','black')
    % quiver(grid(:,1:end-1,1),grid(:,1:end-1,2),nxJFaces,nyJFaces,0.5,'Color','black')


    % for j = 1:20
    %     for i = 1:40 
    %         test = cross([grid(j,i+1,1)-grid(j+1,i,1),grid(j,i+1,2)-grid(j+1,i,2),1 ],[grid(j+1,i+1,1)-grid(j,i,1), grid(j+1,i+1,2)-grid(j,i,2),1]);
    %         Vol(j,i) = norm(test);
    %     end
    % end
    u1 = grid(1:end-1,2:end,1)-grid(2:end,1:end-1,1);
    v1 = grid(1:end-1,2:end,2)-grid(2:end,1:end-1,2);
    u2 = grid(2:end,2:end,1)-grid(1:end-1,1:end-1,1);
    v2 = grid(2:end,2:end,2)-grid(1:end-1,1:end-1,2);

    Vol = 0.5*abs(u1.*v2-u2.*v1);
    % 
    % % Find Normal Vector of Faces
    % normJFaces = atan2(grid_jp1(:,:,2)-grid(:,1:end-1,2),grid_jp1(:,:,1)-grid(:,1:end-1,1))+(pi/2);
    % normIFaces = atan2(grid_ip1(:,:,2)-grid(1:end-1,:,2),grid_ip1(:,:,1)-grid(1:end-1,:,1))-(pi/2);
    % 
    % nxIFaces = cos(normIFaces); nxIFaces(isinf(normIFaces)) = 0;
    % nyIFaces = sin(normIFaces); nyIFaces(isinf(normIFaces)) = 1;
    % nxJFaces = cos(normJFaces); nxJFaces(isinf(normJFaces)) = 0;
    % nyJFaces = sin(normJFaces); nyJFaces(isinf(normJFaces)) = 1;
    % 
    % Concentrate Matrix together
    MIFaces = cat(3,areaIFaces',nxIFaces',nyIFaces');
    MJFaces = cat(3,areaJFaces',nxJFaces',nyJFaces');
    % 
    % % Find volume of Cells
    % areaIFaces_ip1 = areaIFaces(:,2:end);
    % areaIFaces_i = areaIFaces(:,1:end-1);
    % areaJFaces_jp1 = areaJFaces(2:end,:);
    % areaJFaces_j = areaJFaces(1:end-1,:);
    % 
    % IFaces_avg = (areaIFaces_ip1+areaIFaces_i)/2;
    % JFaces_avg = (areaJFaces_jp1+areaJFaces_j)/2;
    % 
    % Vol = IFaces_avg.*JFaces_avg;
    % 
    % % Find cell centers
    cell_centers_i = (grid(:,2:end,:)+grid(:,1:end-1,:))/2;
    % cell_centers_j = (grid(2:end,:,:)+grid(1:end-1,:,:))/2;

    cell_centers = (cell_centers_i(2:end,:,:)+cell_centers_i(1:end-1,:,:))/2;
    % cell_centers = (cell_centers_j(:,2:end,:)+cell_centers_j(:,1:end-1,:))/2;
    cell_centers = permute(cell_centers,[2,1,3]);
    % scatter(cell_centers(:,:,1),cell_centers(:,:,2),5,'red','filled'); axis equal
    % xlabel('x'); ylabel('y')

end