%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Function to limit the 2nd order reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_i = minmod(Qbar,QBC1,QBCEND,QBC2,QBCEND_m1,axis)
    % % Obtain Q p1 and m1 
    Qbar_p1 = resize(Qbar,size(Qbar,axis)-2,Dimension=axis,Side="leading");
    Qbar_m1 = resize(Qbar,size(Qbar,axis)-2,Dimension=axis,Side="trailing");
    Qbar_i = resize(Qbar,size(Qbar,axis)-2,Dimension=axis,Side="both");

    % Find the left and right slopes of Qi
    S_L = (Qbar_i-Qbar_m1);
    S_R = (Qbar_p1-Qbar_i);

    % Limit the slopes
    S_i = sign(S_L).*min(abs(S_L),abs(S_R))*0.95;
    S_i((S_L.*S_R)<=0) = 0;

    % Add 1 and N cells
    S_1 = (QBC2-QBC1);
    S_N = (QBCEND-QBCEND_m1);

    % S_1 = zeros(size(QBC2)); %(QBC2-QBC1);
    % S_N = zeros(size(QBC2)); %(QBCEND-QBCEND_m1);

    S_i = cat(axis,S_1,S_i,S_N);
end