function output = circle(origin,radius,theta_start,theta_stop)
    theta = theta_start:0.01:theta_stop;
    output(1,:) = radius*cosd(theta)+origin(1);
    output(2,:) = radius*sind(theta)+origin(2);
end