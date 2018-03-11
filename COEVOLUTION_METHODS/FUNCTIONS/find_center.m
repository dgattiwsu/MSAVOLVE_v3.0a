function [dist,center] = find_center(data,spline_method,tol_scale,edge_scale,plotoption)

 x = data(:,1);
 y = data(:,2);
 z = data(:,3);
 
 m = size(data,1);
 
 x1_space = linspace(min(x),max(x),10*m);
 for i = 1:length(x1_space)
     x2_space(i) = sum(x <= x1_space(i));
 end
 y1_space = linspace(min(y),max(y),10*m);
 for i = 1:length(y1_space)
     y2_space(i) = sum(y <= y1_space(i));
 end
 z1_space = linspace(min(z),max(z),10*m);
 for i = 1:length(z1_space)
     z2_space(i) = sum(z <= z1_space(i));
 end

switch spline_method
    case 'LEASTSQUARE'
 [x_spline] = least_square_spline(x1_space,x2_space);
 [y_spline] = least_square_spline(y1_space,y2_space);
 [z_spline] = least_square_spline(z1_space,z2_space);
    case 'SMOOTHED' 
 [x_spline] = smoothed_spline(x1_space,x2_space,tol_scale);
 [y_spline] = smoothed_spline(y1_space,y2_space,tol_scale);
 [z_spline] = smoothed_spline(z1_space,z2_space,tol_scale);
end

 
 % close(gcf)
 
 x_length = length(x1_space);
 x_edge = round(x_length/edge_scale);
 x_zeros = false(x_length,1);
 x_zeros(1:x_edge) = true;
 x_zeros((x_length - x_edge):x_length) = true;
 
 y_length = length(y1_space);
 y_edge = round(y_length/edge_scale);
 y_zeros = false(y_length,1);
 y_zeros(1:y_edge) = true;
 y_zeros((y_length - y_edge):y_length) = true;
 
 z_length = length(z1_space);
 z_edge = round(z_length/edge_scale);
 z_zeros = false(z_length,1);
 z_zeros(1:z_edge) = true;
 z_zeros((z_length - z_edge):z_length) = true;
 
 dx_spline = fnder(x_spline);
 dy_spline = fnder(y_spline);
 dz_spline = fnder(z_spline);

 dx1_space = fnval(dx_spline,x1_space);
 dy1_space = fnval(dy_spline,y1_space);
 dz1_space = fnval(dz_spline,z1_space);
 
 dx1_space(x_zeros) = 0;
 dy1_space(y_zeros) = 0;
 dz1_space(z_zeros) = 0;

 switch plotoption
     case 'PLOT'
 figure;
 plot(x1_space,dx1_space,'b');hold on
 plot(y1_space,dy1_space,'r');
 plot(z1_space,dz1_space,'g');hold off
     case 'NOPLOT'
 end
 
 [~,x_center_ind] = max(dx1_space); 
 [~,y_center_ind] = max(dy1_space); 
 [~,z_center_ind] = max(dz1_space);
 
 x_center = x1_space(x_center_ind);
 y_center = y1_space(y_center_ind);
 z_center = z1_space(z_center_ind);
 
 center = [x_center y_center z_center];
 
 dist = zeros(m,1);
 
 for i = 1:m
     coord = data(i,1:3) - center;
     dist(i) = norm(coord);
 end
 
%  [sorted,index] = sort(dist,'ascend');
%  dist = [sorted index];
 
%  [~,ind_1d] = sort(l_dist1d(:,2),'ascend');
%  sort_l_dist1d = [l_dist1d(ind_1d,2) l_dist1d(ind_1d,1)];
 
end

 
