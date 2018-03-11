function [dist,center] = ...
    find_ica_center2(data,ica_A)
 
 center = nanmean(data);
 m = size(data,1);
 dist = zeros(m,1);
 
 for i = 1:m
     coord = data(i,1:3) - center;
     recov_coord = ica_A*coord';
     dist(i) = norm(recov_coord);
 end
 
%  [sorted,index] = sort(dist,'ascend');
%  dist = [sorted index];
 
%  [~,ind_1d] = sort(l_dist1d(:,2),'ascend');
%  sort_l_dist1d = [l_dist1d(ind_1d,2) l_dist1d(ind_1d,1)];
 
end

 
