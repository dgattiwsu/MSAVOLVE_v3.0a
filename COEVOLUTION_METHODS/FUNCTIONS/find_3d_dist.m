function [data1,dist1,center1,data2,dist2,center2,data3,dist3,center3,ica_A] = ...
    find_3d_dist(COV,ica_distance_method,edge_scale,plotoption,stability)
% This function lists in descending order the pairs of variables of the
% covariance matrix that are more strongly correlated.

% First we find the 3 principal components of the data.
 [V,D] = eig(COV);
 [m,k] = size(V);
 [l,i] = sort(diag(D),'descend');
 W=V(:,i);
 R=D(i,i);
 
 for n=1:k
     W(:,n)=sign(mean(W(:,n)))*W(:,n);
 end
 
 % First we find the center of the raw distribution and the distance of the
 % components from this center.
 
 data = W(:,1:3);
 for i = 1:size(data,2)
    data1(:,i) = data(:,i) - mean(data(:,i));
 end
 [dist1,center1] = find_center(data1,'SMOOTHED',10e1,edge_scale,plotoption);

 % Then we perfom a PCA analysis on the same data: 

 cov_data = (data1' * data1)/m;
 [V,D] = eig(cov_data);
 k = size(V,2);
 [l,i] = sort(diag(D),'descend');
 W=V(:,i);
 R=D(i,i); 
 for n=1:k
     W(:,n)=sign(mean(W(:,n)))*W(:,n);
 end
 data2 = data1 * W;
 
 % data2 is data1 projected on the principal components. We could do the
 % entire PCA also with the following: [~,data2,~] = princomp(data1);
 % The only difference is that in the second case the sign of the 
 % eigenvectors would be arbitrary.
 
 [dist2,center2] = find_center(data2,'SMOOTHED',10e1,edge_scale,plotoption);
 
 % Finally we perform an Independent Component Analysis of the same data
 % using FastICA. This analysis requires the FastICA package to be on the 
 % path:
 
 switch stability
     case 'ON'
 [ica_S, ica_A, ica_W] = fastica(data','stabilization','on');
     case 'OFF'       
 [ica_S, ica_A, ica_W] = fastica(data');
 end
 
 data3 = ica_S';
 
 % No rational in fixing the signs for reproducibility: this is ICA, not
 % PCA.
 
%  for n=1:3
%      data3(:,n)=sign(mean(data3(:,n)))*data3(:,n);
%  end

 % We calculate center and distances:
 
 switch ica_distance_method
     case 'ORTHO'
 [dist3,center3] = find_center(data3,'SMOOTHED',10e2,edge_scale,plotoption);
     case 'SKEW'
 [dist3,center3] = find_ica_center(data3,ica_A,'SMOOTHED',10e3,edge_scale,plotoption);
 end
 
%  dist3d_1 = NaN(m,m);
%  dist3d_2 = NaN(m,m);
%  for i = 1:m
%      for j = i:m
%          dist3d_1(i,j) = norm(data1(i,1:3) - data1(j,1:3));
%          dist3d_1(j,i) = dist3d_1(i,j);
%          dist3d_2(i,j) = norm(data2(i,1:3) - data2(j,1:3));
%          dist3d_2(j,i) = dist3d_2(i,j);
%      end
%  end
%  
%  l_array = sort_matrix_ascend_2(dist3d_1,1);
%  nonzero_ind = l_array(:,1) ~= 0;
%  l_dist3d_1 = l_array(nonzero_ind,:);
%  l_array = sort_matrix_ascend_2(dist3d_2,1);
%  nonzero_ind = l_array(:,1) ~= 0;
%  l_dist3d_2 = l_array(nonzero_ind,:);
 
end
