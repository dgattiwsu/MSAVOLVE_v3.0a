function [ordered] = reorder_nmsa_2(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
    dist(i,i) = 0;
end

ordered = zeros(nrows,ncols);

% The 1st sequence is set by default
ordered(1,:) = nmsa(1,:);
ind = 1;

for i = 2:nrows
    % zero the entire test column
    dist(:,ind) = 0;
    % find the index of the closest sequence
    [~,newind] = max(dist(ind,:));
    ordered(i,:) = nmsa(newind,:);
    % the next index will be the index of the sequence just selected
    ind = newind;
end

end

