function [ nmsa ] = ...
    make_ancestor_3( nmsa,branch_start,branch_end,npos,prob_PD )

for i = branch_start:branch_end
    for j=1:npos
    nmsa(i,j) = round(random(prob_PD{j}));
    end
end

end

