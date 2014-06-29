function [newpts, scale, t] = normalise(pts)
% scale is vector of scales for each dimension.
% Undo using bsxfun(@plus, bsxfun(@rdivide, newpts, scale), t);

t = mean(pts);
newpts = bsxfun(@minus, pts, t);

scale = ones(1, size(pts, 2));
for i=1:size(pts, 2)
    scale(i) = sqrt(2) / std(newpts(:,i));
    newpts(:,i) = newpts(:,i) * scale(i);
end

end