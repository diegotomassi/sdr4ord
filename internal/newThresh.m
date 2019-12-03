function [out,outidx] = newThresh(eta,thr)
[p,d] = size(eta);
out = eta;
outidx = ones(p,1);
for j=1:p,
    if norm(eta(j,:))<thr,
        out(j,:) = 0;
        outidx(j) = 0;
    end
end
out = orth(out);