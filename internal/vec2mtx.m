function alpha = vec2mtx(hvec,p,r)
h = reshape(hvec,r,p);
alpha = h';