function [xq,vq] = resampleSpectrum(x,v,xq)
vq = interp1(x(:),v(:),xq(:),'linear');
vq(xq < x(1)) = v(1);
vq(xq > x(end)) = v(end);
vq = vq';
vq = vq(:);
end