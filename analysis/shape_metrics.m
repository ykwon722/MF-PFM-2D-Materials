function m = shape_metrics(phi, dx, rotdeg)
% Shape descriptors of the largest flake (phi>0). Optional rotation (deg) about
% the image centre before measuring/returning the mask (for overlay tests).
if nargin < 2, dx = 1; end
if nargin < 3, rotdeg = 0; end
bw = double(phi) > 0;
if rotdeg ~= 0, bw = imrotate(bw, rotdeg, 'bilinear', 'crop') > 0.5; end
L = bwlabel(bw); st = regionprops(L,'Area','Perimeter','Solidity','Centroid','MajorAxisLength');
[~,i] = max([st.Area]); s = st(i);
[r,c] = find(L==i); dmax = max(hypot(r - s.Centroid(2), c - s.Centroid(1)));
m.area = s.Area*dx^2; m.solidity = s.Solidity; m.circularity = 4*pi*s.Area/s.Perimeter^2;
m.Rmax_over_Req = dmax/sqrt(s.Area/pi); m.mask = (L==i);
end
