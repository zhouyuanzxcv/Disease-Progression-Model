function p = plot_mask(mask, color, create_light, face_alpha)
if nargin < 2
    color = [0,0,1];
end
if nargin < 3
    create_light = 1;
end
if nargin < 4
    face_alpha = 0.8;
end

v = mask;
p = patch(isosurface(v,0.5));
isonormals(v,p)
p.FaceColor = color;
p.FaceAlpha = face_alpha;
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
if create_light
    camlight;
end
lighting gouraud
end
