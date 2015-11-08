function pull_3d_scatter_dots( pull_length, source_coordinates )
% pull dots of a 3D scatter plot towards the camera position so that they
% are not hidden by another object (e.g. mesh object). 

if nargin < 2
    source_coordinates = [0 0 0];
end
if nargin < 1 
    pull_length = 1;
end

h = gca;
cam_pos = h.CameraPosition;

direction_vector = cam_pos - source_coordinates;

direction_vector = direction_vector / norm(direction_vector) * pull_length;

for i = 1:length(h.Children)
    if strcmp(h.Children(i).Type,'scatter')
        h.Children(i).XData = h.Children(i).XData + direction_vector(1);
        h.Children(i).YData = h.Children(i).YData + direction_vector(2);
        h.Children(i).ZData = h.Children(i).ZData + direction_vector(3);
    end
%     if strcmp(h.Children(i).Type,'text')
%         h.Children(i).Position = h.Children(i).Position + direction_vector(1);
%     end
end


end

