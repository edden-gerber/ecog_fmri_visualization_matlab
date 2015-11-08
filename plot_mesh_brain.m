function [ handle ] = plot_mesh_brain( brain_mesh, ...
                                       view_position, ...
                                       vertex_color_values, ...
                                       transparency, ...
                                       color_map )
%PLOT_MESH_BRAIN plots a 3D brain image defined by a mesh object
% 
% Usage:
% PLOT_MESH_BRAIN(Brain), when "Brain" is a structure with a "vertices"
%   field and a "faces" field, plots the 3D image defined by this 
%   structure.
% PLOT_MESH_BRAIN(Brain,P) sets the view position, with P being a 
%   [theta phi] two-element vector. Default value is [0 0] (posterior 
%   view). 
% PLOT_MESH_BRAIN(Brain,P,C) where C is a Nx1 vector, and N is the number 
%   of vertices on the brain mesh, defines the surface color values of the
%   brain according to the figure's color map. If C is an Nx3 matrix, then
%   each row defines the RBG color value of the vertex. The default 
%   colormap is based on the upper ("hot") half of the "jet" colormap, with 
%   value 0 mapped to the gray "brain color" so that the minimal value will 
%   be colored gray rather than another color. 
% PLOT_MESH_BRAIN(Brain,P,C,T) with T being a scalar between 0 and 1
%   indicates the level of transparency of the surface mesh, with 0 being
%   completely opaque (default).
% PLOT_MESH_BRAIN(Brain,P,C,T,cm) with cm being a 64x3 matrix, defines the 
%   colormap for vector C. The default color map is "jet", with the lowest 
%   value replaced by a default brain color (light gray). It's also 
%   possible to use the colormap function as for any Matlab figure. 
% h = PLOT_MESH_BRAIN(...) returns the handle of the created patch object. 
% 
% Any of the optional paramters can be skipped by inputting an null array
% ([]) instead. 
% 
% Use the camera toolbar to change viewing angle. Left-clicking on the 
% figure (with no tool selected) will fix the lighting angle (to 
% "headlight"). It's also possible to do this by calling the external 
% fix_lighting function. 
% 
% Written by Edden M. Gerber, Hebrew University of Jerusalem 2015. 
% edden.gerber@gmail.com
%

DEFAULT_BRAIN_COLOR = [0.85 0.85 0.85]; % Light gray

% Handle input
if nargin < 5 || isempty(color_map)
    % Generate "heat" color map with zero being the default brain color
    color_map = jet(64);
    color_map = [interp1(1:2:63,color_map(33:64,1),1:63)' interp1(1:2:63,color_map(33:64,2),1:63)' interp1(1:2:63,color_map(33:64,3),1:63)'];
    % Add as minimum value the default color: 
    color_map = [DEFAULT_BRAIN_COLOR ; color_map];
end
if nargin < 4 || isempty(transparency)
    transparency = 0;
end
if nargin < 3 || isempty(vertex_color_values)
    % no color information - use the default brain color
    vertex_color_values = zeros(length(brain_mesh.vertices),1);
end
if nargin < 2 || isempty(view_position)
    view_position = [0 0];
end
% If input is row vector, change to column
if size(vertex_color_values,1) == 1
    vertex_color_values = vertex_color_values';
end

if ~isfield(brain_mesh,'faces') || ~isfield(brain_mesh,'vertices')
    error('brain_mesh input structure should include fields "faces" and "vertices".');
end

% Plot
handle = trisurf(brain_mesh.faces, brain_mesh.vertices(:,1), brain_mesh.vertices(:,2), brain_mesh.vertices(:,3),...
    'FaceLighting','gouraud');
set(handle,'FaceVertexCData',vertex_color_values);
colormap(color_map);
set(handle,'FaceAlpha',1-transparency);
shading('interp');
material('dull');
axis('xy'); 
axis('tight'); 
axis('equal'); 
axis('off'); 
hold('all'); 
view(view_position);
l = light(); 
camlight(l,'headlight'); % light source is at the camera position
cameratoolbar('Show'); % use this toolbar for 3D navigation
if all(vertex_color_values==0)
    caxis([0 1]); % to force 0 to be at the bottom of the color scale and not the middle. 
end
if isempty(get(gcf,'WindowButtonDownFcn')) % Set only if this function is not already set for this figure
    set(gcf,'WindowButtonDownFcn',@fix_lighting); % Clicking on the figure will fix the light angle
end
if isempty(get(gcf,'Name')) % Set only if the name is not already set for this figure
    set(gcf,'Name','Click on figure (with no tool selected) to fix lighing angle.'); % Clicking on the figure will fix the light angle
end

% Activate 3D rotation tool
rotate3d on;

end

function fix_lighting(fig_handle,~)
% Set the light source of the selected figure (default: current figure) to
% originate from the viewing direction ("headlight"). If there is no light
% source, a new one is created. 

% Get figure handle
if nargin<1
    fig_handle = gcf;
end

% Find the light object
l_handle = findobj(fig_handle,'Type','light');

% If no light source, create one
if isempty(l_handle)
    l_handle = light;
end

% If there is more than one light source, take the last one
if length(l_handle) > 1
    l_handle = l_handle(end);
end

% Set the light at the camera position
camlight(l_handle,'headlight');

end

