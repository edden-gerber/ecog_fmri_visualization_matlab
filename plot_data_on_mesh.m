function plot_data_on_mesh( coordinates, varargin )
% PLOT_DATA_ON_MESH plots data on a 3D mesh object, e.g. data from
% electrodes/fMRI on a 3D cortical surface. Data is plotted either as a 3D 
% scatter plot or as colored surface patches on the mesh surface itself. 
% The coordinates to be plotted can be given either as XYZ triplets or as
% mesh vertex indexes. The function can only be used when there is a patch 
% object within the current figure. 
% The surface coloring method works by creating a second partially-
% transparent mesh object which overlaps the original. This is done to 
% allow separate coloring schemes for the mesh and the data (e.g. grayscale 
% for cortical curvature and color for overlaid activation map). The 
% original mesh's colormap is frozen so subsequenct changes to the colormap 
% will affect only the plotted data layer. Multiple coloring layers can be 
% created in this way. 
% 
% Note: To plot activation maps on the surface (as for fMRI) rather than 
% discrete points data (as for ECoG), use a patch size value of zero so
% that each vertex will only receive its own color value. This will also
% drastically decrease run time for large maps. Alternatively use the
% function paint_mesh directly (it works the same way but will require less
% input arguments). 
%
% Usage:
% PLOT_DATA_ON_MESH(C) plots circles in the locations specificed
%   by coordinates C. If C is an Nx3 matrix, it is read as XYZ
%   coordinates within the current axis. If it is a vector, it is read as
%   vertex indexes for the mesh object (i.e. the coordinates are those of
%   the indicated mesh vertices). NaN values are ignored. 
% PLOT_DATA_ON_MESH(C,clr) defines that color values of the plotted
%   coordinates. clr can be in any of the following formats:
%   1. Single character (e.g. 'b') - uniform color according to the Matlab 
%   color code scheme. 
%   2. RGB triples - input is an Nx3 matrix where N is the number of
%   plotted coordinates and values range between 0 and 1. 
%   3. Color map - input is an Nx1 numerical vector, and values are mapped 
%   to colors accoring to the figure's colormap. This is the only
%   implemented coloring method for surface patches. 
% PLOT_DATA_ON_MESH(...,'scatter') defines the plotting method as "scatter
%   plot", i.e. indicated coordinates will be plotted using the Matlab's 
%   scatter3 function. This is the default plotting method. 
% PLOT_DATA_ON_MESH(...,'surface') defines the plotting method as "surface
%   patches", i.e. indicated coordinates will be plotted by coloring the 
%   mesh surface in a given radius around them. 
% PLOT_DATA_ON_MESH(...,'markersize',s), where s is either a scalar or an
%   vector of length N where N is the number of plotted coordinates, and 
%   the scatterplot method is used, defines the scatterplot marker size(s).
%   Default value is 100. 
% PLOT_DATA_ON_MESH(...,'scatterarg',arg), where arg is a cell array and 
%   the scatterplot method is used, provides encapsulated optional 
%   arguments to Matlab's scatter3 function to allow full control of its 
%   functionality. For example: 
%   plot_data_on_mesh(...,'scatterarg',{'filled','markertype','*'}). 
% PLOT_DATA_ON_MESH(...,'patchsize',s), where s is a scalar and the
%   "surface" plotting method is used, defines the radius around each
%   plotted coordinate where vertices are colored according to that
%   coordinate's color value. Default value is 3. 
% PLOT_DATA_ON_MESH(...,'overlapmethod',M), when the surface plotting 
%   method is used, defines how overlap between patches is handled. M is a
%   string which can be set to any of the following:
%   1. 'max' - maximal value of the two overlapping coordinates is used, in
%   absolute value but preserving the sign (e.g. max(-3,1)=-3). This is the
%   default value. 
%   2. 'mean' - mean value of the overlapping coordinate is used. 
%   3. 'sum' - summed value of the overlapping coordinate is used. 
% PLOT_DATA_ON_MESH(...,'trasparemcy',T), where the surface plotting method
%   is used and T is a scalar between 0 and 1, defines the degree of
%   transparency of the colored surface patches. Default value is 0. 
% PLOT_DATA_ON_MESH(...,'colormap',cm) where the surface plotting method is
%   used and cm is a 64x3 matrix, defines the surface patches color map. 
%   This can be also set through the colormap function. 
% PLOT_DATA_ON_MESH(...,'shownumbers') will add number labels which will
%   float next to the plotted data (e.g. plotted coordinates will be
%   labeled as 1,2,3...). 
%
% Written by Edden M. Gerber, Hebrew University of Jerusalem 2015,  
% edden.gerber@gmail.com

DEFAULT_PLOT_METHOD = 'scatter';
DEFAULT_MARKER_SIZE = 100;
DEFAULT_MARKER_COLOR = [0 0 1];
DEFAULT_PATCH_SIZE = 3;
DEFAULT_OVERLAP_METHOD = 'max';
DEFAULT_TRANSPARENCY = 0;

% Handle mandatory input
ignore_coord = isnan(coordinates(:,1));
xyz_coord = true;
if size(coordinates,2) == 1 % a vector of mesh vertices instead of 3D coordinates
    h_mesh = findobj(gcf,'Type','Patch');
    if isempty(h_mesh)
        error('No patch object found in the current figure (i.e. no brain plotted).');
    end
    mesh_vert = h_mesh.Vertices;
    vert_idx = coordinates;
    coordinates = mesh_vert(vert_idx(~ignore_coord),:);
    xyz_coord = false;
elseif size(coordinates,2) ~= 3
    error('ERROR: coordinates argument should be a nx3 or a nx1 matrix (3D coordinates or mesh vertex list)');
end

% Handle optional input
method = DEFAULT_PLOT_METHOD; % default
marker_size = DEFAULT_MARKER_SIZE; % default
patch_size = DEFAULT_PATCH_SIZE; % default
surface_overlap_method = DEFAULT_OVERLAP_METHOD;
show_numbers = false;
scatter_color_values = DEFAULT_MARKER_COLOR; % default
surface_color_values = ones(size(coordinates,1),1); % default
transparency_value = DEFAULT_TRANSPARENCY;
color_map = colormap;
scatter_arg = {};
narg = size(varargin,2);
arg = 1;
while arg <= narg
    if ischar(varargin{arg})
        switch lower(varargin{arg})
            case 'scatter'
                method = 'scatter';
                arg = arg + 1;
            case 'surface'
                method = 'surface';
                arg = arg + 1;
            case 'scatterarg'
                if narg > arg && iscell(varargin{arg+1})
                    scatter_arg = varargin{arg+1};
                    arg = arg + 2;
                else
                    error('"ScatterArg" argument should be followed by a cell array of strings');
                end
            case 'patchsize'
                if narg > arg && isscalar(varargin{arg+1})
                    patch_size = varargin{arg+1};
                else
                    error('"PatchSize" argument should be followed by a scalar');
                end
                arg = arg + 2;
            case 'transparency'
                if narg > arg && isscalar(varargin{arg+1})
                    if varargin{arg+1} < 0 || varargin{arg+1} > 1
                        error('"transparency" argument should be followed by a scalar between 0 and 1');
                    else
                        transparency_value = varargin{arg+1};
                    end
                else
                    error('"transparency" argument should be followed by a scalar');
                end
                arg = arg + 2;
            case 'markersize'
                if isnumeric(varargin{arg+1})
                    marker_size = varargin{arg+1};
                else
                    error('"MarkerSize" argument should be numeric');
                end
                arg = arg + 2;
            case 'colormap'
                if isnumeric(varargin{arg+1}) && size(varargin{arg+1},1)==64 && size(varargin{arg+1},2)==3
                    color_map = varargin{arg+1};
                    arg = arg + 2;
                else
                    error('"colormap" argument should be a 64x3 numeric matrix');
                end
            case 'overlapmethod'
                if strcmpi(lower(varargin{arg+1}),'max')
                    surface_overlap_method = 'max';
                    arg = arg + 2;
                elseif strcmpi(lower(varargin{arg+1}),'mean')
                    surface_overlap_method = 'mean';
                    arg = arg + 2;
                elseif strcmpi(lower(varargin{arg+1}),'sum')
                    surface_overlap_method = 'sum';
                    arg = arg + 2;
                else
                    error('"OverlapMethod" argument should be followed by ''max'' or ''mean''.');
                end
            case 'shownumbers'
                show_numbers = true;
                arg = arg + 1;
            otherwise
                if length(varargin{arg}) == 1 % single letter indicates color code
                    scatter_color_values = varargin{arg};
                else
                    error(['Unknown parameter "' varargin{arg} '".']);
                end
                arg = arg + 1;
        end
    elseif islogical(varargin{arg})
        varargin{arg} = double(varargin{arg});
    elseif isnumeric(varargin{arg})
        if size(varargin{arg}(~ignore_coord,:),1) == size(coordinates,1)
            scatter_color_values = varargin{arg}(~ignore_coord,:);
            surface_color_values = varargin{arg}(~ignore_coord,:);
        else
            error('Electrode color data must have the same length as the coordinates');
        end
        arg = arg + 1;
    else
        error('Input argument should be numeric or string.');
    end
end

switch method
    case 'scatter'
        % 3-D scatter plot of electrode coordinates
        scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3),marker_size,scatter_color_values,scatter_arg{:});
        
    case 'surface'
        % First find the mesh data object if we haven't yet
        if ~exist('mesh_vert','var')
            h_mesh = findobj(gcf,'Type','Patch');
            if isempty(h_mesh)
                error('No patch object found in the current figure (i.e. no brain plotted).');
            end
            mesh_vert = h_mesh.Vertices;
        end
        
        % If XYZ coordinates are given, replace the them with the nearest 
        % points on the mesh 
        if xyz_coord
            for e=1:size(coordinates,1)
                dist=[mesh_vert(:,1)-coordinates(e,1) mesh_vert(:,2)-coordinates(e,2) mesh_vert(:,3)-coordinates(e,3)];
                dist=sqrt(sum(dist.^2,2));
                [~,idx] = min(dist);
                coordinates(e,:) = mesh_vert(idx,:);
            end
        end
        
        % Create the color values by masking the electrode color values
        % with a maximal distance mask of each vertex from the electrode
        vertex_color_values = zeros(1,length(mesh_vert));
        if patch_size > 0
            num_overlapping_values = zeros(length(mesh_vert),1);
            for e = 1:size(coordinates,1)
                distx=abs(mesh_vert(:,1)-coordinates(e,1));
                disty=abs(mesh_vert(:,2)-coordinates(e,2));
                distz=abs(mesh_vert(:,3)-coordinates(e,3));
                nearby_vert_list = sqrt(distx.^2+disty.^2+distz.^2) <= patch_size;
            
                num_overlapping_values = num_overlapping_values + nearby_vert_list;
                if strcmp(surface_overlap_method,'sum') % overlapping values are summed
                    vertex_color_values = vertex_color_values + surface_color_values(e) * nearby_vert_list';
                elseif strcmp(surface_overlap_method,'mean') % overlapping values are summed and then averaged (after this loop)
                    vertex_color_values = vertex_color_values + surface_color_values(e) * nearby_vert_list';
                else % max: overlapping values take the max value (in absolute value, but preserving the sign)
                    v = [vertex_color_values ; surface_color_values(e) * nearby_vert_list'];
                    vertex_color_values = sum(v .* [max(abs(v))==abs(v(1,:)) ; max(abs(v))==abs(v(2,:))]);
                end
            end
            if strcmp(surface_overlap_method,'mean')
                num_overlapping_values(num_overlapping_values==0) = 1;
                vertex_color_values = vertex_color_values ./ num_overlapping_values';
            end
        else % patch size is zero - each vertex will only be given its own 
             % color value
            vertex_color_values(vert_idx) = surface_color_values;
        end
        
        % Create a new external semi-transparent surface and paint it:
        transparency_map = ones(length(vertex_color_values),1);
        transparency_map(vertex_color_values ~= 0) = transparency_value;
        
        paint_mesh(vertex_color_values', transparency_map);
        colormap(color_map);
        set(gca,'CLim',[0 max(vertex_color_values)]);
end

% Plot coordinate number labels
if show_numbers
    pushed_out_coord = coordinates * 1.05;
    for e = 1:size(coordinates,1)
        text(pushed_out_coord(e,1),pushed_out_coord(e,2),pushed_out_coord(e,3),num2str(e));
    end
end

end

