%% Load data
load('demo_data');

%%
%% PLOTTING INTRACRANIAL ELECTRODES AND DATA
%%

%% First we will plot the brain (the right hemisphere):

% You can rotate the brain with the 3D rotation tool (selected by default).
% Since you are also rotating the light source, it needs to be manually 
% reset. To do this, de-select any active tool and click on the figure, or
% type "fix_lighting". 

figure;
plot_mesh_brain(brain_data.pial_right);


%% Now let's add electrodes.

% The data table "electrode_data_table" contains the following information
% for each electrode:
% * Is it visually-responsive (true/false)
% * Response onset latency (in ms)
% * XYZ coordinates
% * The index of the vertex on the mesh brain to which it is closest.

d = electrode_data_table;

% We will plot all the electrodes according to their coordinates. 
% Each electrode will be marked by a black circle. If instead of 'k' we 
% included a numerical vector with length equal to the number of
% electrodes, they would have been colored according to the current color
% map. 

plot_data_on_mesh(d.xyz_coordinates,'k','markersize',30);

% Note that some of the plotted circles may be partly "burried" in the
% brain. To fix this, settle on the view you want and then type:
% pull_3d_scatter_dots(d);
% where d is the distance (e.g., 1). The scatter
% objects will be pulled toward the "camera" position. 


%% Add colored surface patches to show the response latency of each responsive electrode

% For each electrodes, a surface patch is defined by the vertices which are
% within a certain distance of it (default is 3). 
plot_data_on_mesh(d.xyz_coordinates(d.responsive,:), d.onset_latency(d.responsive),'surface');

% change color scale to match response latency range
caxis([0 200]);
colorbar;



%%
%% PLOTTING DATA ON THE BRAIN SURFACE (AS FOR FMRI)
%%

%% For this part of the demo we will plot an inflated version of the same 
%% hemisphere
figure;
% The second parameter is the initial viewing position
plot_mesh_brain(brain_data.inflated_right,[90 0]);


%% But the inflated representation is not so clear - let's plot it again 
%% with gray patches to indicate sulci
% Load the curvaure at each vertex
curv = vertex_data.right.vertex_curvature_index;
% The sulci are the negative values
curv = -sign(curv);

% Plot
figure;
plot_mesh_brain(brain_data.inflated_right,[-50 0], curv); 

% Finally, choose a colormap and color range to best represent the data
caxis([-4 1]);
colormap gray


%% Now we can take a map of fMRI activity for a specific part of the 
%% surface, and plot it over this brain with a separate color map

% The variable "fmri_data" hold a simple artificial map of fMRI activity 
% for a circumscribed region-of-interest on the occipital oortex. It is a 
% vector with length equal to the number of mesh vertices, with every 
% vertex outside the ROI set to NaN. 

% To still be able to see the folding structure under the activity map, we
% will use partial transparency (0 is the default no transparency):
transparency = 0.3;

% Plot
paint_mesh(fmri_data,transparency);

% Finally the colormap and range need to be changed for this "layer"
% because they are still identical to the previous ones. 
caxis([0 40]);
colormap hot

% Note that multiple layers can be painted in this way with different color
% maps. 

%% The same activity map on the non-inflated brain
figure;
plot_mesh_brain(brain_data.pial_right,[-50 0]); 
paint_mesh(fmri_data,0.3);
caxis([0 40]);
colormap hot
