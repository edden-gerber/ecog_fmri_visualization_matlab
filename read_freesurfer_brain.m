function brain_data = read_freesurfer_brain( fs_subject_dir, fs_dir )
% READ_FREESURFER_BRAIN encapsulates FreeSurfer's shell and Matlab
% functions to read pial and inflated cortical surfaces into Matlab as a
% mesh object (after an anatomical scan has been fully processed by 
% FreeSurfer). Mesh objects are defined by "vertices", e.g. points in 3D
% space, and "faces" which are triangual surfaces connecting three
% vertices. 
% You can use the function "plot_mesh_brain" to plot the resulting cortical
% surfaces as Matlab objects. 
% This code uses FreeSurfer's Matlab function "read_surf". Make sure you
% have the FreeSurfer Matlab library (freesurfer/matlab) on your path. 
% 
% Usage:
% brain_data = READ_FREESURFER_BRAIN(fs_subject_dir) reads the cortical 
%   surface data from the directory given in fs_subject_dir, which should 
%   be the base folder in which the results of FreeSurfer's algorithm for a 
%   subject's brain scan are saved. The output structure brain_data holds 
%   the mesh data (vertices and faces) for the left and right pial and 
%   inflated surfaces. 
% READ_FREESURFER_BRAIN(fs_subject_dir, fs_dir) defines the base working 
%   folder for the FreeSurfer program (e.g. /usr/local/freesurfer). It is
%   best to modify the default value for this variable so that it is not
%   necessary to manually set it. 
% 
% 
% Written by Edden M. Gerber, Hebrew University of Jerusalem 2015, 
% edden.gerber@gmail.com
%

% Handle optional input
DEFAULT_FS_DIR = '/usr/local/freesurfer'; % It's recommended to change this to the FS directory on your machine.

if nargin < 2 || isempty(fs_dir)
    fs_dir = DEFAULT_FS_DIR;
end


% Get transformation matrices. These transformations are necessary to
% load the surfaces from their original FreeSurfer "surface" coordinate
% system to standard XYZ coordinates. 
% This part of the code uses the "system" function to execute FreeSurfer's
% "mri_info" shell command. 
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' fs_dir '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras ' fs_subject_dir '/mri/orig.mgz']);
transformations.ijk2xyz = affine3d(str2num(cmdout)');
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras-tkr ' fs_subject_dir '/mri/orig.mgz']);
transformations.ijk2xyz_FsMesh = affine3d(str2num(cmdout)');


% Create output structure
brain_data = struct;


% Read pial surfaces to Matlab and apply transformations:
% Right hemisphere
[brain_data.pial_right.vertices,brain_data.pial_right.faces] = read_surf(fullfile(fs_subject_dir,'surf','rh.pial'));
brain_data.pial_right.faces = brain_data.pial_right.faces + 1;
brain_data.pial_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_right.vertices);
brain_data.pial_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_right.vertices);
% Left hemisphere
[brain_data.pial_left.vertices,brain_data.pial_left.faces] = read_surf(fullfile(fs_subject_dir,'surf','lh.pial'));
brain_data.pial_left.faces = brain_data.pial_left.faces + 1;
brain_data.pial_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_left.vertices);
brain_data.pial_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_left.vertices);


% Read inflated surfaces to Matlab and apply transformations:
% Right hemisphere
[brain_data.inflated_right.vertices,brain_data.inflated_right.faces] = read_surf(fullfile(fs_subject_dir,'surf','rh.inflated'));
brain_data.inflated_right.faces = brain_data.inflated_right.faces + 1;
brain_data.inflated_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_right.vertices);
brain_data.inflated_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_right.vertices);
% Left hemisphere
[brain_data.inflated_left.vertices,brain_data.inflated_left.faces] = read_surf(fullfile(fs_subject_dir,'surf','lh.inflated'));
brain_data.inflated_left.faces = brain_data.inflated_left.faces + 1;
brain_data.inflated_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_left.vertices);
brain_data.inflated_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_left.vertices);


% Pull apart the two inflated hemispheres so they do not overlap if plotted together (originally they are both placed in the center of the coordinate space:
[brain_data.inflated_left, brain_data.inflated_right] = pull_apart_inflated_hemispheres(brain_data.inflated_left, brain_data.inflated_right, 10);

end

