% Add current folder plus all subfolders to the path.
addpath(genpath(pwd));
%% Load example workspace
if ~exist('senderIndex') % Load workspace for testing if not already loaded
    load('matlabWorkspace_SAFTCall.mat');                                      % On local Linux profile on the Server
end





%% Generation of normals
bool_visualize_geometry         = 1;
bool_subplots_for_arbitrary     = 1;
which_geography                 = 3; % 1 = dodekaeder(12 Faces)  »|«    2 = icosaeder(20 Faces)  »|«    3 = arbitrary amount of vectors
num_normals_arbitrary           = 14;
[~, ~, normalen]                        = generate_normals(which_geography,num_normals_arbitrary,bool_visualize_geometry,bool_subplots_for_arbitrary);
flags.reflec_characteristics.normals    = normalen;         % normalen as input for mexfunction

test_voxel = [150 150 150];
%% Plot aperture  (with wrong geometry file I think)
f = figure;

plot3(receiverPosGeom(1,receiverIndex),receiverPosGeom(2,receiverIndex),receiverPosGeom(3,receiverIndex),'xg')
hold on;
plot3(senderPosGeom(1,senderIndex),senderPosGeom(2,senderIndex),senderPosGeom(3,senderIndex),'or')
hold on;
sp = plot3(flags.imageInfos.imageStartpoint(1) +test_voxel(1)*flags.imageInfos.IMAGE_RESOLUTION , flags.imageInfos.imageStartpoint(2) +test_voxel(2)*flags.imageInfos.IMAGE_RESOLUTION,flags.imageInfos.imageStartpoint(3) +test_voxel(3)*flags.imageInfos.IMAGE_RESOLUTION,'dm');
set(sp, 'MarkerFaceColor', get(sp,'Color'));
hold on;




%% Plot cones

% Create Cone for decision area
theta = 26;
h = 0.2;

factor = 0.23;
for i = 1:size(normalen,2)
    d(1,:) =  [0, 0, 0];
    d(2,:) =  [normalen(1,i),normalen(2,i),normalen(3,i)];
    
    [X3,Y3,Z3,~,~,~] = cone(theta,d,h);
    hold on;
    sf = surf(X3+flags.imageInfos.imageStartpoint(1)+test_voxel(1)*flags.imageInfos.IMAGE_RESOLUTION ,  Y3+flags.imageInfos.imageStartpoint(2)+test_voxel(2)*flags.imageInfos.IMAGE_RESOLUTION  ,  (Z3+flags.imageInfos.imageStartpoint(3)+test_voxel(3)*flags.imageInfos.IMAGE_RESOLUTION),'HandleVisibility','off' );
    
    RGB = [254 254 254]/256; % Color: Gray-ish
    sf.FaceColor = RGB;
    %sf.EdgeColor = RGB;
    
    %set(sf,'FaceColor','w', 'FaceAlpha',0);
    
    hold on;
    
    q = quiver3(flags.imageInfos.imageStartpoint(1)+test_voxel(1)*flags.imageInfos.IMAGE_RESOLUTION   ,  flags.imageInfos.imageStartpoint(2)+test_voxel(2)*flags.imageInfos.IMAGE_RESOLUTION   ,   flags.imageInfos.imageStartpoint(3)+test_voxel(3)*flags.imageInfos.IMAGE_RESOLUTION   ,factor*normalen(1,i),factor*normalen(2,i), factor*normalen(3,i),'HandleVisibility','off');
    q.MaxHeadSize   = 10;
    q.LineStyle     = '-';
    q.Marker        = '.';
    q.Color         = [0 0 153]/256;
    q.LineWidth     = 1.5;
    
    % Numbers on arrows
    tt = text(flags.imageInfos.imageStartpoint(1)+test_voxel(1)*flags.imageInfos.IMAGE_RESOLUTION + factor*normalen(1,i)     ,       flags.imageInfos.imageStartpoint(2)+test_voxel(2)*flags.imageInfos.IMAGE_RESOLUTION + factor*normalen(2,i)      ,      ( flags.imageInfos.imageStartpoint(3)+test_voxel(3)*flags.imageInfos.IMAGE_RESOLUTION  + factor*normalen(3,i)) , num2str(i-1+1)); % num2str(i-1) because c++ Starts Counting at 0   ||  num2str(i-1+1) because then it looks nice
    tt.FontSize = 14;
    tt.Color         = [0 0 153]/256;
    tt.FontWeight ='bold';
    
    
end


%% Make plot look nice
set(get(gca,'Title'),'String','Cones','FontSize',15);
xlabel('X', 'FontSize', 15);
ylabel('Y', 'FontSize', 15);
zlabel('Z', 'FontSize', 15);
set(gca, 'fontsize', 15);
set(gca,'CameraViewAngleMode','Manual'); % not bouncing during rotation
view(3);
daspect([1 1 1]);
axis equal
camlight;
camlight(-80,-10);
lighting phong; %flat gouraud phong none
material([0.4 0.4 0.01 10 0.8]) % sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
grid on;
f.CurrentAxes.ZDir = 'Reverse';

