% Benedikt Ebener
% Karlsruhe Institute of Technology
% March 2021
% function: generate_normals. Generates a set of directional vectors which are eqally distributed
% 
%inputs :  which_geography:         1 = dodekaeder(12 Faces)  »|«    2 = icosaeder(20 Faces)  »|«    3 = arbitrary amount of normals
%          num_normals_arbitrary (set only for which_geography == 3 ) : amount of normals
%          bool_visualize_geometry: set whether to show plot or nah
%          bool_subplots_for_arbitrary: do you want a single plot for the abitrary case or do you want two seperate plots?
%
%outputs:  flachen: vertice information for plot
%          ecken: used to calculated normals / directional_vectors
%          normalen: directional vectors
% 
% FYI: Generate animation sections can be executed manually via breakpoint or at the end
% to generate nice GIFs of the plots


function [ flaechen, ecken, normalen] = generate_normals(which_geography,num_normals_arbitrary,bool_visualize_geometry,bool_subplots_for_arbitrary)
% which_geography : 1 = dodekaeder(12 Faces)  »|«    2 = icosaeder(20 Faces)  »|«    3 = arbitrary amount of normals
% num_normals_arbitrary (only for which_geography = 3 ) : amount of normals


if (which_geography == 1)
    %% generate DODEKAEDER
    %Input: a = radius
    %Output: Areas(xyz, 5points, 12areas), points(xyz,20points), normals(xyz,12areas)
    
    
    a = 1; % Set radius for Dodekaeder
    
    h = a*(sqrt(5)-1)/2;
    
    ecken=zeros(3,20);
    ecken(:,1)=[     a   a   a  ];
    ecken(:,2)=[     a   a  -a  ];
    ecken(:,3)=[     a  -a   a  ];
    ecken(:,4)=[     a  -a  -a  ];
    ecken(:,5)=[    -a   a   a  ];
    ecken(:,6)=[    -a   a  -a  ];
    ecken(:,7)=[    -a  -a   a  ];
    ecken(:,8)=[    -a  -a  -a  ];
    
    ecken(:, 9)=[ 0     +h     +(h+a)];
    ecken(:,10)=[ 0     +h     -(h+a)];
    ecken(:,11)=[ 0     -h     +(h+a)];
    ecken(:,12)=[ 0     -h     -(h+a)];
    
    ecken(:,13)=[ +h     +(h+a)     0];
    ecken(:,14)=[ +h     -(h+a)     0];
    ecken(:,15)=[ -h     +(h+a)     0];
    ecken(:,16)=[ -h     -(h+a)     0];
    
    ecken(:,17)=[ +(h+a)    0   +h];
    ecken(:,18)=[ +(h+a)    0   -h];
    ecken(:,19)=[ -(h+a)    0   +h];
    ecken(:,20)=[ -(h+a)    0   -h];
    
    flaechen=zeros(3, 5, 12);
    flaechen(1:3, 1:5,  1)=ecken(:,[ 1 13 15  5  9]);
    flaechen(1:3, 1:5,  2)=ecken(:,[17  3 14  4 18]);
    flaechen(1:3, 1:5,  3)=ecken(:,[ 4 12 10  2 18]);
    flaechen(1:3, 1:5,  4)=ecken(:,[ 1 17 18  2 13]);
    flaechen(1:3, 1:5,  5)=ecken(:,[ 9  5 19  7 11]);
    flaechen(1:3, 1:5,  6)=ecken(:,[ 8 20  6 10 12]);
    flaechen(1:3, 1:5,  7)=ecken(:,[ 1  9 11  3 17]);
    flaechen(1:3, 1:5,  8)=ecken(:,[15  6 20 19  5]);
    flaechen(1:3, 1:5,  9)=ecken(:,[ 8 16  7 19 20]);
    flaechen(1:3, 1:5, 10)=ecken(:,[13  2 10  6 15]);
    flaechen(1:3, 1:5, 11)=ecken(:,[ 7 16 14  3 11]);
    flaechen(1:3, 1:5, 12)=ecken(:,[ 8 12  4 14 16]);
    
    normalen=zeros(3,12);
    for i=1:12
        A=flaechen(1:3,1,i);
        B=flaechen(1:3,2,i);
        C=flaechen(1:3,3,i);
        a=A-B;
        b=A-C;
        normalen(1:3,i)=cross(a,b);
        normalen(1:3,i)=normalen(1:3,i)/norm(normalen(1:3,i));
    end
    
    
    
elseif (which_geography == 2)
    %% generate ICOSAEDER
    %Output: Areas(xyz, 3points, 20areas), normals(xyz,20areas)
    
    phi=(1+sqrt(5))/2;
    
    
    X = [
        1/2 1/2 1/2 1/2 phi/2 phi/2 -phi/2 -phi/2 0 0 0 0;
        -1/2 -1/2 -1/2 -1/2 phi/2 phi/2 -phi/2 -phi/2 0 0 0 0;
        0 0 0 0 1/2 1/2 -1/2 -1/2 phi/2 -phi/2 -phi/2 phi/2
        ];
    
    Y = [
        0 0 0 0 1/2 1/2 1/2 1/2 phi/2 phi/2 -phi/2 -phi/2;
        0 0 0 0 -1/2 -1/2 -1/2 -1/2 phi/2 phi/2 -phi/2 -phi/2;
        phi/2 -phi/2 phi/2 -phi/2 0 0 0 0 1/2 1/2 -1/2 -1/2
        ];
    
    Z = [
        phi/2 phi/2 -phi/2 -phi/2 0 0 0 0 1/2 1/2 1/2 1/2;
        phi/2 phi/2 -phi/2 -phi/2 0 0 0 0 -1/2 -1/2 -1/2 -1/2;
        1/2 1/2 -1/2 -1/2 phi/2 -phi/2 -phi/2 phi/2 0 0 0 0
        ];
    
    X2 = [
        -phi/2 -1/2 -1/2 -1/2 phi/2 1/2 1/2 1/2;
        -1/2 -phi/2 -phi/2 -phi/2 1/2 phi/2 phi/2 phi/2;
        0 0 0 0 0 0 0 0
        ];
    
    Y2 = [
        1/2 0 0 0 1/2 0 0 0;
        0 -1/2 -1/2 1/2 0 -1/2 -1/2 1/2;
        phi/2 -phi/2 -phi/2 phi/2 phi/2 -phi/2 -phi/2 phi/2
        ];
    
    Z2 = [
        0 phi/2 -phi/2 -phi/2 0 phi/2 -phi/2 -phi/2;
        phi/2 0 0 0 phi/2 0 0 0;
        1/2 1/2 -1/2 -1/2 1/2 1/2 -1/2 -1/2
        ];
    
    
    flaechen=zeros(3, 3, 20);
    
    flaechen(1,:,1:12) = X;
    flaechen(2,:,1:12) = Y;
    flaechen(3,:,1:12) = Z;
    
    flaechen(1,:,13:20) = X2;
    flaechen(2,:,13:20) = Y2;
    flaechen(3,:,13:20) = Z2;
    
    normalen=zeros(3,20);
    for i=1:20
        A=flaechen(1:3,1,i);
        B=flaechen(1:3,2,i);
        C=flaechen(1:3,3,i);
        a=A-B;
        b=A-C;
        normalen(1:3,i)=cross(a,b);
        normalen(1:3,i)=normalen(1:3,i)/norm(normalen(1:3,i));
    end
    
    for i=[1, 4, 5, 7, 9, 11, 15, 17, 18, 20]
        normalen(:,i)=-normalen(:,i);
    end
    
    
    
elseif (which_geography == 3)
    %% generate arbitrary new normals
    %Output: normals(xyz)
    
    
    if (num_normals_arbitrary < 14)
        num_normals_arbitrary = 14;
    end
    
    
    rng(1); % Set the seed of the RNG to 1 so that the normals are always the same and the results are compareable
    
    % Uniformly distribute 200 charged particles across unit sphere
    [ecken,flaechen,~,~]=ParticleSampleSphere('N',num_normals_arbitrary,'Dtol',1E-5,'Etol',1E-5,'s',1);
    
    % OUTPUT:
    %   - V /ecken      : N-by-3 array of vertex (i.e., particle) coordinates. When
    %                     'asym'=true, -V(i,:) is the antipodal partner of V(i,:).
    %   - Tri /flachen  : M-by-3 list of face-vertex connectivities. When 'asym'=false,
    %                     Tri is triangulation of V. When 'asym'=true, Tri is
    %                     triangulation of 2N particles [V;-V].
    %   - Ue_i          : N-by-1 array of particle (or particle pair) energies.
    %   - Ue            : k-by-1 array of potential energy values, where k-1 is the
    %                     total number of iterations. Ue(1) and Ue(k) correspond to the
    %                     potential energies of the initial and final particle
    %                     configurations, respectively.
    
    
    normalen = ecken';
    
end







%% Plot

if bool_visualize_geometry == 1
    %% Dodecaeder & Icosaeder
    if  (which_geography == 1 || which_geography == 2)
        
        % Calculate centroid for every face
        pre_face_centres = sum(flaechen,2);             % Flaechen has Dimensions: 3xyz x 5 Vertices x 12 faces -->  sum(A,2) is a column vector containing the sum of each row
        face_centres = squeeze(pre_face_centres) / 5;   % returns an array with the same elements as the input array A, but with dimensions of length 1 removed.
        
        % Generate faces vector depending on number of faces of the geography
        faces = zeros(1, size(flaechen,2) ); % Initialize faces depending on number of faces
        for i_faces = 1 : size(flaechen,2)
            faces(i_faces) = i_faces;
        end
        
        figure(76);
        face_index = 1;
        
        for i = 1 : size(flaechen,3)
            
            % Plot faces
            clear vertices; % otherwise error because vertices is tranformed later on
            vertices(1,:) = flaechen(1 , face_index : face_index + ( size(flaechen,2) - 1 )); % ( size(flaechen,2) - 1 ) -> 2nd Dimension of flaechen is the amount of vertices one face has -> we plot the amount of vertices minus 1
            vertices(2,:) = flaechen(2 , face_index : face_index + ( size(flaechen,2) - 1 ));
            vertices(3,:) = flaechen(3 , face_index : face_index + ( size(flaechen,2) - 1 ));
            vertices = vertices';
            
            
            figure(76);
            pp = patch('Faces',faces,'Vertices', vertices,'FaceColor','red');
            pp.FaceAlpha     = '1';
            hold on;
            
            % Plot normals
            factor = 0.75; % Factor to scale normals in Plot
            q = quiver3(face_centres(1,i),face_centres(2,i),face_centres(3,i),factor*normalen(1,i),factor*normalen(2,i),factor*normalen(3,i));
            q.MaxHeadSize   = 10;
            q.LineStyle     = '-';
            q.Marker        = '.';
            q.Color         = [0.7 0.7 0.7];
            q.LineWidth     = 1.5;
            
            % Numbers on arrows
            tt = text(face_centres(1,i)+factor*normalen(1,i) , face_centres(2,i)+factor*normalen(2,i) , face_centres(3,i) + factor*normalen(3,i), num2str(i-1+1) ); % num2str(i-1) because c++ Starts Counting at 0
            tt.FontSize = 14;
            
            
            face_index  = face_index + size(flaechen,2);
        end
        
        %         % Test Vector
        %         Test_vector = 7*[0 -0.342 0.94];
        %         qt = quiver3(0,0,0 , Test_vector(1),Test_vector(2) ,Test_vector(3));
        %         qt.MaxHeadSize   = 7;
        %         qt.LineStyle     = '-';
        %         qt.Marker        = '.';
        %         qt.Color         = [0.7 0.7 0];
        %         qt.LineWidth     = 1.5;
        
        
        % set initial view settings
        view(3);
        daspect([1 1 1]);
        axis equal
        camlight;
        camlight(-80,-10);
        lighting phong; %flat gouraud phong none
        material([0.4 0.4 0.2 10 0.8]) % sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
        set(gca, 'fontsize', 15);
        
        if (which_geography == 1)
            set(get(gca,'Title'),'String','Dodecahedron','FontSize',25)
        elseif (which_geography == 2)
            set(get(gca,'Title'),'String','Icosahedron','FontSize',25)
        end
        xlabel('x', 'FontSize', 15);
        ylabel('y', 'FontSize', 15);
        zlabel('z', 'FontSize', 15);
        set(gca,'CameraViewAngleMode','Manual'); % not bouncing during rotation
        
        
        
        
        % Generate animation
        if 0
            
            fig = figure(76);
            title('')
            
            % Remove TickeLabels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            set(gca,'ztick',[])
            set(gca,'zticklabel',[])
            
            set(gca,'visible','off') % Turn axis off completely
            
            set(gcf,'color','w'); % Set Plot background white
            set(gcf,'position',[44,934,916,873]);
            
            
            set(gca, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', ...
                'CameraPositionMode', 'manual');
            for idx = 1: 1 : 360*3
                
                camorbit((1/3),0)
                %pause(0.11);
                
                drawnow limitrate
                frame = getframe(fig);
                im{idx} = frame2im(frame);
                
                
            end
            
            
            filename = 'icosahedron.gif'; % Specify the output file name
            for idx = 1: 1 : 360*3
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.03);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.03);
                end
            end
            
            
            
            
        end
        
        
        
        
        
        
        % -------------------------------------------------------------
        % Alternative Plot for geometry
        % -------------------------------------------------------------
        
        h = figure();
        for i = 1:size(flaechen,3)
            plot3([flaechen(1,:,i), flaechen(1,1,i)], [flaechen(2,:,i), flaechen(2,1,i)], [flaechen(3,:,i), flaechen(3,1,i)], '-k', 'LineWidth', 0.8); hold on;
            plot3([0, normalen(1,i)], [0, normalen(2,i)], [0, normalen(3,i)], 'b', 'LineWidth', 1.5);
        end
        view(3);
        daspect([1 1 1]);
        set(gca, 'fontsize', 15);
        xlabel('x', 'FontSize', 15);
        ylabel('y', 'FontSize', 15);
        zlabel('z', 'FontSize', 15);
        
        if (which_geography == 1)
            set(get(gca,'Title'),'String','Dodecahedron','FontSize',25)
        elseif (which_geography == 2)
            set(get(gca,'Title'),'String','Icosahedron','FontSize',25)
        end
        
        set(h,'Color','w');
        set(gca,'CameraViewAngleMode','Manual'); % not bouncing during rotation
        
        
        
        
        
        
        %% Arbitrary amount of normals
    elseif which_geography == 3
        
        
        if (bool_subplots_for_arbitrary==1)
            figure(77);
            ax1 = subplot(1,2,1); %subplot(m,n,p) divides the current figure into an m-by-n grid and creates axes in the position specified by p
        else
            figure(75);
        end
        
        
        
        % Visualize mesh based on computed configuration of particles
        fv=struct('faces',flaechen,'vertices',ecken);
        h=patch(fv,'FaceColor','red');
        set(h,'EdgeColor','k','FaceColor','red', 'FaceAlpha',1)
        set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
        %daspect([1 1 1]);
        view(3)
        set(get(gca,'Title'),'String',['Generated points with N = ', num2str( size(normalen,2)) ] ,'FontSize',20)
        set(gca, 'fontsize', 15);
        set(gca,'CameraViewAngleMode','Manual'); % not bouncing during rotation
        xlabel('x', 'FontSize', 20);
        ylabel('y', 'FontSize', 20);
        zlabel('z', 'FontSize', 20);
        grid on;
        axis equal
        hold on
        plot3(normalen(1,:),normalen(2,:),normalen(3,:),'.k','MarkerSize',15)
        %camlight;
        %camlight(-80,-50);
        lighting phong; %flat gouraud phong none
        material([0.4 0.4 0.2 10 0.8]) % sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
        
        factor = 1.1;
        for i = 1 : size(normalen,2)            
            % Numbers on points
            tt = text(0 + factor*normalen(1,i) , 0 + factor*normalen(2,i) , 0 + factor*normalen(3,i), num2str(i-1+1) ); % num2str(i-1) because c++ Starts Counting at 0
            tt.FontSize = 14;
        end
        
                    
                    
        % Generate animation
        if 0
            fig = figure(75);
            title('')
            
            % Remove TickeLabels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            set(gca,'ztick',[])
            set(gca,'zticklabel',[])
            
            set(gca,'visible','off') % Turn axis off completely
            
            set(gcf,'color','w'); % Set Plot background white
            set(gcf,'position',[44,934,916,873]);
            
            
            set(gca, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', ...
                'CameraPositionMode', 'manual');
            for idx = 1: 1 : 360*3
                camorbit((1/3),0)
                %pause(0.11);
                drawnow limitrate
                frame = getframe(fig);
                im{idx} = frame2im(frame);
            end
            
            filename = 'arbitrary_mesh_200_points.gif'; % Specify the output file name
            for idx = 1: 1 : 360*3
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.03);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.03);
                end
            end
        end
        
        
        
        
        
        
        
        if (bool_subplots_for_arbitrary==1)
            figure(77)
            ax2 = subplot(1,2,2);
            hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); % Connect axes and rotation of subplots
        else
            figure(77);
        end
        
        
        %Plot sphere
        RGB = [236 245 255]/256; % Color: Grey
        n = 30;  % Number of faces of sphere. Only optical. No calculation
        [X,Y,Z] = sphere(n);
        %sphe = surf(X,Y,Z, 'FaceAlpha',0 , 'EdgeAlpha',1,'HandleVisibility','off');
        sphe.FaceColor = RGB;
        sphe.EdgeColor = RGB;
        
        % set initial view settings
        view(3);
        %daspect([1 1 1]); %Sets the data aspect ratio for the current axes. For equal data unit lengths in all directions, use [1 1 1].
        axis equal; % or equal
        camlight;
        camlight(-80,-10);
        lighting gouraud; %flat gouraud phong none
        material([0.4 0.4 0.2 10 0.8]) % sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
        set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
        set(gca, 'fontsize', 20);
        set(get(gca,'Title'),'String','Derived Vectors','FontSize',20)
        set(gca,'CameraViewAngleMode','Manual'); % not bouncing during rotation
        xlabel('x', 'FontSize', 20);
        ylabel('y', 'FontSize', 20);
        zlabel('z', 'FontSize', 20);
        grid on;
        hold on;
        
        if 1
            factor_test = 1;
            % Test Vector
            Test_vector = factor_test*[0 0 1];
            qt = quiver3(0,0,0 , Test_vector(1),Test_vector(2) ,Test_vector(3));
            qt.MaxHeadSize   = 3;
            qt.LineStyle     = '-';
            qt.Marker        = '.';
            qt.Color         = [1 0 0];
            qt.LineWidth     = 1.5;
            % Test Text
            ttest = text(Test_vector(1),Test_vector(2),Test_vector(3), 'Test' );
            ttest.FontSize = 14;
            ttest.Color         = [1 0 0];
        end
        
        
        if 0
            % Plot normals
            for i = 1 : size(normalen,2)
                q = quiver3(0,0,0,normalen(1,i),normalen(2,i),normalen(3,i));
                q.MaxHeadSize   = 10;
                q.LineStyle     = '-';
                q.Marker        = '.';
                q.Color         = [0.2 0.2 0.2];
                q.LineWidth     = 1.5;
                
                % Numbers on arrows
                tt = text(0 + normalen(1,i) , 0 + normalen(2,i) , 0 + normalen(3,i), num2str(i-1+1) ); % num2str(i-1) because c++ Starts Counting at 0
                tt.FontSize = 14;
            end
        end
        rotate3d on;
        
        % Generate animation
        if 0
            fig = figure(77);
            title('')
            
            % Remove TickeLabels
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            set(gca,'ztick',[])
            set(gca,'zticklabel',[])
            
            set(gca,'visible','off') % Turn axis off completely
            
            set(gcf,'color','w'); % Set Plot background white
            set(gcf,'position',[44,934,916,873]);
            
            set(gca, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', ...
                'CameraPositionMode', 'manual');
            for idx = 1: 1 : 360*3
                camorbit((1/3),0)
                %pause(0.11);
                drawnow limitrate
                frame = getframe(fig);
                im{idx} = frame2im(frame);
            end
            
            filename = 'arbitrary_normals_200_points.gif'; % Specify the output file name
            for idx = 1: 1 : 360*3
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.03);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.03);
                end
            end
        end
        
        
        
    end
end



end
