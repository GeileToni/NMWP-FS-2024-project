% author: Mauro Morini  
% last modified: 26.04.24
clc;clear;close all;

%% plot resonator and create movie
load("NMWP_Project_resonator_mesh_data_v1.05.mat")
load("uh_T_E4_v1.05.mat")

F = struct('cdata',[],'colormap',[]);
fig = figure;
fig.Visible = 'off';
for i = 1:length(T)
    trisurf(t,p(:,1),p(:,2),uh_T(:,i), 'EdgeColor', 'none')
    title("uh_T, t = " + T(i))
    colorbar ;
    axis equal ; axis off ; axis tight ; %colormap ('jet');
    view(2);
    xlabel("x")
    ylabel("y")
    zlabel("uh")
    %zlim([-0.3,0.1])
    drawnow;
    F(i) = getframe(fig);
    disp("t = "+ T(i))
end
movie(figure,F,2, 10)

% %% save movie
% % create video writer object
%  % writerObj = VideoWriter("wavefront hitting resonator");
%  writerObj = VideoWriter("spacy wave animation");
%  % set the frame rate to one frame per second
%  set(writerObj,'FrameRate',10);
%  % open the writer
%  open(writerObj);
%  writeVideo(writerObj,F)
%  close(writerObj)
