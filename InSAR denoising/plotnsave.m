function plotnsave(u,titletext,filename)
addpath export_fig
figure; imagesc(u);  
% title(titletext); 
colormap('jet')
colorbar; 
imfig = gca;
imfig.FontSize = 12;
filename = [filename '.png'];
export_fig(filename,'-m2')
end