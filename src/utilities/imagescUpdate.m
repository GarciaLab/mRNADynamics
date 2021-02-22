function h = imagescUpdate(ax,im, lims, varargin)

cmap = 'parula';

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if ~isempty(lims)
    h = imagesc(ax, im, double(lims));
else
    h = imagesc(ax, im);
end
colormap(ax, cmap);
ax.Visible = 'off';
ax.DataAspectRatio = [1, 1, 1];
ax.PlotBoxAspectRatio = [1, 1, 1];
    
end