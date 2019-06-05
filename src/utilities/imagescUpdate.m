function h = imagescUpdate(ax,im, lims)

    if ~isempty(lims)
        h = imagesc(ax, im,lims);
    else
        h = imagesc(ax, im);
    end
    colormap(ax, 'gray');
    ax.Visible = 'off';
    ax.DataAspectRatio = [1, 1, 1];
    ax.PlotBoxAspectRatio = [1, 1, 1];
    
end