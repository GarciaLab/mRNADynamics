function imagescUpdate(ax,im, lims)
    if ~isempty(lims)
        imagesc(ax, im,lims);
    else
        imagesc(ax, im);
    end
    colormap(ax, 'gray');
    ax.Visible = 'off';
    ax.DataAspectRatio = [1, 1, 1];
    ax.PlotBoxAspectRatio = [1, 1, 1];
end