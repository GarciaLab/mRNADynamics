function smoothTile(im, sigmaRange, nTiles, varargin)

clrmp = 'hsv';
fun = @(x) x;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

figure('Units', 'normalized', 'Position', [.25, .25, .5, .5]);
tiledlayout('flow', 'TileSpacing', 'none', 'Padding', 'none')

sigmas = linspace(sigmaRange(1),sigmaRange(2), nTiles);
colormap(clrmp);

for s = sigmas
    
    ax = nexttile;
    smth = imgaussfilt(im, s);
    imagesc(ax, fun(smth))
    colorbar;

    ax.Visible = 'off';
    title(ax, ['$\sigma$=', num2str(s)], 'Interpreter', 'latex');
    set(findall(ax, 'type', 'text'), 'Visible', 'on')
    axis image

    drawnow;
    
end

end