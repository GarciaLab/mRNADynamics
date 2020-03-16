function easySliderDemo(varargin)
%make a figure with a slider that moves continuously very simply. 
%give the script a function handle @(x, a) fun(x, a) where a is going to be
%controlled on the slider. 
%Set the limits with 'xLimits', [minX, maxX]



%inputs
xLimits = [0 4*pi];
nPoints = 400;
fun  = @(t, a) sin(a*t); %a is the control parameter


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end



%set up figure and axes
FigH = figure('position',[360 500 400 400]);
axes('XLim', xLimits, 'units','pixels', ...
    'position',[100 50 200 200], 'NextPlot', 'add');
movegui(FigH, 'center')



%plotting
x     = linspace(xLimits(1), xLimits(2), nPoints);
y     = fun(x, 1);
LineH = plot(x,y);



%interactive gui stuff
TextH = uicontrol('style','text',...
    'position',[170 340 40 15]);
SliderH = uicontrol('style','slider','position',[100 280 200 20],...
    'min', xLimits(1), 'max', xLimits(2));
addlistener(SliderH, 'Value', 'PostSet', @callbackfn);




    function callbackfn(source, eventdata)
        sliderVal          = get(eventdata.AffectedObject, 'Value');
        LineH.YData  = fun(x, sliderVal);
        TextH.String = num2str(sliderVal);
    end

end