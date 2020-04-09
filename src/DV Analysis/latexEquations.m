figure;
langmuirStr = 'rate\cdot \frac{\frac{[Dl]}{K_D}}{1 + \frac{[Dl]}{K_D}}';
model = '@(p, x) p(1).*(( (x)./p(2)).^p(3) ./( 1 + ((x)./p(2)).^p(3) ) ) + p(4)';
texModel = texlabel(langmuirStr);
text(0.5, 0.5, texModel)