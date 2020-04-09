function fim = filterWithImageScience(im)

scale = 8;
nonmaxsup = true;
dims = imagescience.image.Dimensions(512, 512)

Image = imagescience.image.Image(dims)

Edges = imagescience.feature.Edges
fim = Edges.run(im, scale, nonmaxsup);

end

