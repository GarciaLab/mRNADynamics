function [filteredIm, sucessFlag]  = filterAttribute(attribute, im)

        filterType = regexp(attribute, '.*(?=_\d)', 'match');
        sigmas = regexp(attribute, '(\d[.]\d)|(\d\d[.]\d)', 'match');
        [filteredIm, sucessFlag] = filterImage(im, filterType{1}, sigmas);

end