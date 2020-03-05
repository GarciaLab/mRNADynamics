function [filteredIm, successFlag]  = filterAttribute(attribute, im)

        filterType = regexp(attribute, '.*(?=_\d)', 'match');
        sigmas = regexp(attribute, '(\d[.]\d)|(\d\d[.]\d)', 'match');
        [filteredIm, successFlag] = filterImage(im, filterType{1}, sigmas);

end