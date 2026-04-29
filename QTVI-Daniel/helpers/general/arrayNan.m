
function c = arrayNan(array,indexs)
    mask = indexs < length(array);
    c = nan(length(indexs),1);
    c(mask) = array(indexs(mask));
end