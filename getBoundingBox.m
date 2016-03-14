function newBox = getBoundingBox(bbox, images, H)
    assert(length(images) == length(H));
    n = length(images);
    newBox = bbox;
    for i=1:n
        [w, h, ~] = size(images{i});
        imgBox = H{i}*[1 w 1 w; 1 1 h h; 1 1 1 1];
        imgBox = [imgBox(1, :) ./ imgBox(3, :); ...
                  imgBox(2, :) ./ imgBox(3, :)];
        newBox = [  ceil(min(newBox(1), min(imgBox(1, :)))) ...
                    ceil(max(newBox(2), max(imgBox(1, :)))) ...
                    ceil(min(newBox(3), min(imgBox(2, :)))) ...
                    ceil(max(newBox(4), max(imgBox(2, :))))];
    end
end