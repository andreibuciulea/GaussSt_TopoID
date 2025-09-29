function d = distanz(XCoords,YCoords)
    N = numel(XCoords);
    matrix_X = repmat(XCoords,1,N);
    matrix_Y = repmat(YCoords,1,N);
    d = sqrt((matrix_X-matrix_X').^2 + (matrix_Y-matrix_Y').^2);
end