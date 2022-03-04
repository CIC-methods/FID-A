function coordinates = createCoordinates(bounds, step)
    bounds = abs(bounds);
    coordinates = -bounds + step/2:step: bounds - step/2;
    %due to numerical instability of floating point numbers, somtimes series
    %is missing the last point.  
    if(abs(coordinates(end) - (bounds - step/2)) > 1e-10)
        coordinates(end + 1) = bounds;
    end
end