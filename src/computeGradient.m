function grad = computeGradient(C,cell_Lind,TME_size,ndims)

if ndims == 2
    [xind,yind] = ind2sub(TME_size,cell_Lind);

    grad = [0;0];

    switch yind
        case 1 % cell is at lowest y value and cannot compute gradient below it
            grad(2) = diff(C(xind,yind+[0,1]));
        case TME_size(2) % cell is at highest y value and cannot compute gradient above it
            grad(2) = diff(C(xind,yind+[-1,0]));
        otherwise % use centered difference and cut in half because we need to /(2*h)
            grad(2) = 0.5*diff(C(xind,yind+[-1,1]));
    end

    switch xind
        case 1 % cell is at lowest x value and cannot compute gradient left of it
            grad(1) = diff(C(xind+[0,1],yind));
        case TME_size(1) % cell is at highest x value and cannot compute gradient right of it
            grad(1) = diff(C(xind+[-1,0],yind));
        otherwise % use centered difference and cut in half because we need to /(2*h)
            grad(1) = 0.5*diff(C(xind+[-1,1],yind));
    end

elseif ndims == 3
    [xind,yind,zind] = ind2sub(TME_size,cell_Lind);
    grad = [0;0;0];

    switch zind
        case 1 % cell is at lowest z value and cannot compute gradient below it
            grad(3) = diff(C(xind,yind,zind+[0,1]));
        case TME_size(3) % cell is at highest z value and cannot compute gradient above it
            grad(3) = diff(C(xind,yind,zind+[-1,0]));
        otherwise % use centered difference and cut in half because we need to /(2*h)
            grad(3) = 0.5*diff(C(xind,yind,zind+[-1,1]));
    end

    switch yind
        case 1 % cell is at lowest y value and cannot compute gradient below it
            grad(2) = diff(C(xind,yind+[0,1],zind));
        case TME_size(2) % cell is at highest y value and cannot compute gradient above it
            grad(2) = diff(C(xind,yind+[-1,0],zind));
        otherwise % use centered difference and cut in half because we need to /(2*h)
            grad(2) = 0.5*diff(C(xind,yind+[-1,1],zind));
    end

    switch xind
        case 1 % cell is at lowest x value and cannot compute gradient left of it
            grad(1) = diff(C(xind+[0,1],yind,zind));
        case TME_size(1) % cell is at highest x value and cannot compute gradient right of it
            grad(1) = diff(C(xind+[-1,0],yind,zind));
        otherwise % use centered difference and cut in half because we need to /(2*h)
            grad(1) = 0.5*diff(C(xind+[-1,1],yind,zind));
    end

end
