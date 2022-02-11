%WORK IN PROGRESS

function out = op_CSISelectVoxels(in, x, y)
    arguments
        in (1,1) struct
        x (1,2) double = [-1,-1]
        y (1,2) double = [-1,-1]
    end
    global selected
    selected = [];
    
    if(x(1) < 0 && y(1) < 0)
        fig = sim_plotCSI(in);
        ax = fig.Children(1);
        set(ax, 'ButtonDownFcn', {@get_button_down, in, ax});
        uibutton(fig, 'Position', [440,250,50,20], 'Text', 'Submit', ...
            'ButtonPushedFcn', @(btn,event) closeFigure(btn,fig));
        waitfor(fig);
        x = [min(selected(:,1)), max(selected(:,1))];
        y = [min(selected(:,1)), max(selected(:,1))];
        
    end
        [out, prev_permute, prev_size] = reshapeDimensions(in, [in.dims.x, in.dims.y]);
        out.fids = getData(out)(x, y, :);
        if(isfield(out, 'specs'))
            out.specs = out.specs(x, y, :);
        end
        prev_size(1) = size(getData(out),1);
        prev_size(2) = size(getData(out),2);
        out = permute_back(out, prev_permute, prev_size);
end

function get_button_down(src, ~, in, ax)
    global selected
    mousePos = src.CurrentPoint;

    %[x,y] coordinates
    mousePos=mousePos(1:2:3);

    x = get_closest(in.xCoordinates,mousePos(1));
    y = get_closest(in.yCoordinates, mousePos(2));
    found = false;
    for i = 1:size(selected,1)
        
        bool = selected(i,:) == [x,y];
        if(bool(1) && bool(2))
            found = true;
            selected(i,:) = [];
            tag = sprintf('%d-%d', x, y);
            h = findobj(ax, 'Tag', tag);
            delete(h);
            break
        end
    end
    if(~found)
        selected(end + 1,:) = [x, y];
        pos = [in.xCoordinates(x) - in.deltaX/2, in.yCoordinates(y) - in.deltaY/2, in.deltaX, in.deltaY];
        tag = sprintf('%d-%d', x, y);
        
        rectangle(ax, 'Position', pos, 'Tag', tag);
    end
    
end

function indx = get_closest(vect, target)
    closest = realmax;
    indx = -1;
    for i = 1:length(vect)
        if(closest > abs(vect(i) - target))
            indx = i;
            closest = abs(vect(i) - target);
        end
    end
end

function closeFigure(src, fig)
    close(fig)
end
