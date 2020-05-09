function [out] = op_CSIFourierTransform(out, settings)

% Does spectral and spatial fourier transform 
% input:    out         = Twix object of CSI data
%           settings    = Set settings.spectralFT to 1 if you want spectral
%                         Fourier Transform
%                         Set settings.spatialFt to 1 if you want spatial
%                         fourier transform
% output:   out         = Twix object with new out.specs field of the
%                         fourier transformed data

%% Initalization and checks
    tic
    %default to doing fourier transform on untranformed dimensions
    if ~exist('settings', 'var')
        if(out.flags.spectralFT == 0)
            settings.spectralFT = 1;
        end
        if(out.flags.spatialFT == 0)
            settings.spatialFT = 1;
        end
        settings.isCartesian = 1;
    end
    if ~isfield(settings, 'spectralFT')
        if(out.flags.spectralFT == 0)
            settings.spectralFT = 1;
        end
    end
    if ~isfield(settings, 'spatialFt')
        if(out.flags.spatialFT == 0)
            settings.spatialFT = 1;
        end
    end
    if ~isfield(settings, 'isCartesian')
        settings.isCartesian = 0;
    end
    
    %checking if fourier transform is done alread in spectral dimension
    if (settings.spectralFT == 1 && out.flags.spectralFT == 1)
        error('spectralFT already done!')
    end
    %checking if fourier transform is done alread in spatial dimension    
    if (settings.spatialFT == 1 && out.flags.spatialFT == 1)
        error('spatialFT already done!')
    end
    
    %if no specs field is made create a new one
    if ~isfield (out, 'specs')
        out.specs = out.fids;
    end
    
   
%%   fourier transform along the spectral dimension 
    if (settings.spectralFT == 1)
        disp('Calculating spectral dimension');
        
        %fourier transform in the spectral domain
        out.specs = fftshift(fft(out.fids,[],out.dims.t),out.dims.t);
        
        %calculating the frequency
        f=(-out.spectralwidth/2)+(out.spectralwidth/(2*out.sz(1))):...
            out.spectralwidth/(out.sz(1)):...
            (out.spectralwidth/2)-(out.spectralwidth/(2*out.sz(1)));
        
        %calculating the ppm
        ppm=-f/(out.Bo*42.577); 
        ppm=ppm+4.65;
        out.ppm = ppm;
        out.flags.spectralFT = 1;
    end
    
%% spatial dimension fourier transform
    if (settings.spatialFT == 1)
        disp('Calculating spatial dimension');
        
        xCoordinates = -out.fovX/2+out.deltaX/2:out.deltaX:out.fovX/2-out.deltaX/2;
        xCoordinates = xCoordinates - out.imageOrigin(1);
        yCoordinates = -out.fovY/2+out.deltaY/2:out.deltaY:out.fovY/2-out.deltaY/2;
        yCoordinates = yCoordinates - out.imageOrigin(2);

        
        if(settings.isCartesian == 1)
            disp('Applying fast fourier transform');
            
            out.specs = fftshift(fft(out.specs, [], out.dims.x),out.dims.x);
            out.specs = fftshift(fft(out.specs, [], out.dims.y),out.dims.y);
        else

            %calculating kSpace trajectory for fourier transform
            [k_x,k_y] = meshgrid(out.k_XCoordinates, out.k_YCoordinates);
            inTraj = [k_x(:), k_y(:)];

            %calculating cartetian image coordinates for fourier transform
            [spatialCoordinates, xCoordinates, yCoordinates] = op_calculateCartetianTrajectory(out);

            %creating fourier transform matrix for the spatial domain
            sftOperator = sft2_Operator(inTraj, spatialCoordinates, 0, ...
                [numel(out.k_XCoordinates), numel(out.k_YCoordinates)]);

            %finding non spatial dimentions
            dimensionNames = fieldnames(out.dims);
            counter = 1;
            for i = 1:numel(dimensionNames)
                if (out.dims.(dimensionNames{i}) ~= 0)
                    dimentionNumbers(counter) = out.dims.(dimensionNames{i}); %#ok<AGROW>
                    counter = counter + 1;
                end
            end
            dimsToFourierTransform = [out.dims.x, out.dims.y];
            dimentionNumbers = setdiff(dimentionNumbers, dimsToFourierTransform);

            fids = out.specs;
            fidsSize = size(out.specs);
            permuteDim = [dimsToFourierTransform, dimentionNumbers];
            permuteBack = ones([1, numel(permuteDim)]);
            for i = 1:numel(permuteBack)
                for j = 1:numel(permuteDim)
                    if permuteDim(j) == i
                        permuteBack(i) = j;
                    end
                end
            end

            %reshaping fids so it can be multiplied along x and y dimensions 
            fids = permute(fids, permuteDim);
            fids = reshape(fids, [prod(fidsSize(dimsToFourierTransform)), prod(fidsSize(dimentionNumbers))]);

            %initalizing space for fourier transformed fids
            final = ones(size(fids));
            %multiplying fids by fourier transform matrix
            for i = 1:size(fids, 2)    
                final(:,i) = sftOperator*fids(:, i);
            end
            %reshaping and returning the dimensions to normal
            final = reshape(final, [fidsSize(dimsToFourierTransform), fidsSize(dimentionNumbers)]);
            out.specs = permute(final, permuteBack);

            %updating twix parameters
            
        end
    end
    toc
    out.flags.spatialFT = 1;
    %since displaying in the scanner coordinates
    %have to flip the x coordinates
    out.xCoordinates = xCoordinates;
    out.yCoordinates = yCoordinates;
    
end

function sft2_Oper = sft2_Operator(InTraj,OutTraj,Ift_flag)
    %
    % sft2 Perform 2-dimensional slow Fourier transform
    %
    % This function was written by Bernhard Strasser, April 2018.
    %
    %
    % The function performs a slow Fourier transform in two spatial dimensions by using the definition of the Fourier Transform
    % FT(f)(a_j) = sum(k=-N/2;N/2-1) f(x_k)exp(-2pi i x_k a_j)
    % (For even N. For odd, change sum slightly)
    %
    %
    % sft2Operator = sft2(A,Ifft_flag)
    %
    % Input: 
    % -         Ift_flag                    ...     Flag of 0 or 1 determining if an inverse Fourier transform or normal should be performed.
    % -         InputSize                   ...     Size of Input data.
    %
    % Output:
    % -         sft2Operator                ...     Fourier transformation matrix which can be applied to the data by Out = sft2Operator*In
    %
    %
    % Feel free to change/reuse/copy the function. 
    % If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
    % Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
    % File dependancy: ?
    %% 0. Preparations

    % Define Exponent of exp
    if(~Ift_flag)          
        Expy = -2*pi*1i;
    else
        Expy = 2*pi*1i;
    end

    % Define Output Size
    NOut = size(OutTraj,1);

    %% 1. FT

    sft2_Oper = zeros([size(OutTraj,1) size(InTraj,1)]);
    for j=1:NOut
        sft2_Oper(j,:) = exp(Expy*(OutTraj(j,1)*InTraj(:,1) ...
            + OutTraj(j,2)*InTraj(:,2)));
    end

end

