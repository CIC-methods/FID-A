% fida_depCheck.m
% Jamie Near, Sunnybrook Research Institute 2024.
% 
% USAGE:
% [flist,plist]=fida_depCheck(in);
% 
% DESCRIPTION:
% Returns the full list of file dependencies and matlab product dependencies 
% within a particular FID-A toolbox (or within all FID-A toolboxes).
% 
% INPUTS:
% in	= The name of the FID-A function that you would like to test for 
%         dependencies, OR the name of the toolbox that you would like to 
%         test for dependencies:
%               -"op":      processing toolbox
%               -"rf":      rf toolbox
%               -"sim":     simulation toolbox
%               -"io":      input-output toolbox
%               -"run":     example run scripts
%               -"all":     all toolboxes (default)
%          
% OUTPUTS:
% flist   = List of all file dependencies
% plist   = List of all Matlab product dependencies

function [flist,plist]=fida_depCheck(in);

if nargin<1
        toolbox = 'all'
end

if exist(in)==2
    [flist,plist]=matlab.codetools.requiredFilesAndProducts(in);
else

    switch in
        case 'all'
            tbs={'simulationTools','rfPulseTools','processingTools','inputOutput','exampleRunScripts'};
        case 'sim'
            tbs={'simulationTools'};
        case 'rf'
            tbs={'rfPulseTools'};
        case 'op'
            tbs={'processingTools'};
        case 'io'
            tbs={'inputOutput'};
        case 'run'
            tbs={'exampleRunScripts'};
        otherwise
            error('ERROR:  toolbox not recognized. Options are: ''all'', ''sim'', ''rf'', ''op'', ''io'', and ''run''');
    end

    %Start by finding the FID-A root directory:
    a=which('op_averaging.m');
    fields=strsplit(a,'/'); %Not sure if this will work on a Windows machine.  Might need to adjust later.
    fields=fields(2:end-2);
    tempStr='';
    for n=1:length(fields);
        tempStr=[tempStr '/' fields{n}];
    end
    FIDAroot=tempStr;

    %Now navigate to the FID-A root directory:
    cd(FIDAroot);

    %Now decend into the desired toolboxes and make a list of all of the files
    %that need to be checked for dependencies:
    filelist={};
    k=1;
    for n=1:length(tbs)
        cd(tbs{n});
        a=dir('./*');
        for m=1:length(a);
            if strcmp(a(m).name,'.') || strcmp(a(m).name,'..')
                %do nothing
                disp('Ignoring ''.'' and ''..'' dirs.');
            elseif a(m).isdir
                cd(a(m).name);
                b=dir('./*');
                %Decend into 1st level sub-directory:
                for l=1:length(b)
                    if strcmp(b(l).name,'.') || strcmp(b(l).name,'..')
                        %do nothing
                        disp('Ignoring ''.'' and ''..'' dirs.');
                    elseif b(l).isdir
                        %do nothing
                        disp('Ignoring 2nd level sub-directory');
                    elseif strcmp(b(l).name(end-1:end),'.m')
                        filelist{k}=[FIDAroot '/' tbs{n} '/' a(m).name '/' b(l).name];
                        k=k+1;
                    end
                end
                cd('..');
            elseif strcmp(a(m).name(end-1:end),'.m')
                filelist{k}=[FIDAroot '/' tbs{n} '/' a(m).name];
                k=k+1;
            end
        end
        cd('..');
    end

    [flist,plist]=matlab.codetools.requiredFilesAndProducts(filelist);


end





