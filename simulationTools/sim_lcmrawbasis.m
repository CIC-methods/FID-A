%sim_lcmrawbasis.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [RF,out]=sim_lcmrawbasis(n,sw,Bfield,linewidth,metab,tau1,tau2,addref,makeraw,seq)
% 
% DESCRIPTION:
% Generate an LCModel .RAW file to be used as an individual metabolite basis 
% spectrum in an LCModel basis set.  The relevant characteristics of the
% acquisition can be specified (pulse sequence, number of points, spectral
% width, etc)
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% tau1      = first echo time in [s] (if seq='st', tau1 = TE)
% tau2      = second echo time in [s].  (Used in Press, but not used in SE.
%     (If seq='st', tau2=TM).
% addref    = add reference at 0ppm (for use in LCModel makebasis) ['y' or 'n']
% makeraw   = make output file for lcmodel ['y' or 'n']
% seq       = pulse sequence ['se' for Spin Echo or 'p' for Press]
% metab     = one of the following choices
%   'H2O'    = Water
%   'Ala'    = Alanine
%   'Asp'    = Aspartate
%   'PCh'    = PhosphoCholine
%   'Cr'     = Creatine
%   'PCr'    = PhosphoCreatine
%   'GABA'   = Gamma-aminobutyric acid (kaiser)
%   'GABA3'  = Gamma-aminobutyric acid (de Graaf)
%   'Gln'    = Glutamine
%   'Glu'    = Glutamate
%   'GSH'    = Glutathione
%   'Gly'    = Glycine
%   'Ins'    = Myo-inositol
%   'Lac'    = Lactate
%   'NAA'    = N-acetyl aspartate
%   'Scyllo' = Scyllo-inositol
%   'Tau'    = Taurine
%   'Asc'    = Ascorbate (Vitamin C)
%   'bHB'    = beta-Hydroxybutyrate
%   'bHG'    = beta-Hydroxyglutarate
%   'Glc'    = Glucose
%   'NAAG'   = N-acetyl aspartyl glutamate
%   'GPC'    = Glycero-phosphocholine
%   'PE'     = Phosphoryl ethanolamine
%   'Ser'    = Serine

function [RF,out]=sim_lcmrawbasis(n,sw,Bfield,linewidth,metab,tau1,tau2,addref,makeraw,seq)


load('spinSystems.mat');

switch metab
    case 'H2O'
        sys=sysH2O;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        out=op_ampScale(out,2);
        
    case 'Ala'
        sys=sysAla;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
    case 'Asp'
        sys=sysAsp;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'PCh'
        sys1=sysPCh_ch2ch2;
        disp(['simulating metabolite ' metab ' (part 1 of 2) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysPCh_trimethyl;
        disp(['simulating metabolite ' metab ' (part 2 of 2) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        out2=op_ampScale(out2,9);
        out=op_addScans(out1,out2);
    case 'Cr'
        sys=sysCr;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'PCr'
        sys=sysPCr;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'GABA'
        sys=sysGABA;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Gln'
        sys1=sysGln_nh2;
        disp(['simulating metabolite ' metab ' (part 1 of 2) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysGln_ch2ch2ch;
        disp(['simulating metabolite ' metab ' (part 2 of 2) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out=op_addScans(out1,out2);
    case 'Glu'
        sys=sysGlu;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'GSH'
        sys1=sysGSH_cyst;
        disp(['simulating metabolite ' metab ' (part 1 of 3) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysGSH_glut;
        disp(['simulating metabolite ' metab ' (part 2 of 3) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys3=sysGSH_glyc;
        disp(['simulating metabolite ' metab ' (part 3 of 3) with ' num2str(length(sys3.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out3 = sim_spinecho(n,sw,Bfield,linewidth,sys3,tau1);
            case 'p'
                out3 = sim_press(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            case 'st'
                out3 = sim_steam(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out=op_addScans(out1,out2);
        out=op_addScans(out,out3);
        
    case 'Gly'
        sys=sysGly;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Ins'
        sys=sysIns;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Lac'
        sys=sysLac;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'NAA'
        sys=sysNAA;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Scyllo'
        sys=sysScyllo;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out=op_ampScale(out,6);
        
    case 'Tau'
        sys=sysTau;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Asc'
        sys=sysAsc;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'bHB'
        sys=sysbHB;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'bHG'
        sys=sysbHG;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Glc'
        sys1=sysGlc_alpha;
        disp(['simulating metabolite ' metab ' (part 1 of 2) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysGlc_beta;
        disp(['simulating metabolite ' metab ' (part 2 of 2) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out=op_addScans(out1,out2);
        
    case 'NAAG'
        sys1=sysNAAG_acetyl;
        disp(['simulating metabolite ' metab ' (part 1 of 3) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysNAAG_aspartyl;
        disp(['simulating metabolite ' metab ' (part 2 of 3) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys3=sysNAAG_glu;
        disp(['simulating metabolite ' metab ' (part 3 of 3) with ' num2str(length(sys3.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out3 = sim_spinecho(n,sw,Bfield,linewidth,sys3,tau1);
            case 'p'
                out3 = sim_press(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            case 'st'
                out3 = sim_steam(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out=op_addScans(out1,out2);
        out=op_addScans(out,out3);
        
    case 'GPC'
        sys1=sysGPC_gp;
        disp(['simulating metabolite ' metab ' (part 1 of 3) with ' num2str(length(sys1.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out1 = sim_spinecho(n,sw,Bfield,linewidth,sys1,tau1);
            case 'p'
                out1 = sim_press(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            case 'st'
                out1 = sim_steam(n,sw,Bfield,linewidth,sys1,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys2=sysGPC_pCh1;
        disp(['simulating metabolite ' metab ' (part 2 of 3) with ' num2str(length(sys2.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out2 = sim_spinecho(n,sw,Bfield,linewidth,sys2,tau1);
            case 'p'
                out2 = sim_press(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            case 'st'
                out2 = sim_steam(n,sw,Bfield,linewidth,sys2,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        sys3=sysGPC_pCh2;
        disp(['simulating metabolite ' metab ' (part 3 of 3) with ' num2str(length(sys3.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out3 = sim_spinecho(n,sw,Bfield,linewidth,sys3,tau1);
            case 'p'
                out3 = sim_press(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            case 'st'
                out3 = sim_steam(n,sw,Bfield,linewidth,sys3,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
        
        out3=op_ampScale(out3,3);
        out=op_addScans(out1,out2);
        out=op_addScans(out,out3);
        
    case 'PE'
        sys=sysPE;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    case 'Ser'
        sys=sysSer;
        disp(['simulating metabolite ' metab ' with ' num2str(length(sys.shifts)) ' spins...  Please Wait...']);
        switch seq
            case 'se'
                out = sim_spinecho(n,sw,Bfield,linewidth,sys,tau1);
            case 'p'
                out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2);
            case 'st'
                out = sim_steam(n,sw,Bfield,linewidth,sys,tau1,tau2);
            otherwise
                disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
        end
    otherwise
        disp('ERROR:  Unknown metabolite');

end

if addref=='y'||addref=='Y'
    addref=1;
elseif addref=='n'||addref=='N'
    addref=0;
end

if addref
    sysRef.J=0;
    sysRef.shifts=0;
    switch seq
        case 'se'
            ref = sim_spinecho(n,sw,Bfield,linewidth,sysRef,tau1);
        case 'p'
            ref = sim_press(n,sw,Bfield,linewidth,sysRef,tau1,tau2);
        case 'st'
            ref = sim_steam(n,sw,Bfield,linewidth,sysRef,tau1,tau2);
        otherwise
            disp(['ERROR:  Sequence ' seq 'not recognized!!!']);
    end
    out=op_addScans(out,ref);
end

if makeraw=='y'||makeraw=='Y'
    makeraw=1;
elseif makeraw=='n'||makeraw=='N'
    makeraw=0;
end

if makeraw
    RF=io_writelcmraw(out_conj,[metab '.RAW'],metab);
else
    RF=[];
end


   