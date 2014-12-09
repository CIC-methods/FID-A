function [par] = readprocpar(fpath, ext)

% READPROCPAR reads in parameter used in FID acquisition
%
%   PAR = READPROCPAR(FPATH) reads the parameters in the the
%   procpar file in the FID directory FPATH
%
% see also WRITEPROCPAR, READFID, WRITEFID


% Written by L Martyn Klassen
% Copyright 2003 Robarts Research Institute

if nargin < 1
   error('CFMM:readprocpar','One input argument required.');
end

if nargin < 2
   ext = [];
end

[data,fpath] = readprocpar2(fpath, ext);

% Initialize outputs
par.path     = fpath;
par.dp       = [];
par.seqfil   = [];
par.nf       = [];
par.ni       = [];
par.sw       = [];
par.ns       = [];
par.ne1      = [];
par.ne2      = [];
par.rfspoil  = [];
par.rfphase  = [];
par.nt       = [];
par.np       = [];
par.nv       = [];
par.nv2      = [];
par.nv3      = [];
par.ne       = [];
par.polarity = [];
par.evenecho = [];
par.tr       = [];
par.te       = [];
par.esp      = [];
par.espincr  = [];
par.nprof    = [];
par.tproj    = [];
par.phi      = [];
par.psi      = [];
par.theta    = [];
par.vpsi     = [];
par.vphi     = [];
par.vtheta   = [];
par.array{1} = [];
par.arraydim = [];
par.seqcon   = [];
par.lro      = [];
par.lpe      = [];
par.lpe2     = [];
par.pro      = [];
par.ppe      = [];
par.ppe2     = [];
par.pss      = [];
par.thk      = [];
par.thk2     = [];
par.pos1     = [];
par.pos2     = [];
par.pos3     = [];
par.vox1     = [];
par.vox2     = [];
par.vox3     = [];
par.nD       = [];
par.cntr     = [];
par.gain     = [];
par.shimset  = [];
par.z1       = [];
par.z2       = [];
par.z3       = [];
par.z4       = [];
par.z5       = [];
par.z6       = [];
par.z7       = [];
par.z8       = [];
par.x1       = [];
par.y1       = [];
par.xz       = [];
par.yz       = [];
par.xy       = [];
par.x3       = [];
par.y3       = [];
par.x4       = [];
par.y4       = [];
par.z1c      = [];
par.z2c      = [];
par.z3c      = [];
par.z4c      = [];
par.xz2      = [];
par.yz2      = [];
par.xz2      = [];
par.yz2      = [];
par.zxy      = [];
par.z3x      = [];
par.z3y      = [];
par.zx3      = [];
par.zy3      = [];
par.z4x        = [];
par.z4y        = [];
par.z5x        = [];
par.z5y        = [];
par.x2y2       = [];
par.z2xy       = [];
par.z3xy       = [];
par.z2x3       = [];
par.z2y3       = [];
par.z3x3       = [];
par.z3y3       = [];
par.z4xy       = [];
par.zx2y2      = [];
par.z2x2y2     = [];
par.z3x2y2     = [];
par.z4x2y2     = [];
par.petable    = [];
par.nrcvrs     = [];
par.trise      = [];
par.at         = [];
par.gro        = [];
par.gmax       = [];
par.intlv      = [];
par.rcvrs      = [];
par.celem      = [];
par.arrayelemts = [];
par.contrast   = [];
par.tep        = [];
par.date       = [];
par.ti         = [];
par.gss2       = [];
par.gss        = [];
par.tpwri      = [];
par.tpwr1      = [];
par.tpwr2      = [];
par.orient     = [];
par.tof        = [];
par.resto      = [];
par.grox       = [];
par.groy       = [];
par.fov        = [];
par.res        = [];
par.npix       = [];
par.nseg       = [];
par.nzseg      = [];
par.waveform   = [];
par.SR         = [];
par.gradfrac   = [];
par.sfrq       = [];
par.B0         = [];
par.dtmap      = [];
par.nnav       = [];
par.tnav       = [];
par.fast       = [];
par.bt         = [];
par.nhomo      = [];
par.fpmult     = [];
par.d1         = [];
par.ss         = [];
par.ssc        = [];
par.r1         = [];
par.r2         = [];
par.ps_coils   = [];
par.coil_array = [];
par.nav        = [];
par.fliplist   = [];
par.varflip    = [];
par.nfreq      = [];
par.freq       = [];
par.flip       = [];
par.flip1      = [];
par.flipprep   = [];
par.seg        = [];
par.state      = [];
par.rfdelay    = [];
par.gro        = [];
par.gimp       = [];
par.SR         = [];
par.readaxis   = [];
par.timescale  = [];
par.etl        = [];
par.grof       = [];
par.Po         = [];
par.Psl        = [];
par.console    = [];
par.shimcoils  = [];
par.spiral_gmax      = [];
par.spiral_gamma     = [];
par.spiral_delay     = [];
par.spiral_tep       = [];
par.dtg              = [];
par.nturns           = [];
par.direction        = [];
par.ninterleave      = [];
par.tn               = [];
par.randomseed       = [];
par.interleave_order = [];
par.profile          = [];
par.image            = [];
par.spiral_version   = [];
par.spiral_filter    = [];
par.spiral_alpha     = [];
par.spiral_density   = [];
par.navigator        = [];
par.tfirst           = [];
par.tpe              = [];
par.ky_order         = [];
par.alternate        = [];
par.offlineAverages  = [];
par.threshold        = [];
par.cluster          = [];
par.weightfit        = [];
par.offlineAverage   = [];
par.bipolar          = [];
par.tpwr             = [];
par.tpwrf            = [];
par.dpwr             = [];
par.dpwrf            = [];
par.echoes           = [];
par.combinedAcq      = [];
par.reductionFactor  = [];
par.referenceUnits   = [];
par.referenceLines   = [];

names = fieldnames(par);
params = fieldnames(data);

for i = 1:length(params);
   par.(params{i}) = data.(params{i}).value;
end

% Make the returned values fully backwards compatible.
for i = 1:length(names)
   if isfield(data,names{i})
      field = data.(names{i});
      switch field.basictype
         case 1
            par.(names{i}) = field.value;
         case {0 2}
            par.(names{i}) = field.value{1};
      end
   end
end

if isfield(data, 'rcvrs')
   par.nrcvrs = data.rcvrs.number;
end

if isfield(data, 'array')
   par.array = data.array.loops;
end

if isfield(data, 'shimcoils')
   par.shimcoils = data.shimcoils.value;
end

if isfield(data, 'ps_coils')
   par.ps_coils = data.ps_coils.value;
end

if isfield(data, 'coil_array')
   par.coil_array = data.coil_array.value;
end

par.path = fpath;

return
