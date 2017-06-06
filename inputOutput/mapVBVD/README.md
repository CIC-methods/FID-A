# mapVBVD

Reads Siemens raw .dat file from VB/VD MRI raw data.

Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de)

Input: filename or simply meas. id, e.g. mapVBVD(122) (if file is in same path)
    Keep empty (mapVBVD() or mapVBVD([])) to open file open dialog.
    Optional arguments: see below

Output: structure of twix_map_obj with elements (if available):
 * .image          image scan
 * .noise          noise scan
 * .phasecor       phase correction scan
 * .phasestab      phase stabilization scan
 * .phasestabRef0  phasestab. ref. (MDH_REFPHASESTABSCAN && !MDH_PHASESTABSCAN)
 * .phasestabRef1  phasestab. ref. (MDH_REFPHASESTABSCAN &&  MDH_PHASESTABSCAN)
 * .refscan        parallel imaging reference scan
 * .refscanPC      phase correction scan for reference data
 * .refscanPS      phase stabilization scan for reference data
 * .refscanPSRef0  phasestab. ref scan for reference data
 * .refscanPSRef1  phasestab. ref scan for reference data
 * .RTfeedback     realtime feedback data
 * .vop            vop rf data (for pTX systems)


## Order of raw data
1. Columns
2. Channels/Coils
3. Lines
4. Partitions
5. Slices
6. Averages
7. (Cardiac-) Phases
8. Contrasts/Echoes
9. Measurements
10. Sets
11. Segments
12. Ida
13. Idb
14. Idc
15. Idd
16. Ide


## Optional parameters/flags
 * removeOS           removes oversampling (factor 2) in read direction
 * doAverage          performs average (resulting avg-dim has thus size 1)
 * ignoreSeg          ignores segment mdh index (works basically the same as
                      the average flag)
 * rampSampRegrid     optional on-the-fly regridding of ramp-sampled readout
 * doRawDataCorrect   enables raw data correction if used in the acquisition
                      (only works for VB atm)

These and more flags can also be set later, e.g "twix.image.flagRemoveOS = 1"


## Examples
```matlab
twix = mapVBVD(measID);

% return all image-data
image_data = twix.image();
% return all image-data with all singular dimensions removed/squeezed:
image_data = twix.image{''}; % '' necessary due to a matlab limitation
% return only data for line numbers 1 and 5; all dims higher than 4 are
% grouped into dim 5):
image_data = twix.image(:,:,[1 5],:,:);
% return only data for coil channels 2 to 6; all dims higher than 4 are
% grouped into dim 5); but work with the squeezed data order
% => use '{}' instead of '()':
image_data = twix.image{:,2:6,:,:,:};
```
So basically it works like regular matlab array slicing (but 'end' is
not supported). Note that there are still a few bugs with more complicated
array slicing (e.g. indices have to be in increasing order; repeated indices not
supported).

### get unsorted raw data (in acq. order)
```matlab
image_data = twix.image.unsorted(); % no slicing supported atm
```

## General remarks

The raw data can be obtained by calling e.g. twix.image() or for
squeezed data twix.image{''} (the '' are needed due to a limitation
of matlab's overloading capabilities).
Slicing is supported as well, e.g. twix_obj.image(:,:,1,:) will return
only the first line of the full data set (all later dimensions are
squeezed into one). Thus, slicing of the "memory-mapped" data objects
works exactly the same as regular matlab array slicing - with one
exception: The keyword 'end' is not supported.

Overloading of the '()' and '{}' operators works by overloading matlab's
built-in 'subsref' function. Matlab calls subsref whenever the operators
'()', '{}', or '.' are called. In the latter case, the overloaded subsref
just calls the built-in subsref function since we don't want to change
the behaviour for '.'-calls. However, this has one minor consequence:
There's no way (afaik) to know whether the original call was terminated
with a semicolon. Thus, a call to e.g. twix.image.NLin will produce
no output with or without semicolon termination. 'a = twix.image.NLin'
will however produce the expected result.
