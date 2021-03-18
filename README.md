# Matlab Codes for efish NPIX analysis

Here are some codes to extract and analyze data recorded with NPIX and SpikeGLX and to analyze behavioral responses captured with Spike2.
---------------------------------------------------------

**last updates: 03/18/2021**
- changed filter cutoff setting in some codes that handle the EODf filtering in the envelope codes. Filter cutoffs are now related to the respective envelope frequency.


update 03/16/2021 changes: 
 - envelope gain calculations for behavior and neurons have been revised
 - chirpremoval uses a crosshair to threshhold chirps instead of using a fixed threshhold
 - NPIX_PreAnalysis: envelope part now can handle multiple stimulations of the envelope block (i.e., if                        all frequencies have been repeated multiple times)

Instructions:
Donwload the directory and add to Matlab path (including the subfunctions).

More codes will be added... 
