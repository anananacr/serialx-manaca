;Pilatus 2M
photon_energy = 12688
adu_per_eV = 0.0001
clen =  0.12500 
coffset=0 
res = 5814.0  ; 172 micron pixel size

; Define rigid group quadrant for a single panel, asic group and collections for geoptimiser
rigid_group_q0 = 0
rigid_group_a0 = 0
rigid_group_collection_quadrants = q0
rigid_group_collection_asics = a0


; corner_{x,y} set the position of the corner of the detector (in pixels)
; relative to the beam

0/min_fs = 0
0/max_fs = 1474
0/min_ss = 0
0/max_ss = 1678
0/corner_x = -736
0/corner_y = -858
0/fs = x
0/ss = y

bad_beamstop/min_x = -736
bad_beamstop/max_x = 22
bad_beamstop/min_y = -66
bad_beamstop/max_y = 24
