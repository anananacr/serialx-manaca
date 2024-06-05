# Serial crystallography at Manacá
Automatic data processing pipeline for Serial Crystallography on Manacá beamline (Sirius, LNLS, Brazil). Scripts communicates with CrystFEL[1] packages and handle snapshot/still-diffraction images.


- Python 3.8.12
- CrystFEL: 0.10.1

Tested on dataset referenced in Crystfel's tutorial (0.10.1) https://gitlab.desy.de/thomas.white/crystfel/-/blob/3619f795/doc/articles/tutorial.rst 

Commands:

- First, edit proc_config.py with your correct parameters

python3 runcrystfel.py -m index_no_cell>out.txt&

python3 runcrystfel.py -m index_all>out.txt&

python3 runcrystfel.py -m merge_all>out.txt&

python3 runcrystfel.py -m fom_plots

python3 runcrystfel.py -m export_mtz


References:
1. T. A. White, R. A. Kirian, A. V. Martin, A. Aquila, K. Nass, A. Barty and H. N. Chapman. "CrystFEL: a software suite for snapshot serial crystallography". J. Appl. Cryst. 45 (2012), p335–341. doi:10.1107/S0021889812002312
