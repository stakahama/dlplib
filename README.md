dlplib
===

This library includes scripts to be used with the molecular dynamics package, DL\_POLY.

Current inventory:

- `acutil.py`: convert molecular file formats, runn antechamber, and read the resulting *.ac* file 
- `fieldsutil.py`: find angles and dihedrals from atom and bond connectivity information
- `gaffutil.py`: look up GAFF parameters (from *gaff.dat*) according to atom/bond type perceptions determined by antechamber 

See `example.py`. Demo files included. Converting from *.mol* files downloaded from ChemSpider led to errors when read by antechamber; *.sdf* from PubChem seemed to give successful results.

TODO:

- Incorporate missing parameters suggested by antechamber's `parmchk` (output: *.frcmod*). For now, include manually.
- Add function to calculate Lennard Jones interaction parameters.

Note: antechamber and openbabel needs to be installed separately.
