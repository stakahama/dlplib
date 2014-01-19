dlplib
===

This library includes scripts to be used with the molecular dynamics package, DL\_POLY.

Current inventory:

- `acutil.py`: convert molecular file formats, running antechamber, and reading the resulting *.ac* file 
- `fieldsutil.py`: find angles and dihedrals from atom and bond connectivity information
- `gaffutil.py`: look up GAFF parameters (from *gaff.dat*) according to atom/bond type perceptions determined by antechamber 

See `example.py`. Demo files included.

TODO:

- Incorporate missing parameters suggested by antechamber's `parmchk`. For now, include manually.
- Add function to calculate Lennard Jones interaction parameters.

Note: antechamber and openbabel needs to be installed separately.
