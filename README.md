# pymargay
PYthon based MAnga Resolved GAlaxY analysis

Authors: Himansh Rathore and Kavin Kumar

Manuscript: Himansh Rathore, Kavin Kumar, Preetish K Mishra, Yogesh Wadadekar, Omkar Bait, Star-forming S0 Galaxies in SDSS-MaNGA: fading spirals or rejuvenated S0s?, Monthly Notices of the Royal Astronomical Society, Volume 513, Issue 1, June 2022, Pages 389–404, https://doi.org/10.1093/mnras/stac871

arXiv pre-print: arXiv:2203.14283 [astro-ph.GA]

Description:

The python module 'pymargay.py' can be used to perform various analysis with resolved Pipe3D data products of MaNGA galaxies. These include obtaining resolved maps of various quanities, making resolved main sequence type of plots, generating radial profiles of various quantities and quite a few other functionalities.

Documentation:

The first thing to do is to supply the addresses to catalogs which will be needed:

'address1' -> global properties table, cross of Pipe3D + Salim
'address2' -> pymorph r-band photometric data
'address3' -> Postion angle correction to pymorph measurements
'address4' -> MaNGA target catalog (Wake et al. 2017)

These catalogs should be in an astropy readable format. If necessary, specify the format as a keyword argument in Table.read()

'address5' -> location of Pipe3D datacubes (in fits format). This address should be a function of MaNGA ID

Description of the functions are given in 'pymargay.py' itself, in the form of comments


