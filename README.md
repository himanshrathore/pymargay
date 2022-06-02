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

For the function "get_delta_sigma_sfr_profile()", the following is the general idea:

  Alternate way of studying suppression and elevation of star formation rate at resolved scales (based on Ellison et al. 2018)

  On y axis: delta_sigma_sfr, on x axis: radial_distance

  1) Divide the control sample galaxies into various bins of total stellar mass.
  2) Suppose one such bin is 9.9 to 10.4 log solar masses, and it has say 150 objects.
  3) Pick a main sample galaxy, suppose it has a stellar mass of 10.2. Properties of this main sample galaxy shall be compared with properties of the 150 control sample galaxies in 2)
  4) Now suppose you want to compare sigma_sfr at a radial distance of 0.5Kpc.
  5) You know the sigma_sfr at a radial distance of 0.5Kpc for the main sample galaxy and the 150 control sample galaxies. So define delta_sigma_sfr = sigma_sfr for this main sample galaxy - sigma_sfr of control sample galaxies, you will have 150 such values of delta_sigma_sfr. Take mean and std of these.
  6) Similarly do it for other radial distances, i.e. 1.5Kpc, 2.5Kpc, ...
  7) You will get the delta_sigma_sfr radial profile for this main sample galaxy 

################################################################################################################################################################
