# SEDfitting
Scripts use for fitting astronomical data to galaxy models (composites of simple stellar populations)

The scripts are designed to take 4-dimensional datasets composed from data from many different instruments/telescopes and match all dimensions with sets for 50K models to figure out information about the stellar population age of a galaxy, its star formation rate and make robust measurements of total stellar mass, gas mass and dust mass.


Main scripts:

fitSED_STARS+DUST.pl - SED fitting code for multiwavelength images ranging from UV to far-infrared light. Stellar and dust models must be generated first using generatePEGASE.pl and generateDUST.pl. This code fits the stellar component first and uses the residual flux in the optical through near-infrared bands to fit the dust component along with the far-infrared. All input images are convolved to a common PSF and platescale.

genearteDrainegrid.pl - code to generate Dust models ala Draine et al 2007

generatePEGASE.pl - code to generate stellar models following the PEGASE.2 codebase
