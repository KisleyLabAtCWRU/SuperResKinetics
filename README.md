DOI of 20240823 release: doi: 

# SuperResKinetics
Super-resolution imaging and single molecule kinetics MATLAB analysis

Developed in 2013 by Jixin Chen, Lydia Kisley, Joey Tauzin, and Bo Shuang in Christy Landes' research group (LRG) at Rice University. 
Expanded for 3D HILO image reconstruction by Ricardo Monge Neria in Lydia Kisley's research group at Case Western Reserve University. 

Please see <b>LRG_SuperRes_Kinetics_GuideFINAL.pdf</b> for full description and instructions for use in 2D.
Please see <b>RM_code_summary.doc</b> for description of code for use in 3D and analysis of chromatography stationary phase particles.

The main script for user interfacing is <b>LRG_SuperRes_Main.m</b>

Super-resolution imaging and single molecule kinetics have the powerful potential to understand a variety of interfacial interactions. This analysis takes wide field single molecule data and extracts spatial locations of where molecules are located and how long they are at particular locations. The output includes super-resolution images and temporal statistics for kinetic information. We have achieved ~30 nm spatial  and ~10 ms temporal resolutions based on experimental conditions (higher resolutions can be achieved with improved photophysics/optical equipment). The 3D extension takes the 2D analysis from each z-slice and reconstructs a 3D image. Our analysis has been applied to understand DNA sequencing, chromatography, protein/polymer interactions, biofouling, and corrosion, among others.

If used with 3D reconstruction, please cite:

Monge Neria, R.; Kisley, L. Single-Molecule Imaging in Commercial Stationary Phase Particles Using Highly Inclined and Laminated Optical Sheet Microscopy. Anal. Chem. 2023, 95, 4, 2245–2252. https://doi.org/10.1021/acs.analchem.2c03753

If used soley with 2D, please cite:

Chen, J.; Bremauntz, A.; Kisley, L.; Shuang, B.; Landes, C. F. Super-Resolution mbPAINT for Optical Localization of Single-Stranded DNA. ACS Appl. Mater. Interfaces 2013, 5, 9338-9343. https://doi.org/10.1021/am403984k

Additional references where this method has been used:

Kisley, L.; Chen, J.; Mansur, A. P.; Shuang, B.; Kourentzi, K.; Poongavanam, M. -V.; Chen, W. -H.; Dhamane, S.; Willson, R.C.; Landes, C.F. Unified superresolution experiments and stochastic theory provide mechanistic insight into protein ion-exchange adsorptive separations. <i>Proc. Natl. Acad. Sci. USA.</i> <b>2014</b>, <i>111</i>, 2075-2080. DOI: 10.1073/pnas.1318405111

Kisley, L.; Chen, J.; Mansur, A. P.; Dominguez-Medina, S.; Kulla, E.; Kang, M.; Shuang, B.;  Kourentzi, K.; Poongavanam, M. -V.; Dhamane, S.; Willson, R.C.; Landes, C.F. High ionic strength narrows the population of sites participating in protein ion-exchange adsorption: A single-molecule study. <i>J. Chromatogr. A</i> <b>2014</b>, <i>1343</i>, 135-142. 

Kisley, L.; Patil, U.; Dhamane, S.; Kourentzi, K.; Tauzin, L. J.; Willson, R. C.; Landes, C. F. Competitive multi-component anion exchange adsorption of proteins at the single-molecule level. <i>Analyst</i> <b>2017</b>, <i>142</i>, 3127-3131.

Saini, A.; Messenger, H.; Kisley, L. Fluorophores “turned-on” by corrosion reactions can be detected at the single-molecule level. <i>ACS Appl. Mater. Interfaces</i> <b>2021</b>, <i>13</i>, 2000-2006.
