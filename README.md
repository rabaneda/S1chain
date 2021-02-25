# S1chain
Processing chain for Sentinel-1 images (GRD, VV-polarization) in order to retrieve wind and wave parameters, and for object detection.

# Prerequisites
Sentinel-1 Toolbox (S1TBX) or SNAP by the European Space Agency (ESA) must be installed on your local computer for preprocessing purposes.
https://step.esa.int/main/
https://github.com/senbox-org

Its "gpt" command must also be available on your path environment. So, when installing the ESA toolbox allow set SNAP as environmnet variable. Otherwise, you can set it manually afterwards. For more info, check:
http://step.esa.int/docs/tutorials/SNAP_CommandLine_Tutorial.pdf

# Python usage
There is only one script to run: "source/options.py". The rest of scripts will be called later on automatically.

There are some compulsory and optional variables. Please open the script for further instructions. Modify these variables before running the script.

# Acknowledgements

This processing chain was developed at the University of Hull (UK) with funding from the DCS4COP project (http://dcs4cop.eu/).

Thanks to the Nansen Center (NERSC) for the OpenWind repository and the Royal Netherlands Meteorological Institute (KNMI) for the CMOD5 and CMOD7 software.

# References

Verhoef, A., M. Portabella, A. Stoffelen and H. Hersbach, CMOD5.n - the CMOD5 GMF for neutral winds Document external project: 2008, SAF/OSI/CDOP/KNMI/TEC/TN/165, EUMETSAT, 2008.

Corcione, V., Grieco, G., Portabella, M., Nunziata, F., Migliaccio, M.; A novel azimuth cutoff implementation to retrieve sea surface wind speed from SAR imagery; IEEE Transactions on geoscience and remote sensing, 57; 2019.

Valenzuela, G.R.; Theories for the interaction of electromagnetic and oceanic waves - A review; J. of Boundary-Layer Meteorology, 13; 1977.

Koch, W.; Directional analysis of SAR images aiming at wind direction; IEEE Transactions on geosciences and remote sensing, 42; 2004.

Rana, F.M., Adamo, M., Pasquariello, G., De Carolis, G., Morelli, S.; LG-Mod: A modified local gradient (LG) method to retrieve SAR sea surface wind directions in marine coastal areas; J. of Sensors; 2016.

Rana, F.M., Adamo, M., Lucas, R., Blonda, P.; Sea sirface wind retrieval in coastal areas by means of Sentinel-1 and numerical weather prediction model data; J. of Remote Sensing of Environmnet, 225; 2019.

Shao, W., Zhang, Z., Li, X., Li, H.; Ocean wave parameters retrieval from Sentinel-1 SAR imagery; J. of Remote Sensing, 8; 2016.

Shao, W., Hu, Y., Yang, J., Nunziata, F., Sun, J., Li, H., Zuo, J.; An empirical algorithm to retrieve significant wave height from Sentinel-1 Synthetic Aperture Radar Imagery Collected under Cyclonic Conditions; J. of Remote Sensing, 10, 2018.

Greidnaus, H., Alvarez, M., Santamaria, C., Thoorens, F.X., Kuorti, N., Argentieri, P.; The SUMO ship detector algorithm for satellite radar images; J. of Remote Sensing, 9; 2017.

Santamaria, C., Alvarez, M., Greidanus, H., Syrris, V., Soille, P., Argentieri, P.; Mass processing of Sentinel-1 images for maritime surveillance; J. of Remote Sensing, 9; 2017.
