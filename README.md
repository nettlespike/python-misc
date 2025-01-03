# Miscellaneous Projects
Various snippets of code exploring geospatial packages in Python.
<br>

### CalculatingShorelineMovement.py
This code is used in ArcGIS Pro's modelbuilder to determine the average distance between two polylines. It takes into account the positive and negative distances between, however, it is highly inaccurate due to use of Nearest Point function, and not calculating distance via grid lines.

### GCP_Correction.py
This code is used in ArcGIS Pro's modelbuilder to extract the shoreline in raster mask, and convert it to a polyline featureclass.

### shoreline_extraction.py
This code requires license to the PCI Geomatica products and package. It receives a multispectral satellite imagery of a shoreline, and uses a simple NDWI method to create a mask of the water (16S Bit Signed). Then the mask is run through a sieve function, before being combined with a band ratio mask of (B > R) & (B > NIR) through a logical OR function. The accuracy is significantly affected by ice, sediment, shadow, and turbid water. 

