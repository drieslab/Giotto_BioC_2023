# ------------------------------------------------------#
#                                                       #
#  Giotto Suite: A multi-scale and technology-agnostic  #
#  spatial omics analysis framework                     #
#                                                       #
#  Jiaji George Chen, Joselyn Cristina Chavez-Fuentes,  #
#  Matthew O'Brien, Irzam Sarfraz, Iqra Amin,           #
#  Eddie Ruiz, Pratishtha Guckhool, Guo-Cheng Yuan,     #
#  Ruben Dries                                          #
#                                                       #
#  8/2/2023                                             #
#                                                       #
#                                                       #
# ------------------------------------------------------#

# Tutorial: Giotto object structure and manipulation


# start with basic aggregation
# and workflow diagram
# what do you do when there are datasets with more than one set of polys?
# turns out this particular dataset actually does come with multiple layers of annotations that are different from each other
#
# This can be seen both in the values that are carried in the points and also in the polys provided

# Ensure Giotto Suite is installed. -------------------------- #
i_p = installed.packages()
if(!"Giotto" %in% i_p) devtools::install_github("drieslab/Giotto@suite")
# Ensure Giotto Data is installed
if(!"GiottoData" %in% i_p) devtools::install_github("drieslab/GiottoData")


library(data.table)
library(Giotto)
library(GiottoData)

# Ensure the Python environment for Giotto has been installed
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment
  installGiottoEnvironment()
}
# ------------------------------------------------------------ #



# Create a Giotto object - General method with multiple layers
# The Giotto object holds spatial-omic data and facilitates its analysis and
# visualization.
#
# Data used: a subset of a Vizgen MERSCOPE mouse brain dataset to show flexiblity
# of framework



## provide path to vizgen folder
data_path = system.file('/Mini_datasets/Vizgen/Raw/', package = 'GiottoData')



## 0.1 path to transcripts ####
# --------------------------- #

## each transcript has x, y and z coordinate
tx_path = paste0(data_path, '/', 'vizgen_transcripts.gz')
tx_dt = data.table::fread(tx_path)



## 0.2 path to cell boundaries folder ####
# -------------------------------------- #

## vizgen already provides segmentation information in .hdf5 files
## Here I have already converted the hdf5 files to a simple data.table format

boundary_path = paste0(data_path, '/cell_boundaries/')

z0_polygon_DT = fread(paste0(boundary_path, '/', 'z0_polygons.gz'))
z1_polygon_DT = fread(paste0(boundary_path, '/', 'z1_polygons.gz'))

z0_polygons = createGiottoPolygonsFromDfr(name = 'z0', segmdfr = z0_polygon_DT)
z1_polygons = createGiottoPolygonsFromDfr(name = 'z1', segmdfr = z1_polygon_DT)

plot(z0_polygons)
plot(z1_polygons)


# 1. create dataset with transcript and polygon information ####
# ------------------------------------------------------------------------ #
# Giotto objects should need some form of expression information in order to
# perform analyses
# This can be an expression matrix (with spatial locations)
# or simple feature detections with some polygon annotations
#
# Note the nested inputs that define feature type and spatial unit
viz = createGiottoObject(feat_info = list('rna' = tx_dt[,.(global_x, -global_y, gene, global_z)]),
                         spatial_info = list('z0' = z0_polygons,
                                             'z1' = z1_polygons))

viz # directly returns quick summary of attached info


# add some visualizations

# calculate centroid for each polygon ( = cell)
# this can/will be used when aggregating for example counts to cells
viz = addSpatialCentroidLocations(gobject = viz,
                                  poly_info = paste0('z',0:1),
                                  provenance = list('z0', 'z1'),
                                  return_gobject = TRUE)

showGiottoSpatLocs(viz) # show functions provide a more detailed look at specific slots



#Subcellular data can then be immediately plotted

spatInSituPlotPoints(viz,
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     largeImage_name = 'dapi_z0',
                     point_size = 0.5,
                     plot_method = 'ggplot',
                     show_polygon = TRUE,
                     use_overlap = F,
                     polygon_feat_type = 'z1',
                     polygon_color = 'cyan',
                     polygon_bg_color = 'cyan',
                     polygon_line_size = 0.2,
                     coord_fix_ratio = TRUE,
                     background_color = 'black')



# 3. aggregate information to matrix: polygons and transcripts ####
# --------------------------------------------------------------- #

# we will use the z0 polygon information
# we can set activeSpatUnit() or specify this for each command
# activeSpatUnit(viz) = 'z0' # now you don't need to think about setting spat_unit each time

# note the feat_subset_column and feat_subset_ids
# explanation giotto points and polyons and attached values
tx_z0 = getFeatureInfo(viz)

viz = calculateOverlapRaster(viz,
                             spatial_info = 'z0',
                             feat_info = 'rna',
                             feat_subset_column = 'global_z',
                             feat_subset_ids = 0)

viz = overlapToMatrix(viz,
                      poly_info = 'z0',
                      feat_info = 'rna',
                      name = 'raw')

viz = calculateOverlapRaster(viz,
                             spatial_info = 'z1',
                             feat_info = 'rna',
                             feat_subset_column = 'global_z',
                             feat_subset_ids = 1)

viz = overlapToMatrix(viz,
                      poly_info = 'z1',
                      feat_info = 'rna',
                      name = 'raw')

showGiottoExpression(viz)


# combine the calculated data for z1 and z0 into a new spatial unit called
# 'aggregate'
viz = aggregateStacks(gobject = viz,
                      spat_units = c('z0', 'z1'),
                      feat_type = 'rna',
                      values = 'raw',
                      summarize_expression = 'sum',
                      summarize_locations = 'mean',
                      new_spat_unit = 'aggregate')


# From here normal analysis can continue
activeSpatUnit(viz) = 'aggregate'
viz = filterGiotto(viz, expression_threshold = 1, feat_det_in_min_cells = 1, min_det_feats_per_cell = 1)
viz = normalizeGiotto(viz)
viz = addStatistics(viz)

spatInSituPlotPoints(viz,
                     show_polygon = T,
                     polygon_feat_type = 'aggregate',
                     spat_unit = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill_as_factor = F,
                     polygon_fill = 'nr_feats',
                     polygon_alpha = 1)

viz = createSpatialNetwork(viz, method = 'kNN', k = 8)
viz = spatialAutoCorLocal(viz, method = 'moran', data_to_use = 'expression')

spatPlot2D(viz,
           spat_enr_names = 'expr_moran',
           cell_color = featIDs(viz)[1],
           color_as_factor = FALSE)
# issues with meta-based enrs

activeFeatType(viz) = 'rna'
activeSpatUnit(viz) = 'cell'
## Generate polygons and binning
# Bin values using tessellate.
tx = getFeatureInfo(viz, feat_type = 'rna', return_giottoPoints = TRUE)
e = ext(tx) # get extent (spatial bounding box)
hexbin = tessellate(e, shape = 'hexagon', radius = 10, name = 'hexbin')

viz = setPolygonInfo(viz, hexbin)

pseudo_vis = makePseudoVisium(e)
viz = setPolygonInfo(viz, pseudo_vis, name = 'pseudo_visium')

viz = spatQueryGiottoPolygons(viz,
                              filters = list(pseudo_visium = 'all', hexbin = 'all'))




viz = calculateOverlapRaster(viz,
                             spatial_info = 'hexbin',
                             feat_info = 'rna')

viz = overlapToMatrix(viz,
                      poly_info = 'hexbin',
                      feat_info = 'rna',
                      name = 'raw')
activeSpatUnit(viz) = 'hexbin'
viz = filterGiotto(viz, expression_threshold = 1, feat_det_in_min_cells = 1, min_det_feats_per_cell = 1)
viz = normalizeGiotto(viz)
viz = addStatistics(viz)

viz = addSpatialCentroidLocations(viz,
                                  poly_info = 'hexbin',
                                  init_metadata = FALSE)

spatInSituPlotPoints(viz,
                     polygon_feat_type = 'hexbin',
                     polygon_line_size = 0.2,
                     polygon_fill = 'nr_feats',
                     polygon_fill_as_factor = FALSE,
                     point_size = 0.5,
                     feats = list(rna = featIDs(viz)[1:4]))


# linked spatial queries
# further explanation about the dataset and z stacks
# workflow diagram



## Completed Giotto Objects (short)
## which technologies
## easy to play with
##
viz2 = GiottoData::loadGiottoMini('vizgen')
viz2

# spatial manipulation of certain objects



# autocor

# create pre-gen html with images
# follow this up with pre-loaded objects
















