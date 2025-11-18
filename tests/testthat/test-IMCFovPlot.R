# tests/testthat/test-IMCFovPlot.R

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install.packages('testthat') to run the tests.")
}

testthat::test_that("IMCFovPlot (centroids, counts, max) returns a ggplot with expected title tokens", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV001")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  # fake mask (TIFF) so IMCFovPlot can read background size
  img <- magick::image_blank(width = 64, height = 64, color = "black")
  magick::image_write(img, path = file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  # 2 cells, 2 labeled (kept) features
  df <- data.frame(
    cell_id     = 1:2,
    centroid_x  = c(10, 30),
    centroid_y  = c(20, 40),
    cell_area   = c(50, 60),
    "Cd3e_141Pr"  = c(2, 5),
    "Ms4a1_151Eu" = c(3, 7),
    check.names  = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # Build object (segmentation disabled; we just need centroids + mask path)
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    make_sparse = TRUE,                 # avoid "coercing to dgCMatrix" message
    add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE,
    verbose = FALSE
  ))
  obj <- objs[["FOV001"]]
  
  # Plot with centroids/continuous composite
  p <- IMCFovPlot(
    obj,
    features   = c("141Pr", "151Eu"),
    layer      = "counts",
    boundaries = "centroids",
    combine    = "max",
    normalize  = "none",
    transform  = "none",
    verbose    = FALSE
  )
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "boundary=centroids")
  testthat::expect_match(p$labels$title, "combine=max")
  testthat::expect_match(p$labels$title, "layer=counts")
})

testthat::test_that("IMCFovPlot falls back to an available layer and emits a message", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV002")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  img <- magick::image_blank(64, 64, "black")
  magick::image_write(img, file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2,
    centroid_x = c(10, 12), centroid_y = c(10, 14),
    "A_141Pr" = c(1, 3), "B_151Eu" = c(2, 5),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV002"]]
  
  # Request a missing layer -> should message and choose 'counts'
  testthat::expect_message(
    IMCFovPlot(
      obj,
      features   = c("141Pr", "151Eu"),
      layer      = "scale.data",  # not present
      boundaries = "centroids",
      combine    = "sum",
      verbose    = TRUE
    ),
    regexp = "Layer 'scale.data' not found; using 'counts'"
  )
})

testthat::test_that("IMCFovPlot errors when an assay or feature is missing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV003")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(32, 32, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2, centroid_x = c(5, 15), centroid_y = c(8, 20),
    "A_141Pr" = c(3, 4), "B_151Eu" = c(6, 1),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV003"]]
  
  # missing assay
  testthat::expect_error(
    IMCFovPlot(obj, features = c("141Pr"), assay = "NOT_THERE"),
    regexp = "Assay not found"
  )
  
  # missing feature
  testthat::expect_error(
    IMCFovPlot(obj, features = c("DOES_NOT_EXIST"), layer = "counts"),
    regexp = "Some features missing from layer 'counts'"
  )
})

testthat::test_that("IMCFovPlot uses all features if 'features' is omitted and messages count", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV004")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(40, 40, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:3,
    centroid_x = c(10, 20, 30), centroid_y = c(5, 10, 15),
    "A_141Pr" = c(1, 0, 2), "B_151Eu" = c(3, 4, 5),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV004"]]
  
  # Call WITHOUT a 'features' argument
  testthat::expect_message(
    IMCFovPlot(
      obj,
      boundaries = "centroids",
      layer      = "counts",
      combine    = "mean",
      verbose    = TRUE
    ),
    regexp = "No 'features' provided; using all \\(n=2\\)"
  )
})

testthat::test_that("IMCFovPlot export_path writes a file and returns ggplot invisibly", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV005")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(50, 50, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2,
    centroid_x = c(2, 8), centroid_y = c(3, 9),
    "A_141Pr" = c(1, 2), "B_151Eu" = c(1, 3),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV005"]]
  
  out <- file.path(tempdir(), "imc_overlay.png")
  p <- IMCFovPlot(
    obj,
    features   = c("141Pr", "151Eu"),
    boundaries = "centroids",
    combine    = "sum",
    layer      = "counts",
    export_path = out,
    dpi        = 72,
    verbose    = TRUE
  )
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(file.exists(out))
})

testthat::test_that("IMCFovPlot 'rgb' requires exactly 3 features", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV006")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(40, 40, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2,
    centroid_x = c(1, 2), centroid_y = c(3, 4),
    "A_141Pr" = c(1, 2), "B_151Eu" = c(2, 3), "C_160Gd" = c(4, 5),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV006"]]
  
  # too few
  testthat::expect_error(
    IMCFovPlot(obj, features = c("141Pr", "151Eu"), combine = "rgb"),
    regexp = "requires exactly 3 features"
  )
  # correct length 3 works
  p <- IMCFovPlot(obj, features = c("141Pr", "151Eu", "160Gd"), combine = "rgb",
                  boundaries = "centroids", layer = "counts", verbose = FALSE)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("IMCFovPlot applies limits on continuous scale (centroids path)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV007")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(40, 40, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2,
    centroid_x = c(0, 1), centroid_y = c(0, 1),
    "A_141Pr" = c(0, 5), "B_151Eu" = c(1, 2),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV007"]]
  
  p <- IMCFovPlot(obj, features = c("141Pr", "151Eu"),
                  boundaries = "centroids", combine = "max",
                  layer = "counts", limits = c(0, 10), verbose = FALSE)
  testthat::expect_s3_class(p, "ggplot")
  # find the first colour scale and check its limits
  scs <- p$scales$scales
  lims <- NULL
  for (sc in scs) {
    if (!is.null(sc$aesthetics) && any(grepl("colour|color", sc$aesthetics))) {
      lims <- sc$limits
      break
    }
  }
  testthat::expect_equal(lims, c(0, 10))
})

testthat::test_that("IMCFovPlot (segmentation path) returns ggplot when segmentation exists", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if_not_installed("terra")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV008")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  
  # integer label map (3 regions) as the mask (segmentation can be built)
  m <- matrix(0L, 48, 48)
  m[ 5:15,  5:15] <- 1L
  m[20:30, 20:30] <- 2L
  m[32:44, 10:20] <- 3L
  terra::writeRaster(terra::rast(m), file.path(cs_leaf, "Cell_mask.tiff"),
                     overwrite = TRUE, datatype = "INT2U")
  
  df <- data.frame(
    cell_id = 1:3,
    centroid_x = c(10, 25, 15),
    centroid_y = c(10, 25, 40),
    "A_141Pr" = c(1, 2, 3),
    "B_151Eu" = c(3, 2, 1),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # Build with segmentation so IMCFovPlot can fetch coords from the object/tool
  obj <- LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = TRUE, verbose = FALSE
  )[["FOV008"]]
  
  p <- IMCFovPlot(
    obj,
    features   = c("141Pr", "151Eu"),
    layer      = "counts",
    boundaries = "segmentation",
    combine    = "sum",
    verbose    = FALSE
  )
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_match(p$labels$title, "boundary=segmentation")
})

testthat::test_that("IMCFovPlot respects spatial_rotate_k in the title suffix", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV009")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(40, 40, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:2,
    centroid_x = c(0, 2), centroid_y = c(0, 0),
    "A_141Pr" = c(1, 2), "B_151Eu" = c(3, 4),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV009"]]
  
  p <- IMCFovPlot(
    obj,
    features   = c("141Pr", "151Eu"),
    boundaries = "centroids",
    combine    = "mean",
    layer      = "counts",
    spatial_rotate_k = 2L,
    verbose = FALSE
  )
  testthat::expect_match(p$labels$title, "rot=2")
})

testthat::test_that("IMCFovPlot 'stack' uses manual color scale for centroids when colors are provided", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("withr")
  testthat::skip_if_not_installed("rlang")
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("magick")
  testthat::skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- file.path(top, "CellSegmentation", "ProjP", "FOV010")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  magick::image_write(magick::image_blank(40, 40, "black"), file.path(cs_leaf, "Cell_mask.tiff"), format = "tiff")
  
  df <- data.frame(
    cell_id = 1:3,
    centroid_x = c(1, 2, 3), centroid_y = c(1, 2, 3),
    "A_141Pr" = c(1, 2, 3), "B_151Eu" = c(3, 2, 1),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  obj <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[["FOV010"]]
  
  cols <- c("141Pr" = "#1b9e77", "151Eu" = "#d95f02")
  p <- IMCFovPlot(
    obj,
    features   = c("141Pr", "151Eu"),
    boundaries = "centroids",
    combine    = "stack",
    layer      = "counts",
    colors     = cols,
    alpha      = 0.8,
    verbose    = FALSE
  )
  testthat::expect_s3_class(p, "ggplot")
  # Check there is at least one manual colour scale added
  found_manual_colour <- FALSE
  for (sc in p$scales$scales) {
    if (!is.null(sc$aesthetics) && any(grepl("colour|color", sc$aesthetics)) &&
        inherits(sc, "ScaleDiscrete")) {
      found_manual_colour <- TRUE
      break
    }
  }
  testthat::expect_true(found_manual_colour)
})

