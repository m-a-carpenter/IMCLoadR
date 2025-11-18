# tests/testthat/test-LoadIMCSegmented.R

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install.packages('testthat') to run the tests.")
}

# helper: accept both dense and sparse matrices for layer checks
expect_matrix_like <- function(x) {
  expect_true(is.matrix(x) || methods::is(x, "Matrix"),
              info = paste("Expected base matrix or Matrix::Matrix; got",
                           paste(class(x), collapse = ", ")))
}

# helper: create a tiny workspace (CSV only)
.make_workspace <- function(top, project, fov, df) {
  cs_leaf <- file.path(top, "CellSegmentation", project, fov)
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  invisible(cs_leaf)
}

# helper: write an integer label map to Cell_mask.tiff using terra
.write_label_mask <- function(cs_leaf, mat_int) {
  mask_path <- file.path(cs_leaf, "Cell_mask.tiff")
  r <- terra::rast(mat_int)
  terra::writeRaster(r, filename = mask_path, overwrite = TRUE, datatype = "INT2U")
  mask_path
}

test_that("LoadIMCSegmented builds minimal Seurat v5 object with expected structure", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0",
          "Requires SeuratObject >= 5.0.0")
  
  withr::local_tempdir() -> top
  cs_root <- file.path(top, "CellSegmentation")
  dir.create(file.path(cs_root, "ProjA", "FOV001"), recursive = TRUE, showWarnings = FALSE)
  
  # 3 cells; two labeled features kept, NA and zero-var dropped, unlabeled dropped
  df <- data.frame(
    cell_id     = 1:3,
    centroid_x  = c(10, 20, 30),
    centroid_y  = c( 5, 15, 25),
    cell_area   = c(50, 60, 70),
    "Cd3e_141Pr"    = c(1, 2, 3),
    "Ms4a1_151Eu"   = c(0, 5, 0),
    "BadNA_152Sm"   = c(7, NA, 9),  # dropped by NA hygiene
    "ZeroVar_153Eu" = c(0, 0, 0),   # dropped by zero-var hygiene
    "Noise"         = c(100, 200, 300), # unlabeled -> dropped
    check.names = FALSE
  )
  csv_path <- file.path(cs_root, "ProjA", "FOV001", "SingleCellObjects.csv")
  utils::write.csv(df, csv_path, row.names = FALSE)
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    assay    = "IMC",
    make_sparse = FALSE,
    verbose  = FALSE,
    add_asinh_layer = TRUE,
    asinh_cofactor  = 5,
    put_matrix_in_data_layer = FALSE,
    drop_na_features = TRUE,
    drop_zero_var_features = TRUE,
    spatial_rotate_k = 1L,
    build_segmentation_from_mask = FALSE,
    default_boundary = "centroids"
  ))
  
  # list shape & naming
  expect_type(objs, "list")
  expect_equal(length(objs), 1L)
  expect_true("FOV001" %in% names(objs))
  obj <- objs[["FOV001"]]
  expect_s4_class(obj, "Seurat")
  
  # assay & layers
  expect_true("IMC" %in% SeuratObject::Assays(obj))
  ass <- obj[["IMC"]]
  
  counts <- SeuratObject::LayerData(ass, layer = "counts")
  expect_matrix_like(counts)
  expect_setequal(rownames(counts), c("141Pr", "151Eu"))
  
  asinh_layer <- SeuratObject::LayerData(ass, layer = "asinh")
  expect_matrix_like(asinh_layer)
  expect_equal(dim(asinh_layer), dim(counts))
  
  # meta.data
  md <- obj@meta.data
  expect_true(all(c("cell_area_um2", "fov", "project") %in% colnames(md)))
  expect_equal(unique(md$fov), "FOV001")
  expect_equal(unique(md$project), "ProjA")
  expect_equal(nrow(md), ncol(counts))
  expect_identical(rownames(md), colnames(counts))
  
  # spatial reduction & FOV/default boundary
  expect_true("spatial" %in% SeuratObject::Reductions(obj))
  emb <- SeuratObject::Embeddings(obj[["spatial"]])
  expect_equal(ncol(emb), 2L)
  expect_setequal(colnames(emb), c("SP_1", "SP_2"))
  expect_equal(nrow(emb), ncol(counts))
  
  fov_key <- SeuratObject::DefaultFOV(obj)
  expect_true(nchar(fov_key) > 0)
  expect_identical(SeuratObject::DefaultBoundary(obj[[fov_key]]), "centroids")
  
  # tool data (via S4 slot)
  tool <- obj@tools[["LoadIMCSegmented"]]
  expect_type(tool, "list")
  expect_true(all(c(
    "paths", "feature_meta", "scale_factors",
    "spatial_orientation", "transforms", "provenance", "segmentation_coords"
  ) %in% names(tool)))
  
  expect_identical(tool$paths$fov, "FOV001")
  expect_identical(tool$paths$project, "ProjA")
  expect_identical(normalizePath(tool$paths$single_cell_csv_path, winslash = "/", mustWork = TRUE),
                   normalizePath(csv_path, winslash = "/", mustWork = TRUE))
  
  # transform record for asinh
  expect_true("asinh" %in% names(tool$transforms))
  expect_identical(tool$transforms$asinh$source_layer, "counts")
  expect_identical(tool$transforms$asinh$cofactor, 5)
  
  # feature meta mirrors post-hygiene features and mapping
  fm <- tool$feature_meta[["IMC"]]
  expect_s3_class(fm, "data.frame")
  expect_setequal(rownames(fm), c("141Pr", "151Eu"))
  expect_true(all(fm$channel_group %in% c("labeled", "unlabeled")))
  expect_identical(fm["141Pr", "original_channel", drop = TRUE], "Cd3e_141Pr")
  expect_identical(fm["151Eu", "original_channel", drop = TRUE], "Ms4a1_151Eu")
  
  # spatial orientation recorded
  expect_equal(tool$spatial_orientation$rotate_k, 1L)
  expect_true(is.numeric(tool$spatial_orientation$center$cx))
  expect_true(is.numeric(tool$spatial_orientation$center$cy))
})

test_that("put_matrix_in_data_layer duplicates counts into data layer when requested", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- .make_workspace(
    top, "ProjA", "FOV001",
    data.frame(
      cell_id     = 1:2,
      centroid_x  = c(1, 2),
      centroid_y  = c(3, 4),
      cell_area   = c(10, 20),
      "Cd3e_141Pr"  = c(5, 6),
      "Ms4a1_151Eu" = c(7, 8),
      check.names = FALSE
    )
  )
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    make_sparse = FALSE,
    verbose = FALSE,
    add_asinh_layer = FALSE,
    put_matrix_in_data_layer = TRUE,
    build_segmentation_from_mask = FALSE
  ))
  ass <- objs[[1]][["IMC"]]
  counts <- SeuratObject::LayerData(ass, "counts")
  data   <- SeuratObject::LayerData(ass, "data")
  expect_equal(as.matrix(data), as.matrix(counts))
})

test_that("default_boundary = 'segmentation' falls back to centroids when segmentation is not built", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  .make_workspace(
    top, "ProjA", "FOV001",
    data.frame(
      cell_id     = 1:2,
      centroid_x  = c(1, 2),
      centroid_y  = c(3, 4),
      cell_area   = c(10, 20),
      # two non-constant labeled features to avoid 1-feature Assay5 edge-case
      "Cd3e_141Pr"  = c(5, 6),
      "Ms4a1_151Eu" = c(1, 3),
      check.names = FALSE
    )
  )
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    verbose = FALSE,
    make_sparse = FALSE,
    add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE,
    default_boundary = "segmentation"   # ask for segmentation, but it's disabled
  ))
  
  obj <- objs[[1]]
  fov_key <- SeuratObject::DefaultFOV(obj)
  expect_identical(SeuratObject::DefaultBoundary(obj[[fov_key]]), "centroids")
})

test_that("unlabeled channels are retained when keep_unlabeled_channels = TRUE", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  .make_workspace(
    top, "ProjA", "FOV001",
    data.frame(
      cell_id     = 1:2,
      centroid_x  = c(1, 2),
      centroid_y  = c(3, 4),
      cell_area   = c(10, 20),
      "Cd3e_141Pr" = c(9, 8),   # labeled -> friendly "141Pr"
      "Noise"      = c(1, 2),   # unlabeled; should be kept
      check.names = FALSE
    )
  )
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = TRUE,
    verbose = FALSE,
    make_sparse = FALSE,
    add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE
  ))
  counts <- SeuratObject::LayerData(objs[[1]][["IMC"]], "counts")
  expect_true(all(c("141Pr", "Noise") %in% rownames(counts)))
  
  fm <- objs[[1]]@tools[["LoadIMCSegmented"]]$feature_meta[["IMC"]]
  expect_identical(fm["Noise", "channel_group", drop = TRUE], "unlabeled")
})

test_that("segmentation polygons are built and selected as default boundary when requested", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- .make_workspace(
    top, "ProjA", "FOV001",
    data.frame(
      cell_id     = 1:3,
      centroid_x  = c(10, 35, 15),
      centroid_y  = c(10, 35, 50),
      cell_area   = c(25, 25, 25),
      "Cd3e_141Pr"  = c(5, 2, 3),
      "Ms4a1_151Eu" = c(0, 1, 0),
      check.names = FALSE
    )
  )
  
  # label map with labels 1,2,3; background 0
  m <- matrix(0L, nrow = 64, ncol = 64)
  m[ 5:15,  5:15] <- 1L
  m[25:35, 25:35] <- 2L
  m[45:60, 10:20] <- 3L
  mask_path <- .write_label_mask(cs_leaf, m)
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    make_sparse = FALSE,
    verbose = FALSE,
    add_asinh_layer = FALSE,
    build_segmentation_from_mask = TRUE,
    seg_simplify_tol = 1,
    spatial_rotate_k = 1L,
    default_boundary = "segmentation"
  ))
  
  obj <- objs[["FOV001"]]
  expect_s4_class(obj, "Seurat")
  
  fov_key <- SeuratObject::DefaultFOV(obj)
  expect_true(nchar(fov_key) > 0)
  expect_identical(SeuratObject::DefaultBoundary(obj[[fov_key]]), "segmentation")
  expect_true(!is.null(obj[[fov_key]][["segmentation"]]))
  
  tool <- obj@tools[["LoadIMCSegmented"]]
  expect_true("segmentation_coords" %in% names(tool))
  expect_true(fov_key %in% names(tool$segmentation_coords))
  seg_df <- tool$segmentation_coords[[fov_key]]
  expect_s3_class(seg_df, "data.frame")
  expect_true(all(c("cell", "x", "y") %in% colnames(seg_df)))
  
  # All cells accounted for
  expect_true(all(colnames(obj) %in% unique(seg_df$cell)))
  expect_equal(length(unique(seg_df$cell)), ncol(obj))
  
  # Mask path recorded
  expect_identical(
    normalizePath(tool$paths$cell_mask_path, winslash = "/", mustWork = TRUE),
    normalizePath(mask_path, winslash = "/", mustWork = TRUE)
  )
  
  # Alignment: segmentation bbox ~= spatial embedding bbox
  emb <- SeuratObject::Embeddings(obj[["spatial"]])
  xr_sp <- range(emb[, "SP_1"], na.rm = TRUE)
  yr_sp <- range(emb[, "SP_2"], na.rm = TRUE)
  xr_seg <- range(seg_df$x, na.rm = TRUE)
  yr_seg <- range(seg_df$y, na.rm = TRUE)
  
  tol <- 5.0  # allow small differences from rasterization/simplification/rotation
  expect_lt(max(abs(xr_sp - xr_seg)), tol)
  expect_lt(max(abs(yr_sp - yr_seg)), tol)
})

test_that("segmentation can coexist with centroids when default_boundary = 'centroids'", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  cs_leaf <- .make_workspace(
    top, "ProjB", "FOV010",
    data.frame(
      cell_id     = 1:2,
      centroid_x  = c(20, 40),
      centroid_y  = c(20, 40),
      cell_area   = c(30, 30),
      "Cd3e_141Pr"  = c(7, 1),
      "Ms4a1_151Eu" = c(4, 2),  # ensure >= 2 kept features
      check.names = FALSE
    )
  )
  
  m <- matrix(0L, nrow = 48, ncol = 48)
  m[ 8:16,  8:16] <- 1L
  m[30:40, 30:40] <- 2L
  .write_label_mask(cs_leaf, m)
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels = FALSE,
    make_sparse = FALSE,
    verbose = FALSE,
    add_asinh_layer = FALSE,
    build_segmentation_from_mask = TRUE,
    spatial_rotate_k = 0L,
    default_boundary = "centroids"
  ))
  
  obj <- objs[["FOV010"]]
  fov_key <- SeuratObject::DefaultFOV(obj)
  
  expect_true(!is.null(obj[[fov_key]][["segmentation"]]))
  expect_identical(SeuratObject::DefaultBoundary(obj[[fov_key]]), "centroids")
  
  seg_df <- obj@tools[["LoadIMCSegmented"]]$segmentation_coords[[fov_key]]
  expect_s3_class(seg_df, "data.frame")
  expect_equal(length(unique(seg_df$cell)), ncol(obj))
})

test_that("meaningful errors are raised for missing inputs", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if_not_installed("rlang")
  
  # Nonexistent top_dir
  expect_error(
    suppressWarnings(LoadIMCSegmented("does/not/exist")),
    regexp = "Missing top_dir"
  )
  
  # Present top_dir but no SingleCellObjects.csv
  withr::local_tempdir() -> top
  dir.create(file.path(top, "CellSegmentation"), recursive = TRUE, showWarnings = FALSE)
  expect_error(
    suppressWarnings(LoadIMCSegmented(
      top_dir = top,
      build_segmentation_from_mask = FALSE,
      verbose = FALSE
    )),
    regexp = "No SingleCellObjects\\.csv files found"
  )
  
  # Missing required meta columns
  #   -> Omit 'cell_id' (keep centroid_x / centroid_y present) to avoid warnings from range()/mean()
  cs_root <- file.path(top, "CellSegmentation", "ProjA", "FOV001")
  dir.create(cs_root, recursive = TRUE, showWarnings = FALSE)
  bad <- data.frame(
    # cell_id intentionally omitted
    centroid_x = c(1, 2),
    centroid_y = c(3, 4),
    "Cd3e_141Pr"  = c(5, 6),
    "Ms4a1_151Eu" = c(1, 2),
    check.names = FALSE
  )
  utils::write.csv(bad, file.path(cs_root, "SingleCellObjects.csv"), row.names = FALSE)
  
  expect_error(
    suppressWarnings(LoadIMCSegmented(
      top_dir = top,
      build_segmentation_from_mask = FALSE,
      verbose = FALSE
    )),
    regexp = "Missing required meta column\\(s\\)"
  )
})

test_that("make_sparse toggles class of the data layer when put_matrix_in_data_layer = TRUE", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  # Small 2x2 dataset
  cs_leaf <- file.path(top, "CellSegmentation", "ProjS", "FOVS")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    cell_id = 1:2, centroid_x = c(0, 2), centroid_y = c(0, 0),
    "A_141Pr" = c(1, 2), "B_151Eu" = c(3, 4), check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # FALSE: expect a base matrix stored in the data layer
  objA <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    put_matrix_in_data_layer = TRUE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  dataA <- SeuratObject::LayerData(objA[["IMC"]], "data")
  expect_true(is.matrix(dataA))           # dense
  
  # TRUE: expect a sparse Matrix stored in the data layer
  skip_if_not_installed("Matrix")
  objB <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    put_matrix_in_data_layer = TRUE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  dataB <- SeuratObject::LayerData(objB[["IMC"]], "data")
  expect_true(methods::is(dataB, "Matrix"))  # sparse
  
  # Same numeric content
  expect_equal(as.matrix(dataA), as.matrix(dataB))
})

test_that("verbose controls emitted messages", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  
  withr::local_tempdir() -> top
  
  # Create both roots so the loader doesn't emit the "missing tiff" message.
  dir.create(file.path(top, "CellSegmentation", "ProjV", "FOVV"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(top, "tiff_multi_page_export", "ProjV"), recursive = TRUE, showWarnings = FALSE)
  
  # Use >= 2 cells and make both features non-constant so they survive hygiene.
  df <- data.frame(
    cell_id    = 1:2,
    centroid_x = c(1, 2),
    centroid_y = c(3, 4),
    cell_area  = c(10, 20),
    "A_141Pr"  = c(5, 6),   # varying
    "B_151Eu"  = c(6, 8),   # varying
    check.names = FALSE
  )
  utils::write.csv(df, file.path(top, "CellSegmentation", "ProjV", "FOVV", "SingleCellObjects.csv"),
                   row.names = FALSE)
  
  # verbose = TRUE -> should emit "Found 1 FOV csv(s)." (we suppress warnings so only messages are checked)
  expect_message(
    suppressWarnings(LoadIMCSegmented(
      top_dir = top,
      keep_unlabeled_channels = FALSE,
      make_sparse = TRUE,                 # avoids "coercing to dgCMatrix" warning
      add_asinh_layer = FALSE,
      build_segmentation_from_mask = FALSE,
      verbose = TRUE
    )),
    regexp = "Found\\s+1\\s+FOV csv\\(s\\)\\."
  )
  
  # verbose = FALSE -> no messages
  expect_no_message(
    suppressWarnings(LoadIMCSegmented(
      top_dir = top,
      keep_unlabeled_channels = FALSE,
      make_sparse = TRUE,
      add_asinh_layer = FALSE,
      build_segmentation_from_mask = FALSE,
      verbose = FALSE
    ))
  )
})


test_that("add_asinh_layer creates correct values and is absent when disabled; cofactor respected", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  # Two features with simple values to check math
  cs_leaf <- file.path(top, "CellSegmentation", "ProjH", "FOVH")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    cell_id = 1:2, centroid_x = c(0, 1), centroid_y = c(0, 1),
    "A_141Pr" = c(0, 10), "B_151Eu" = c(4, 8), check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # Enabled with cofactor = 2
  obj1 <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    add_asinh_layer = TRUE, asinh_cofactor = 2, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  ass1 <- obj1[["IMC"]]
  counts <- SeuratObject::LayerData(ass1, "counts")
  asinh1 <- SeuratObject::LayerData(ass1, "asinh")
  expect_equal(as.matrix(asinh1), asinh(as.matrix(counts) / 2), tolerance = 1e-12)
  # Tool should record cofactor
  expect_identical(obj1@tools[["LoadIMCSegmented"]]$transforms$asinh$cofactor, 2)
  
  # Disabled -> layer absent
  obj2 <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  ass2 <- obj2[["IMC"]]
  expect_false("asinh" %in% SeuratObject::Layers(ass2))
  
  # Different cofactor changes values
  obj3 <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    add_asinh_layer = TRUE, asinh_cofactor = 5, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  asinh3 <- SeuratObject::LayerData(obj3[["IMC"]], "asinh")
  expect_false(isTRUE(all.equal(as.matrix(asinh1), as.matrix(asinh3))))
  expect_equal(as.matrix(asinh3), asinh(as.matrix(counts) / 5), tolerance = 1e-12)
})

test_that("drop_zero_var_features = FALSE keeps zero-variance features", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  cs_leaf <- file.path(top, "CellSegmentation", "ProjZV", "FOVZV")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    cell_id = 1:3, centroid_x = c(0,1,2), centroid_y = c(0,1,2),
    "Const_141Pr" = c(1,1,1),   # zero variance
    "Var1_151Eu"  = c(0,2,4),   # variable
    "Var2_155Gd"  = c(3,1,5),   # variable (ensures >=2 features remain after drop)
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # Keep zero-var: expect all three friendly features present
  obj_keep <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    drop_zero_var_features = FALSE, drop_na_features = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  feats_keep <- rownames(SeuratObject::LayerData(obj_keep[["IMC"]], "counts"))
  expect_setequal(feats_keep, c("141Pr", "151Eu", "155Gd"))
  
  # Drop zero-var: expect the constant feature removed, two variables kept
  obj_drop <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    drop_zero_var_features = TRUE, drop_na_features = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))[[1]]
  feats_drop <- rownames(SeuratObject::LayerData(obj_drop[["IMC"]], "counts"))
  expect_setequal(feats_drop, c("151Eu", "155Gd"))
})


test_that("feature_name_fun is applied to friendly channel names", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  cs_leaf <- file.path(top, "CellSegmentation", "ProjF", "FOVF")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    cell_id = 1:2, centroid_x = c(0,1), centroid_y = c(0,1),
    "Cd3e_141Pr" = c(1,2), "Ms4a1_151Eu" = c(3,4), check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  objs <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = FALSE,
    feature_name_fun = function(x) paste0("fn_", toupper(x)),
    add_asinh_layer = FALSE, build_segmentation_from_mask = FALSE, verbose = FALSE
  ))
  feats <- rownames(SeuratObject::LayerData(objs[[1]][["IMC"]], "counts"))
  
  # Treat "_" and "-" as equivalent separators
  feats_norm <- sub("^fn[-_]", "fn_", feats)
  expect_setequal(feats_norm, c("fn_141PR", "fn_151EU"))
})


test_that("min_cells filters features; min_features=1 is accepted when satisfied", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  cs_leaf <- file.path(top, "CellSegmentation", "ProjM", "FOVM")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  
  # A present in 1 cell (should drop at min.cells=2);
  # B and C present in BOTH cells and VARIABLE across cells (avoid zero-var drop)
  df <- data.frame(
    cell_id    = 1:2,
    centroid_x = c(0, 1),
    centroid_y = c(0, 1),
    "A_141Pr"  = c(1, 0),  # present in 1 cell → drop by min.cells=2
    "B_151Eu"  = c(2, 1),  # present in 2 cells, variance > 0
    "C_160Gd"  = c(3, 4),  # present in 2 cells, variance > 0
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # min.cells = 2 → drop A, keep B & C; keep zero-var filtering OFF to isolate min.cells behavior
  objA <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels   = FALSE,
    make_sparse               = FALSE,
    add_asinh_layer           = FALSE,
    min_cells                 = 2,
    min_features              = 0,
    drop_zero_var_features    = FALSE,  # <- important
    build_segmentation_from_mask = FALSE,
    verbose                   = FALSE
  ))[[1]]
  featsA <- rownames(SeuratObject::LayerData(objA[["IMC"]], "counts"))
  expect_setequal(featsA, c("151Eu", "160Gd"))
  
  # min.features = 1 → both cells retained (each has ≥ 1 detected feature)
  objB <- suppressWarnings(LoadIMCSegmented(
    top_dir = top,
    keep_unlabeled_channels   = FALSE,
    make_sparse               = FALSE,
    add_asinh_layer           = FALSE,
    min_cells                 = 0,
    min_features              = 1,
    drop_zero_var_features    = FALSE,  # <- important
    build_segmentation_from_mask = FALSE,
    verbose                   = FALSE
  ))[[1]]
  expect_equal(ncol(SeuratObject::LayerData(objB[["IMC"]], "counts")), 2L)
})



test_that("spatial_rotate_k=0 is identity; k=2 rotates 180°; 1 vs 1L identical", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  cs_leaf <- file.path(top, "CellSegmentation", "ProjR", "FOVR")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  # Two cells along x-axis: midrange center is (1, 0)
  df <- data.frame(
    cell_id = 1:2, centroid_x = c(0, 2), centroid_y = c(0, 0),
    "A_141Pr" = c(1, 2),  # var > 0 (avoid feature drop)
    "B_151Eu" = c(2, 3),  # var > 0
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  cell_names <- paste("FOVR", df$cell_id, sep = "_")
  
  # k = 0 (identity)
  obj0 <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE, make_sparse = FALSE,
    spatial_rotate_k = 0L
  ))[[1]]
  emb0 <- SeuratObject::Embeddings(obj0[["spatial"]])[cell_names, , drop = FALSE]
  expect_equal(unname(emb0[, "SP_1"]), df$centroid_x)
  expect_equal(unname(emb0[, "SP_2"]), df$centroid_y)
  
  # k = 2 (180°) around midrange center
  obj2 <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE, make_sparse = FALSE,
    spatial_rotate_k = 2L
  ))[[1]]
  emb2 <- SeuratObject::Embeddings(obj2[["spatial"]])[cell_names, , drop = FALSE]
  cx <- mean(range(df$centroid_x)); cy <- mean(range(df$centroid_y))
  expect_equal(unname(emb2[, "SP_1"]), 2*cx - df$centroid_x, tolerance = 1e-10)
  expect_equal(unname(emb2[, "SP_2"]), 2*cy - df$centroid_y, tolerance = 1e-10)
  
  # numeric 1 vs integer 1L produce identical embeddings
  obj1d <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE, make_sparse = FALSE,
    spatial_rotate_k = 1          # numeric
  ))[[1]]
  obj1i <- suppressWarnings(LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, add_asinh_layer = FALSE,
    build_segmentation_from_mask = FALSE, verbose = FALSE, make_sparse = FALSE,
    spatial_rotate_k = 1L         # integer
  ))[[1]]
  emb1d <- SeuratObject::Embeddings(obj1d[["spatial"]])[cell_names, , drop = FALSE]
  emb1i <- SeuratObject::Embeddings(obj1i[["spatial"]])[cell_names, , drop = FALSE]
  expect_equal(emb1d, emb1i, tolerance = 1e-10)
})


test_that("segmentation builds for seg_simplify_tol = 0 and higher values; centroid_nsides accepts various values", {
  skip_on_cran()
  skip_if_not_installed("withr"); skip_if_not_installed("rlang")
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  skip_if_not_installed("terra");  skip_if_not_installed("sf")
  skip_if(utils::packageVersion("SeuratObject") < "5.0.0")
  withr::local_tempdir() -> top
  
  # Workspace with 3 cells
  cs_leaf <- file.path(top, "CellSegmentation", "ProjSeg", "FOVSeg")
  dir.create(cs_leaf, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    cell_id = 1:3,
    centroid_x = c(10, 35, 15),
    centroid_y = c(10, 35, 50),
    "A_141Pr" = c(1,2,3),
    "B_151Eu" = c(3,2,1),
    check.names = FALSE
  )
  utils::write.csv(df, file.path(cs_leaf, "SingleCellObjects.csv"), row.names = FALSE)
  
  # Simple label map with three distinct regions
  m <- matrix(0L, nrow = 64, ncol = 64)
  m[ 5:15,  5:15] <- 1L
  m[25:35, 25:35] <- 2L
  m[45:60, 10:20] <- 3L
  r <- terra::rast(m); terra::writeRaster(r, file.path(cs_leaf, "Cell_mask.tiff"), overwrite = TRUE, datatype = "INT2U")
  
  # seg_simplify_tol = 0
  obj0 <- LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = TRUE, seg_simplify_tol = 0,
    centroid_nsides = 6, default_boundary = "segmentation", verbose = FALSE
  )[[1]]
  fov0 <- SeuratObject::DefaultFOV(obj0)
  expect_true(!is.null(obj0[[fov0]][["segmentation"]]))
  
  # seg_simplify_tol = 4 (still builds)
  obj4 <- LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = TRUE, seg_simplify_tol = 4,
    centroid_nsides = 6, default_boundary = "segmentation", verbose = FALSE
  )[[1]]
  fov4 <- SeuratObject::DefaultFOV(obj4)
  expect_true(!is.null(obj4[[fov4]][["segmentation"]]))
  
  # centroid_nsides = 3 and 8 (centroid boundary still present and usable)
  obj_tri <- LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = TRUE, seg_simplify_tol = 2,
    centroid_nsides = 3, default_boundary = "segmentation", verbose = FALSE
  )[[1]]
  obj_oct <- LoadIMCSegmented(
    top_dir = top, keep_unlabeled_channels = FALSE, make_sparse = TRUE,
    add_asinh_layer = FALSE, build_segmentation_from_mask = TRUE, seg_simplify_tol = 2,
    centroid_nsides = 8, default_boundary = "segmentation", verbose = FALSE
  )[[1]]
  
  # In both cases, 'centroids' boundary exists (we don't assert geometry; rendering happens at plot time)
  expect_true("centroids" %in% SeuratObject::Boundaries(obj_tri[[SeuratObject::DefaultFOV(obj_tri)]]))
  expect_true("centroids" %in% SeuratObject::Boundaries(obj_oct[[SeuratObject::DefaultFOV(obj_oct)]]))
})

