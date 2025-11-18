if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("SP_1","SP_2","x","y","cell","ftr","val","alpha_val"))
}

#' Build Seurat (v5) objects from MCDV IMC exports (per-FOV, no merge)
#'
#' @description
#' Ingest a top-level MCDV workspace and build **one Seurat v5 object per FOV**.
#' The function minimizes redundancy (counts as the default layer, optional \code{asinh}),
#' attaches spatial information, and (optionally) reconstructs segmentation polygons
#' from the label mask, aligned to the centroid coordinate frame.
#'
#' @param top_dir Character. Path to the MCDV workspace (e.g., \code{"Workspace_008"}).
#' @param keep_unlabeled_channels Logical. Keep channels without "_" in their name
#'   (default \code{FALSE} drops unlabeled channels).
#' @param assay Character. Assay name to create (default \code{"IMC"}).
#' @param make_sparse Logical. Store expression as a sparse matrix when possible (default \code{TRUE}).
#' @param verbose Logical. Emit progress messages (default \code{TRUE}).
#'
#' @section Layers & normalization:
#' The loader creates a `counts` layer and can derive an `asinh` layer (cofactor configurable).
#' @param add_asinh_layer Logical. Add an \code{asinh} layer derived from \code{counts} (default \code{TRUE}).
#' @param asinh_cofactor Numeric. Cofactor for \code{asinh} (default \code{5}).
#' @param put_matrix_in_data_layer Logical. Also copy the input matrix to a Seurat \code{data} layer (default \code{FALSE}).
#'
#' @section Feature hygiene:
#' Logical options for dropping NA and/or 0 variance features, applying user specified naming conditions, and setting minimum cells or features required.
#' @param drop_na_features Logical. Drop features containing any NA (default \code{TRUE}).
#' @param drop_zero_var_features Logical. Drop zero-variance features (default \code{TRUE}).
#' @param feature_name_fun Optional function applied to **friendly** feature names (after channel mapping).
#' @param min_cells,min_features Passed through to \code{Seurat::CreateSeuratObject} (defaults 0, 0).
#'
#' @section Orientation:
#' Optional rotation of image to designate tissue orientation alignment.
#' @param spatial_rotate_k Integer in \{0,1,2,3\}. Rotate by 90°k around the FOV center (default \code{0}).
#'
#' @section Segmentation (optional):
#' Segmentation parameters.
#' @param build_segmentation_from_mask Logical. Rebuild polygons from \code{Cell_mask.tiff} if it looks like a label map (default \code{TRUE}).
#' @param seg_simplify_tol Numeric. dTolerance (pixels) for polygon simplification (default \code{2}).
#' @param centroid_nsides Numeric. Number of sides for centroid glyphs (default \code{6}; use \code{Inf} for circles).
#' @param default_boundary One of \code{"centroids"} or \code{"segmentation"}; which to set as active on the FOV (default \code{"centroids"}).
#'
#' @section Output & provenance:
#' Output options for saving generated objects.
#' @param save_h5seurat Optional file path template; if non-NULL, save each object as H5Seurat
#'   using \code{SeuratDisk}, with \code{\{FOV\}} substituted by the FOV name (e.g., \code{"~/out/\{FOV\}.h5seurat"}).
#'
#' @return A named \code{list} of Seurat objects (one per FOV). Each object contains:
#' \itemize{
#'   \item an \code{IMC} assay with \code{counts} (and optional \code{asinh}, \code{data})
#'   \item a \code{spatial} reduction (\code{SP_1}, \code{SP_2}) rotated by 90°×k
#'   \item a \code{FOV} with a \code{centroids} boundary; optional \code{segmentation}
#'   \item tool data accessible via \code{SeuratObject::Tool(obj)}
#'         (entry keyed to this function) containing \code{paths}, \code{feature_meta},
#'         \code{scale_factors}, \code{spatial_orientation}, \code{transforms}, \code{provenance},
#'         and optional \code{segmentation_coords}.
#' }
#' @author Meredith A. Carpenter
#' @references Giesen C et al. (2014). Highly multiplexed imaging of tumor tissues with subcellular resolution by mass cytometry. Cell 159(7):1660–1675. doi:10.1016/j.cell.2014.11.023
#' @references Keren L et al. (2018). A structured tumor-immune microenvironment revealed by multiplexed ion beam imaging. Cell 174(6):1373–1387.e19. doi:10.1016/j.cell.2018.08.039
#' @references Le Rochais M et al. (2022). Application of High-Throughput Imaging Mass Cytometry Hyperion in Cancer Research. Front Immunol 13:859414. doi:10.3389/fimmu.2022.859414
#' @references Fluidigm Corporation (2019). MCD Viewer v1.0.560.6 User Guide (PN 400317-4, Rev. 3). https://www.imc.unibe.ch/unibe/portal/fak_medizin/micro_imc/content/e987276/e1000503/e1000513/MCDViewerv1.0.560.6UserGuide.pdf
#'
#' @examples
#' \dontrun{
#' objs <- load_segmented_IMC(
#'   top_dir = "~/inst/extdata/Workspace_008",
#'   keep_unlabeled_channels = FALSE,
#'   add_asinh_layer = TRUE,
#'   spatial_rotate_k = 1L,  # rotate 90 degrees
#'   build_segmentation_from_mask = TRUE,
#'   default_boundary = "segmentation"
#' )
#'
#' first <- objs[[1]]
#' DefaultAssay(first)
#' SeuratObject::Images(first)
#' ImageDimPlot(first, fov = SeuratObject::DefaultFOV(first), boundaries = "segmentation")
#' }
#'
#' @seealso \code{\link{IMCFovPlot}}
#' @export
#' @importFrom utils glob2rx tail
#' @importFrom stats var ave setNames
#' @importFrom grDevices as.raster
LoadIMCSegmented <- function(
    top_dir,
    keep_unlabeled_channels = FALSE,
    assay = "IMC",
    make_sparse = TRUE,
    verbose = TRUE,
    # layers
    add_asinh_layer = TRUE,
    asinh_cofactor  = 5,
    put_matrix_in_data_layer = FALSE,
    # feature hygiene
    drop_na_features = TRUE,
    drop_zero_var_features = TRUE,
    feature_name_fun = NULL,
    min_cells = 0, min_features = 0,
    # orientation (rotation only)
    spatial_rotate_k = 0L,
    # segmentation
    build_segmentation_from_mask = TRUE,
    seg_simplify_tol = 2,
    default_boundary = c("centroids","segmentation"),
    centroid_nsides = 6,
    # output
    save_h5seurat = NULL
) {
  default_boundary <- match.arg(default_boundary)
  if (!dir.exists(top_dir)) stop("Missing top_dir: ", top_dir)

  if (!requireNamespace("rlang", quietly = TRUE))
    stop("Please install 'rlang' to run dependency checks.", call. = FALSE)

  rlang::check_installed(c("Seurat", "SeuratObject"))
  if (isTRUE(build_segmentation_from_mask))
    rlang::check_installed(c("terra", "sf"))

  # ----- paths & discovery -----
  cs_root   <- file.path(top_dir, "CellSegmentation")
  tiff_root <- file.path(top_dir, "tiff_multi_page_export")
  if (!dir.exists(cs_root))  stop("Missing directory: ", cs_root)
  if (!dir.exists(tiff_root)) .msg(TRUE, "Warning: missing %s (image paths may be NA)", tiff_root)

  csv_files <- list.files(cs_root, pattern = "^SingleCellObjects\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(csv_files)) stop("No SingleCellObjects.csv files found under: ", cs_root)
  .msg(verbose, "Found %d FOV csv(s).", length(csv_files))

  find_mask_tiff <- function(project, fov) {
    p <- file.path(cs_root, project, fov, "Cell_mask.tiff"); if (file.exists(p)) p else NA_character_
  }
  find_panorama_png <- function(project, fov) {
    proj_dir <- file.path(tiff_root, project); if (!dir.exists(proj_dir)) return(NA_character_)
    patt <- utils::glob2rx(paste0("*", fov, "*Panorama*.png"))
    hits <- list.files(proj_dir, pattern = patt, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) hits[[1]] else NA_character_
  }
  find_ome_tiff <- function(project, fov) {
    p <- file.path(tiff_root, project, fov, paste0(fov, ".ome.tiff"))
    if (file.exists(p)) return(p)
    proj_dir <- file.path(tiff_root, project); if (!dir.exists(proj_dir)) return(NA_character_)
    patt <- utils::glob2rx(paste0(fov, ".ome.tiff"))
    hits <- list.files(proj_dir, pattern = patt, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) hits[[1]] else NA_character_
  }

  # ----- global channel map -----
  all_feature_cols <- character(0)
  for (f in csv_files) {
    df <- .read_csv_fast(f)
    area_col <- names(df)[grepl("^cell_area", names(df), ignore.case = TRUE)]
    meta_cols <- unique(c("cell_id","centroid_x","centroid_y", area_col))
    feat_cols <- setdiff(names(df), meta_cols[meta_cols %in% names(df)])
    all_feature_cols <- union(all_feature_cols, feat_cols)
  }
  labeled_all   <- all_feature_cols[grepl("_", all_feature_cols, fixed = TRUE)]
  unlabeled_all <- setdiff(all_feature_cols, labeled_all)
  if (!length(labeled_all) && !keep_unlabeled_channels) {
    stop("No labeled channels (with '_') found globally, and keep_unlabeled_channels = FALSE.")
  }
  friendly_all <- make.unique(sub("^[^_]*_", "", labeled_all), sep = "___dup")
  if (!is.null(feature_name_fun)) friendly_all <- vapply(friendly_all, feature_name_fun, "", USE.NAMES = FALSE)
  friendly_map <- setNames(friendly_all, labeled_all)
  if (keep_unlabeled_channels && length(unlabeled_all)) {
    friendly_map <- c(friendly_map, setNames(unlabeled_all, unlabeled_all))
  }

  # ----- per-FOV build -----
  obj_list <- list()

  for (f in csv_files) {
    df <- .read_csv_fast(f)
    fov     <- basename(dirname(f))
    project <- basename(dirname(dirname(f)))

    # paths & mask info early (for rotation center)
    mask_path <- find_mask_tiff(project, fov)
    px_meta  <- tryCatch({
      # optional dependency; only use if available
      if (!is.na(mask_path) && file.exists(mask_path) && requireNamespace("magick", quietly = TRUE)) {
        magick::image_info(magick::image_read(mask_path))[c("width","height")]
      } else NULL
    }, error = function(...) NULL)
    cx0 <- if (!is.null(px_meta)) px_meta$width/2  else mean(range(df$centroid_x, na.rm = TRUE))
    cy0 <- if (!is.null(px_meta)) px_meta$height/2 else mean(range(df$centroid_y, na.rm = TRUE))

    area_col <- names(df)[grepl("^cell_area", names(df), ignore.case = TRUE)]
    meta_cols <- unique(c("cell_id","centroid_x","centroid_y", area_col))
    missing_meta <- setdiff(c("cell_id","centroid_x","centroid_y"), names(df))
    if (length(missing_meta)) stop("Missing required meta column(s) in ", f, ": ", paste(missing_meta, collapse=", "))

    feat_here <- setdiff(names(df), meta_cols)
    kept_here <- intersect(feat_here, names(friendly_map))
    if (!length(kept_here)) { .msg(verbose, "Skipping %s: no kept channels after filtering.", fov); next }

    expr <- df[, kept_here, drop = FALSE]
    colnames(expr) <- unname(friendly_map[kept_here])
    expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(x)))
    cell_names <- paste(fov, df$cell_id, sep = "_")
    rownames(expr) <- cell_names

    # Dense first (features x cells)
    expr_mat <- t(as.matrix(expr))

    # Hygiene on dense
    if (isTRUE(drop_na_features) || isTRUE(drop_zero_var_features)) {
      keep <- rep(TRUE, nrow(expr_mat))
      if (isTRUE(drop_na_features))       keep <- keep & (rowSums(is.na(expr_mat)) == 0)
      if (isTRUE(drop_zero_var_features)) {
        v <- apply(expr_mat, 1, stats::var, na.rm = TRUE); v[!is.finite(v)] <- 0; keep <- keep & (v > 0)
      }
      keep[is.na(keep)] <- FALSE
      expr_mat <- expr_mat[keep, , drop = FALSE]
    }

    # Optional sparse
    if (isTRUE(make_sparse) && requireNamespace("Matrix", quietly = TRUE)) {
      expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
    }

    # Seurat object (v5 assay with layers)
    obj <- Seurat::CreateSeuratObject(
      counts = expr_mat, assay = assay, project = fov,
      min.cells = min_cells, min.features = min_features
    )

    # Optional data layer duplication
    if (isTRUE(put_matrix_in_data_layer)) {
      ass <- obj[[assay]]
      SeuratObject::LayerData(ass, layer = "data") <- expr_mat
      obj[[assay]] <- ass
    }

    # Optional asinh layer + record transform in Tool (later)
    .imc_tool <- list(transforms = NULL)   # accumulate tool data here
    if (isTRUE(add_asinh_layer)) {
      ass <- obj[[assay]]
      M <- SeuratObject::LayerData(ass, layer = "counts")
      SeuratObject::LayerData(ass, layer = "asinh") <- .sparse_asinh(M, cofactor = asinh_cofactor)
      obj[[assay]] <- ass
      .imc_tool$transforms <- list(asinh = list(source_layer = "counts", cofactor = asinh_cofactor))
    }

    # meta.data to AddMetaData (avoid direct slot writes)
    meta <- df[, intersect(names(df), meta_cols), drop = FALSE]
    rownames(meta) <- cell_names
    if (length(area_col) == 1) {
      names(meta)[names(meta) == area_col] <- "cell_area_um2"
    } else if (!"cell_area_um2" %in% names(meta)) {
      meta$cell_area_um2 <- NA_real_
    }
    meta$fov <- fov; meta$project <- project
    obj <- SeuratObject::AddMetaData(obj, metadata = meta[colnames(obj), , drop = FALSE])

    # feature meta (store in Tool later)
    feats_present <- rownames(obj[[assay]])
    feature_meta <- data.frame(
      feature = feats_present,
      original_channel = vapply(feats_present, function(fn) {
        inv <- names(friendly_map)[friendly_map == fn]
        if (length(inv) == 1) inv else if (length(inv) > 1) paste(inv, collapse=";") else NA_character_
      }, character(1)),
      friendly_name = feats_present,
      channel_group = ifelse(feats_present %in% unname(friendly_map[labeled_all]), "labeled", "unlabeled"),
      stringsAsFactors = FALSE,
      row.names = feats_present
    )

    # spatial reduction from raw centroids with ROTATION only
    rot <- .rotate_xy(df$centroid_x, df$centroid_y, k = spatial_rotate_k, cx = cx0, cy = cy0)
    emb <- cbind(SP_1 = rot$x, SP_2 = rot$y)
    rownames(emb) <- cell_names
    obj[["spatial"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "SP_", assay = assay)

    # FOV + centroids
    cent_coords <- data.frame(x = emb[, "SP_1"], y = emb[, "SP_2"], row.names = colnames(obj))
    cent <- SeuratObject::CreateCentroids(coords = cent_coords, nsides = centroid_nsides)
    safe_key <- gsub("[^A-Za-z0-9.]+", ".", fov)
    before <- tryCatch(SeuratObject::Images(obj), error = function(...) character(0))
    obj[[safe_key]] <- SeuratObject::CreateFOV(coords = cent, assay = assay, name = "centroids")
    after  <- tryCatch(SeuratObject::Images(obj), error = function(...) character(0))
    actual_key <- setdiff(after, before); if (length(actual_key) != 1L) actual_key <- tail(after, 1L)
    SeuratObject::DefaultFOV(obj) <- actual_key
    SeuratObject::DefaultBoundary(obj[[actual_key]]) <- "centroids"

    # Segmentation (rotate first, then bbox-align)
    seg <- NULL
    seg_coords_list <- NULL
    if (isTRUE(build_segmentation_from_mask) &&
        .mask_is_labelmap(mask_path, n_cells = ncol(obj), verbose = verbose)) {
      target_bbox <- list(xr = range(cent_coords$x, na.rm = TRUE),
                          yr = range(cent_coords$y, na.rm = TRUE))
      seg <- .mask_to_segmentation_aligned(
        mask_path,
        cell_names   = colnames(obj),
        target_bbox  = target_bbox,
        simplify_tol = seg_simplify_tol,
        verbose      = verbose,
        rotate_k_first = spatial_rotate_k,
        cx_first = cx0, cy_first = cy0
      )
      if (!is.null(seg)) {
        obj[[actual_key]][["segmentation"]] <- seg
        seg_df <- attr(seg, "coords_df")
        if (is.data.frame(seg_df)) {
          seg_coords_list <- list()
          seg_coords_list[[actual_key]] <- seg_df[, c("cell","x","y")]
        }
      } else if (isTRUE(verbose)) {
        message("Segmentation build failed; keeping centroids only for FOV ", fov)
      }
    } else if (isTRUE(verbose)) {
      message("Skipping segmentation for FOV ", fov, " (no label map or disabled).")
    }

    # paths + scale factors + provenance (store in Tool)
    ome_tiff <- find_ome_tiff(project, fov)
    ome_xml  <- if (!is.na(ome_tiff)) sub("\\.ome\\.tiff$", ".ome.xml", ome_tiff) else NA_character_
    tool_paths <- list(
      top_dir              = normalizePath(top_dir, winslash = "/", mustWork = FALSE),
      project              = project,
      fov                  = fov,
      single_cell_csv_path = f,
      cell_mask_path       = mask_path,
      panorama_png_path    = find_panorama_png(project, fov),
      ome_tiff_path        = ome_tiff
    )
    tool_scale <- list(
      pixel  = list(width = if (!is.null(px_meta)) px_meta$width else NA_real_,
                    height= if (!is.null(px_meta)) px_meta$height else NA_real_),
      micron = .read_ome_pixel_size(`%||%`(ome_xml, ""))
    )
    tool_orient <- list(
      rotate_k = as.integer(spatial_rotate_k),
      center   = list(cx = as.numeric(cx0), cy = as.numeric(cy0))
    )
    tool_prov <- list(
      loader = "load_segmented_IMC",
      seurat_version = as.character(utils::packageVersion("Seurat")),
      seuratobject_version = as.character(utils::packageVersion("SeuratObject")),
      terra_version = if (requireNamespace("terra", quietly = TRUE)) as.character(utils::packageVersion("terra")) else NA,
      sf_version    = if (requireNamespace("sf", quietly = TRUE))    as.character(utils::packageVersion("sf"))    else NA
    )

    # default boundary preference
    if (default_boundary == "segmentation" && !is.null(seg)) {
      SeuratObject::DefaultBoundary(obj[[actual_key]]) <- "segmentation"
    }

    # Save optional H5Seurat
    if (!is.null(save_h5seurat)) {
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        warning("save_h5seurat specified but SeuratDisk not installed; skipping save.")
      } else {
        path <- sub("\\{FOV\\}", fov, save_h5seurat, fixed = TRUE)
        dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
        SeuratDisk::SaveH5Seurat(obj, filename = path, overwrite = TRUE)
      }
    }

    # ---- write tool data (one call; entry auto-named to calling function) ----
    SeuratObject::Tool(obj) <- list(
      paths               = tool_paths,
      feature_meta        = setNames(list(feature_meta), assay),
      scale_factors       = tool_scale,
      spatial_orientation = tool_orient,
      transforms          = .imc_tool$transforms,
      provenance          = tool_prov,
      segmentation_coords = seg_coords_list
    )

    obj <- SeuratObject::LogSeuratCommand(obj)  # optional: log call for reproducibility
    obj_list[[fov]] <- obj
    .msg(verbose, "Built FOV '%s': %d cells, %d features.", fov, ncol(obj), nrow(obj))
  }

  if (!length(obj_list)) stop("No per-FOV Seurat v5 objects were created; check inputs and filters.")
  obj_list
}
