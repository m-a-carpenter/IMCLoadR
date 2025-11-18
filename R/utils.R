# R/utils-imc.R

#' Internal helpers for IMC wrappers
#' @keywords internal
#' @noRd

# Null-coalescing for scalars
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Conditional message
.msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

# Fast CSV (data.table if available; otherwise base utils)
.read_csv_fast <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# Pull attributes from a terra Spat* object; robust to missing/odd extents
.get_attrs <- function(v) {
  out <- try(terra::values(v), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  a <- terra::as.data.frame(v, geom = "XY")
  a[, setdiff(names(a), c("x","y","part","hole","piece","id")), drop = FALSE]
}

# Heuristic: does a mask look like a labelmap (many integerish unique labels)?
.mask_is_labelmap <- function(mask_path, n_cells, max_check = 2e6, verbose = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) return(FALSE)
  if (!is.character(mask_path) || is.na(mask_path) || !file.exists(mask_path)) return(FALSE)
  r <- suppressWarnings(try(terra::rast(mask_path), silent = TRUE))
  if (inherits(r, "try-error")) return(FALSE)
  vals <- suppressWarnings(try(terra::values(r, mat = FALSE), silent = TRUE))
  if (inherits(vals, "try-error") || is.null(vals)) return(FALSE)
  if (!is.null(max_check) && length(vals) > max_check) {
    set.seed(1L); vals <- vals[sample.int(length(vals), max_check)]
  }
  vals <- vals[is.finite(vals)]; if (!length(vals)) return(FALSE)
  nunq <- length(unique(round(vals)))
  frac_int <- mean(abs(vals - round(vals)) < 1e-6)
  ncells_r <- try(terra::ncell(r), silent = TRUE); if (inherits(ncells_r,"try-error")) ncells_r <- NA_integer_
  ok_count <- nunq >= 0.5 * n_cells && nunq <= 3 * n_cells
  ok_int   <- frac_int > 0.99
  ok_dense <- is.na(ncells_r) || nunq < 0.25 * ncells_r
  is_lab <- isTRUE(ok_count && ok_int && ok_dense)
  .msg(verbose, "Mask '%s': unique=%s, cells=%s, integerish=%.3f, dense_ok=%s -> labelmap=%s",
       basename(mask_path), nunq, n_cells, frac_int, ok_dense, is_lab)
  is_lab
}

# 90 deg k rotation about a center
.rotate_xy <- function(x, y, k = 0L, cx = NULL, cy = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  k <- as.integer(k) %% 4L
  if (k == 0L) return(list(x = x, y = y))
  if (is.null(cx) || !is.finite(cx)) cx <- mean(range(x, na.rm = TRUE))
  if (is.null(cy) || !is.finite(cy)) cy <- mean(range(y, na.rm = TRUE))
  xr <- x - cx; yr <- y - cy
  if      (k == 1L) { x <-  yr + cx; y <- -xr + cy }
  else if (k == 2L) { x <- -xr + cx; y <- -yr + cy }
  else if (k == 3L) { x <- -yr + cx; y <-  xr + cy }
  list(x = x, y = y)
}

# asinh transform that preserves sparsity for dgCMatrix
.sparse_asinh <- function(M, cofactor = 5) {
  if (inherits(M, "dgCMatrix")) { M2 <- M; M2@x <- asinh(M2@x / cofactor); return(M2) }
  asinh(M / cofactor)
}

# Read OME-XML pixel size (um/pixel), if available
.read_ome_pixel_size <- function(ome_xml_path) {
  if (!requireNamespace("xml2", quietly = TRUE)) return(NULL)
  if (!is.character(ome_xml_path) || !nzchar(ome_xml_path) || !file.exists(ome_xml_path)) return(NULL)
  x <- xml2::read_xml(ome_xml_path)
  px <- suppressWarnings(as.numeric(xml2::xml_attr(xml2::xml_find_first(x, ".//Pixels"), "PhysicalSizeX")))
  py <- suppressWarnings(as.numeric(xml2::xml_attr(xml2::xml_find_first(x, ".//Pixels"), "PhysicalSizeY")))
  if (is.finite(px) && is.finite(py)) list(um_per_pixel_x = px, um_per_pixel_y = py) else NULL
}

# Ensure a terra raster has a valid extent (some export tools omit it)
.ensure_raster_extent <- function(r) {
  ex_ok <- FALSE
  ex <- try(terra::ext(r), silent = TRUE)
  if (!inherits(ex, "try-error")) {
    vals <- c(
      tryCatch(terra::xmin(ex), error = function(...) NA_real_),
      tryCatch(terra::xmax(ex), error = function(...) NA_real_),
      tryCatch(terra::ymin(ex), error = function(...) NA_real_),
      tryCatch(terra::ymax(ex), error = function(...) NA_real_)
    )
    ex_ok <- all(is.finite(vals))
  }
  if (!ex_ok) terra::ext(r) <- terra::ext(0, terra::ncol(r), 0, terra::nrow(r))
  r
}

# Build a SeuratObject Segmentation boundary from a label mask, aligned to a target bbox
.mask_to_segmentation_aligned <- function(mask_path,
                                          cell_names,
                                          target_bbox,         # list(xr=c(min,max), yr=c(min,max))
                                          simplify_tol = 0,
                                          verbose = TRUE,
                                          rotate_k_first = 0L,
                                          cx_first = NULL, cy_first = NULL) {
  if (!requireNamespace("terra", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    if (verbose) message("terra/sf not available; skipping segmentation build."); return(NULL)
  }
  if (!is.character(mask_path) || is.na(mask_path) || !file.exists(mask_path)) {
    if (verbose) message("Cell mask not found; skipping segmentation build."); return(NULL)
  }
  
  ids_num <- suppressWarnings(as.integer(sub("^.*_", "", cell_names)))
  keep <- which(is.finite(ids_num)); if (!length(keep)) return(NULL)
  ids_num <- unique(ids_num[keep]); id2cell <- setNames(cell_names[keep], as.character(ids_num))
  
  r <- suppressWarnings(try(terra::rast(mask_path), silent = TRUE))
  if (inherits(r, "try-error")) return(NULL)
  ex_try <- try(terra::ext(r), silent = TRUE)
  if (inherits(ex_try, "try-error")) {
    terra::ext(r) <- terra::ext(0, terra::ncol(r), 0, terra::nrow(r))
  } else {
    b <- c(tryCatch(terra::xmin(ex_try), error=function(...) NA_real_),
           tryCatch(terra::xmax(ex_try), error=function(...) NA_real_),
           tryCatch(terra::ymin(ex_try), error=function(...) NA_real_),
           tryCatch(terra::ymax(ex_try), error=function(...) NA_real_))
    if (any(!is.finite(b))) terra::ext(r) <- terra::ext(0, terra::ncol(r), 0, terra::nrow(r))
  }
  
  v <- terra::as.polygons(r, dissolve = TRUE)
  if (terra::nrow(v) == 0) { if (verbose) message("Mask polygonization returned no polygons."); return(NULL) }
  attrs <- .get_attrs(v); if (!ncol(attrs)) { if (verbose) message("No attributes on polygons."); return(NULL) }
  num_cols <- names(attrs)[vapply(attrs, is.numeric, TRUE)]
  id_col <- if (length(num_cols)) num_cols[1] else names(attrs)[1]
  vals_num <- suppressWarnings(as.integer(round(as.numeric(attrs[[id_col]]))))
  keep_rows <- which(is.finite(vals_num) & vals_num %in% ids_num)
  if (!length(keep_rows)) { if (verbose) message("No overlapping labels between mask and cell IDs."); return(NULL) }
  v <- v[keep_rows, ]; vals_num <- vals_num[keep_rows]
  
  s <- sf::st_as_sf(v)
  if (simplify_tol > 0) s$geometry <- sf::st_simplify(s$geometry, dTolerance = simplify_tol, preserveTopology = TRUE)
  
  coords_list <- vector("list", nrow(s))
  for (i in seq_len(nrow(s))) {
    lab <- vals_num[i]; cell <- id2cell[[as.character(lab)]]
    if (is.null(cell)) next
    mat <- sf::st_coordinates(s$geometry[i])
    if (!is.matrix(mat) || ncol(mat) < 2) next
    if ("L2" %in% colnames(mat)) mat <- mat[mat[, "L2"] == 1, , drop = FALSE]
    coords_list[[i]] <- data.frame(cell = cell, x = mat[, "X"], y = mat[, "Y"], row.names = NULL)
  }
  coords_df <- do.call(rbind, coords_list)
  if (is.null(coords_df) || !nrow(coords_df)) { if (verbose) message("Polygon extraction yielded no coordinates."); return(NULL) }
  
  # rotate polygons (same center & k as centroids) BEFORE alignment
  rot <- .rotate_xy(coords_df$x, coords_df$y, k = rotate_k_first, cx = cx_first, cy = cy_first)
  coords_df$x <- rot$x; coords_df$y <- rot$y
  
  # align segmentation bbox to target (centroids) bbox
  xr_s <- range(coords_df$x, na.rm = TRUE); yr_s <- range(coords_df$y, na.rm = TRUE)
  xr_t <- target_bbox$xr;                 yr_t <- target_bbox$yr
  sx <- (xr_t[2] - xr_t[1]) / max(xr_s[2] - xr_s[1], .Machine$double.eps)
  sy <- (yr_t[2] - yr_t[1]) / max(yr_s[2] - yr_s[1], .Machine$double.eps)
  coords_df$x <- (coords_df$x - xr_s[1]) * sx + xr_t[1]
  coords_df$y <- (coords_df$y - yr_s[1]) * sy + yr_t[1]
  
  seg <- SeuratObject::CreateSegmentation(coords = coords_df)
  attr(seg, "coords_df") <- coords_df
  seg
}

# Get segmentation coords for an object:
.get_seg_df_for_obj <- function(obj, mask_path, rebuild = TRUE, verbose = TRUE) {
  fov_key <- tryCatch(SeuratObject::DefaultFOV(obj), error = function(e) NULL)
  
  # 1) named tool key first
  seg_df <- NULL
  tool_named <- suppressWarnings(try(SeuratObject::Tool(obj, "LoadIMCSegmented"), silent = TRUE))
  if (!inherits(tool_named, "try-error") && is.list(tool_named)) {
    seg_list <- tool_named$segmentation_coords
    if (is.list(seg_list) && !is.null(fov_key)) {
      cand <- seg_list[[fov_key]]
      if (is.data.frame(cand) && all(c("cell","x","y") %in% names(cand))) {
        return(cand[, c("cell","x","y")])
      }
    }
  }
  
  # 2) flat Tool(obj), nested "LoadIMCSegmented"
  tool_flat <- suppressWarnings(try(SeuratObject::Tool(obj), silent = TRUE))
  if (!inherits(tool_flat, "try-error") && is.list(tool_flat)) {
    seg_list <- tool_flat$LoadIMCSegmented$segmentation_coords
    if (is.list(seg_list) && !is.null(fov_key)) {
      cand <- seg_list[[fov_key]]
      if (is.data.frame(cand) && all(c("cell","x","y") %in% names(cand)))
        return(cand[, c("cell","x","y")])
    }
    # 2b) very old flat location
    seg_list <- tool_flat$segmentation_coords
    if (is.list(seg_list) && !is.null(fov_key)) {
      cand <- seg_list[[fov_key]]
      if (is.data.frame(cand) && all(c("cell","x","y") %in% names(cand)))
        return(cand[, c("cell","x","y")])
    }
  }
  
  # 3) boundary stored on FOV (if present)
  if (!is.null(fov_key)) {
    seg_df <- try(as.data.frame(obj[[fov_key]][["segmentation"]]), silent = TRUE)
    if (!inherits(seg_df, "try-error") && is.data.frame(seg_df)) {
      nm <- names(seg_df)
      if (!"cell" %in% nm) {
        cand <- intersect(c("cell","id","cell_id","Cell","ID","name"), nm)
        if (length(cand)) names(seg_df)[match(cand[1], nm)] <- "cell"
      }
      if (all(c("cell","x","y") %in% names(seg_df))) return(seg_df[, c("cell","x","y")])
    }
  }
  
  # 4) rebuild from mask if allowed (terra/sf)
  if (!rebuild) return(NULL)
  if (!requireNamespace("terra", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    if (isTRUE(verbose)) message("terra/sf not available; cannot rebuild segmentation.")
    return(NULL)
  }
  if (!is.character(mask_path) || is.na(mask_path) || !file.exists(mask_path)) return(NULL)
  
  r <- suppressWarnings(try(terra::rast(mask_path), silent = TRUE))
  if (inherits(r, "try-error")) return(NULL)
  r <- .ensure_raster_extent(r)
  v <- terra::as.polygons(r, dissolve = TRUE)
  if (terra::nrow(v) == 0) return(NULL)
  
  attrs <- .get_attrs(v); if (!ncol(attrs)) return(NULL)
  num_cols <- names(attrs)[vapply(attrs, is.numeric, TRUE)]
  id_col <- if (length(num_cols)) num_cols[1] else names(attrs)[1]
  vals_num <- suppressWarnings(as.integer(round(as.numeric(attrs[[id_col]]))))
  
  ids_num <- suppressWarnings(as.integer(sub("^.*_", "", colnames(obj))))
  ids_num <- unique(ids_num[is.finite(ids_num)])
  keep_rows <- which(is.finite(vals_num) & vals_num %in% ids_num)
  if (!length(keep_rows)) return(NULL)
  v <- v[keep_rows, ]; vals_num <- vals_num[keep_rows]
  
  s <- sf::st_as_sf(v)
  id2cell <- setNames(colnames(obj), as.character(suppressWarnings(as.integer(sub("^.*_", "", colnames(obj))))))
  lst <- vector("list", nrow(s))
  for (i in seq_len(nrow(s))) {
    lab <- vals_num[i]; cell <- id2cell[[as.character(lab)]]
    if (is.null(cell)) next
    mat <- sf::st_coordinates(s$geometry[i])
    if (!is.matrix(mat) || ncol(mat) < 2) next
    if ("L2" %in% colnames(mat)) mat <- mat[mat[, "L2"] == 1, , drop = FALSE]
    lst[[i]] <- data.frame(cell = cell, x = mat[, "X"], y = mat[, "Y"])
  }
  seg_df <- do.call(rbind, lst)
  if (!is.null(seg_df) && nrow(seg_df)) seg_df[, c("cell","x","y")] else NULL
}
