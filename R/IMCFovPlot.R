#' Plot feature overlays on an IMC FOV (Seurat v5)
#' 
#' @description
#' Overlays cell intensities onto the MCDV mask/background image using either
#' segmentation polygons or centroids. Supports composite modes
#' (\code{"max"}, \code{"mean"}, \code{"sum"}, \code{"rgb"}, \code{"stack"})
#' and an optional visualization-only rotation like the loader.
#'
#' @details
#' \strong{Data path:}
#' The matrix is taken from the requested \code{layer}. If that layer is absent,
#' the function chooses the first available layer from
#' \code{c("asinh","scale.data","data","counts")} and, when \code{verbose = TRUE},
#' emits a message indicating the fallback.
#'
#' \strong{Transform/normalize order:}
#' Values are transformed first (via \code{transform = "log1p"} or \code{"sqrt"}), then
#' optionally min–max normalized \emph{per feature across cells}
#' (\code{normalize = "minmax"}), and finally aggregated by \code{combine}.
#' For \code{combine = "rgb"}, each of the three features is independently
#' rescaled to 0,1 before composing \eqn{rgb(r,g,b, alpha)}.
#'
#' \strong{Mask discovery:}
#' The background image path (\code{cell_mask_path}) is discovered in this order:
#' \enumerate{
#'   \item \code{SeuratObject::Tool(obj, "LoadIMCSegmented")$paths$cell_mask_path}
#'   \item \code{SeuratObject::Tool(obj)$LoadIMCSegmented$paths$cell_mask_path} (nested legacy)
#'   \item \code{SeuratObject::Tool(obj)$paths$cell_mask_path} (flat legacy)
#'   \item \code{obj@misc$paths$cell_mask_path}
#'   \item \code{obj[[]]$cell_mask_path} (meta fallback; first non-NA)
#' }
#' The image dimensions determine the plotting canvas and the size used when exporting.
#'
#' \strong{Segmentation rebuild (optional):}
#' If \code{boundaries = "segmentation"} but polygons are missing, setting
#' \code{rebuild_seg_if_missing = TRUE} attempts to rebuild them from the mask.
#' This requires the \pkg{terra} and \pkg{sf} packages to be installed.
#'
#' \strong{Stack colours:}
#' For centroids, \code{combine = "stack"} uses \code{scale_color_manual()}; for segmentation
#' it uses \code{scale_fill_manual()}. Provide a \emph{named} vector \code{colors}
#' (names must match \code{features}) to guarantee the mapping.
#'
#' \strong{Rotation center:}
#' The rotation center is taken from \code{rotate_center} when provided;
#' otherwise from the loader tool's stored center; otherwise the image center.
#'
#' @param obj A Seurat object for a single FOV (as created by \code{\link{LoadIMCSegmented}}).
#' @param features Character vector of feature names in \code{assay} (friendly names produced by the loader, e.g., \code{"141Pr"}).
#' @param assay Assay to pull data from (default: \code{Seurat::DefaultAssay(obj)}).
#' @param layer Layer to use (\code{"asinh"}, \code{"scale.data"}, \code{"data"}, \code{"counts"}).
#' @param boundaries One of \code{"segmentation"} or \code{"centroids"}.
#' @param combine Composite mode: \code{"max"}, \code{"mean"}, \code{"sum"}, \code{"rgb"}, \code{"stack"}.
#' @param normalize Either \code{"none"} or \code{"minmax"} (per-feature).
#' @param transform Either \code{"none"}, \code{"log1p"}, \code{"sqrt"}.
#' @param pt_size Point size for centroid mode (ignored for segmentation polygons).
#' @param alpha Alpha for non-\code{"stack"} composites.
#' @param colors Optional colours for \code{"stack"} mode; a \emph{named} vector (names == \code{features}) is recommended.
#' @param alpha_range For \code{"stack"}: \code{c(min, max)} alpha; defaults to \code{c(0, alpha)}.
#' @param stack_rescale_alpha_by_n If \code{TRUE}, reduces the maximum alpha by \code{length(features)}.
#' @param stack_alpha_gamma Gamma exponent for alpha scaling in \code{"stack"}.
#' @param limits Optional numeric limits passed to the continuous colour/fill scale.
#' @param export_path Optional file path for \code{ggplot2::ggsave}; returns the plot invisibly after saving.
#'   The file type is inferred from the extension; output size is \code{(image_width / dpi) × (image_height / dpi)} inches.
#' @param dpi DPI for export when \code{export_path} is used.
#' @param verbose Logical; emit progress messages.
#' @param rebuild_seg_if_missing Try rebuilding segmentation from the mask if not found (needs \pkg{terra} and \pkg{sf}).
#' @param spatial_rotate_k Integer in \{0,1,2,3\}; rotate coordinates by 90°\eqn{\times k}.
#' @param rotate_center Optional numeric length‑2 vector \code{c(cx, cy)} for rotation center.
#' @return A \code{ggplot} object (printed normally; returned invisibly if \code{export_path} is used).
#'
#' @author Meredith A. Carpenter
#' @examples
#' \dontrun{
#' # Build per-FOV Seurat objects (requires Seurat v5, SeuratObject v5)
#' objs <- LoadIMCSegmented(
#'   top_dir = "~/Workspace_008",
#'   keep_unlabeled_channels = FALSE,
#'   add_asinh_layer = TRUE,
#'   build_segmentation_from_mask = TRUE,   # needs terra/sf if TRUE
#'   verbose = TRUE
#' )
#'
#' fov <- objs[[1]]
#' # Choose a couple of features from the IMC assay
#' feats2 <- rownames(fov[["IMC"]])[1:2]
#' feats3 <- rownames(fov[["IMC"]])[1:3]
#'
#' # Centroids overlay, continuous composite (max) from counts layer
#' p1 <- IMCFovPlot(
#'   fov,
#'   features   = feats2,
#'   assay      = "IMC",
#'   layer      = "counts",
#'   boundaries = "centroids",
#'   combine    = "max",
#'   normalize  = "minmax",
#'   spatial_rotate_k = 0L
#' )
#' print(p1)
#'
#' # Segmentation overlay (sum) from asinh layer, rotated 90 degrees
#' p2 <- IMCFovPlot(
#'   fov,
#'   features   = feats2,
#'   assay      = "IMC",
#'   layer      = "asinh",
#'   boundaries = "segmentation",
#'   combine    = "sum",
#'   spatial_rotate_k = 1L
#' )
#' print(p2)
#'
#' # RGB composite (exactly 3 features)
#' p3 <- IMCFovPlot(
#'   fov,
#'   features   = feats3,   # R,G,B in order
#'   assay      = "IMC",
#'   layer      = "asinh",
#'   boundaries = "centroids",
#'   combine    = "rgb"
#' )
#' print(p3)
#'
#' # Save a stack composite to disk
#' IMCFovPlot(
#'   fov,
#'   features   = feats2,
#'   assay      = "IMC",
#'   boundaries = "centroids",
#'   combine    = "stack",
#'   export_path = tempfile(fileext = ".png"),
#'   dpi = 300
#' )
#' }
#' @seealso \code{\link{LoadIMCSegmented}}
#' @export
#' @importFrom grid rasterGrob unit
#' @importFrom grDevices rgb
IMCFovPlot <- function(
    obj,
    features,
    assay      = Seurat::DefaultAssay(obj),
    layer      = c("asinh","scale.data","data","counts"),
    boundaries = c("segmentation","centroids"),
    combine    = c("max","mean","sum","rgb","stack"),
    normalize  = c("none","minmax"),
    transform  = c("none","log1p","sqrt"),
    pt_size    = 0.6,
    alpha      = 0.9,
    colors     = NULL,
    alpha_range = NULL,
    stack_rescale_alpha_by_n = FALSE,
    stack_alpha_gamma = 1,
    limits     = NULL,
    export_path = NULL,
    dpi        = 300,
    verbose    = TRUE,
    rebuild_seg_if_missing = TRUE,
    spatial_rotate_k = 0L,
    rotate_center = c(NA_real_, NA_real_)
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("magick",  quietly = TRUE)) stop("Package 'magick' is required.")
  
  layer      <- match.arg(layer)
  boundaries <- match.arg(boundaries)
  combine    <- match.arg(combine)
  normalize  <- match.arg(normalize)
  transform  <- match.arg(transform)
  
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
  if (missing(features) || is.null(features)) {
    features <- rownames(obj[[assay]])
    if (isTRUE(verbose)) message("No 'features' provided; using all (n=", length(features), ").")
  }
  features <- unique(features)
  if (!all(c("centroid_x","centroid_y") %in% colnames(obj[[]])))
    stop("meta.data must contain 'centroid_x' and 'centroid_y' (in obj[[]]).")
  
  # ---------------------------
  # Named tool read (preferred)
  # ---------------------------
  tool_named <- suppressWarnings(try(SeuratObject::Tool(obj, "LoadIMCSegmented"), silent = TRUE))
  if (inherits(tool_named, "try-error")) tool_named <- NULL
  
  # Back-compat: flat Tool(obj) (list) and legacy @tools / @misc
  tool_flat  <- suppressWarnings(try(SeuratObject::Tool(obj), silent = TRUE))
  if (inherits(tool_flat, "try-error") || !is.list(tool_flat)) tool_flat <- NULL
  
  # mask path: prefer named tool to flat tool to @misc to meta fallback
  mask_path <- NULL
  if (is.list(tool_named) && !is.null(tool_named$paths$cell_mask_path)) {
    mask_path <- tool_named$paths$cell_mask_path
  }
  if ((is.null(mask_path) || !file.exists(mask_path)) && is.list(tool_flat) &&
      !is.null(tool_flat$LoadIMCSegmented$paths$cell_mask_path)) {
    mask_path <- tool_flat$LoadIMCSegmented$paths$cell_mask_path
  }
  if ((is.null(mask_path) || !file.exists(mask_path)) && is.list(tool_flat) &&
      !is.null(tool_flat$paths$cell_mask_path)) {
    mask_path <- tool_flat$paths$cell_mask_path
  }
  if (is.null(mask_path) || !file.exists(mask_path)) {
    # legacy @misc
    mp <- suppressWarnings(try(obj@misc$paths$cell_mask_path, silent = TRUE))
    if (!inherits(mp, "try-error") && is.character(mp) && length(mp) == 1 && file.exists(mp))
      mask_path <- mp
  }
  if (is.null(mask_path) || !file.exists(mask_path)) {
    # per-cell meta fallback
    mp <- obj[[]]$cell_mask_path
    mp <- mp[!is.na(mp)]
    if (length(mp) && file.exists(mp[1])) {
      mask_path <- mp[1]
      if (isTRUE(verbose)) message("Using fallback meta.data$cell_mask_path: ", mask_path)
    } else {
      stop("Mask not found in Tool(obj, 'LoadIMCSegmented'), Tool(obj), obj@misc, or meta.")
    }
  }
  
  img <- magick::image_read(mask_path); ii <- magick::image_info(img)
  W <- ii$width; H <- ii$height
  
  # rotation center: prefer named tool to flat tool to image center
  cx_tool <- NA_real_; cy_tool <- NA_real_
  if (is.list(tool_named) && !is.null(tool_named$spatial_orientation$center)) {
    cx_tool <- suppressWarnings(as.numeric(tool_named$spatial_orientation$center$cx))
    cy_tool <- suppressWarnings(as.numeric(tool_named$spatial_orientation$center$cy))
  } else if (is.list(tool_flat)) {
    # nested first
    c_try <- tool_flat$LoadIMCSegmented$spatial_orientation$center
    if (!is.null(c_try)) {
      cx_tool <- suppressWarnings(as.numeric(c_try$cx))
      cy_tool <- suppressWarnings(as.numeric(c_try$cy))
    }
    # flat legacy
    if (!is.finite(cx_tool) || !is.finite(cy_tool)) {
      c_try <- tool_flat$spatial_orientation$center
      if (!is.null(c_try)) {
        if (!is.finite(cx_tool)) cx_tool <- suppressWarnings(as.numeric(c_try$cx))
        if (!is.finite(cy_tool)) cy_tool <- suppressWarnings(as.numeric(c_try$cy))
      }
    }
  }
  cx0 <- if (!is.na(rotate_center[1])) rotate_center[1] else { if (is.finite(cx_tool)) cx_tool else W/2 }
  cy0 <- if (!is.na(rotate_center[2])) rotate_center[2] else { if (is.finite(cy_tool)) cy_tool else H/2 }
  
  # choose layer and pull matrix
  ass <- obj[[assay]]
  present <- tryCatch(SeuratObject::Layers(ass), error = function(e) character(0))
  if (!(layer %in% present)) {
    pri <- c("asinh","scale.data","data","counts")
    alt <- pri[pri %in% present]
    if (!length(alt)) stop("No usable layers found in assay '", assay, "'.")
    if (isTRUE(verbose)) message("Layer '", layer, "' not found; using '", alt[1], "'.")
    layer <- alt[1]
  }
  M <- SeuratObject::LayerData(ass, layer = layer)
  if (!all(features %in% rownames(M))) {
    miss <- setdiff(features, rownames(M))
    stop("Some features missing from layer '", layer, "': ", paste(miss, collapse = ", "))
  }
  V <- as.matrix(M[features, colnames(obj), drop = FALSE])
  
  # optional transform & per-feature min-max
  if (transform == "log1p") V <- log1p(V)
  else if (transform == "sqrt") V <- sqrt(pmax(V, 0))
  
  if (normalize == "minmax" && combine %in% c("max","mean","sum","stack")) {
    rng <- apply(V, 1, function(x) {
      r <- range(x, na.rm = TRUE)
      if (!is.finite(r[2]) || r[2] <= r[1]) c(0,1) else r
    })
    for (i in seq_len(nrow(V))) {
      r1 <- rng[1,i]; r2 <- rng[2,i]
      if (is.finite(r2) && r2 > r1) V[i,] <- (V[i,] - r1) / (r2 - r1)
    }
  }
  
  # build background canvas
  img_rast <- as.raster(img)
  bg <- grid::rasterGrob(img_rast, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  p <- ggplot2::ggplot() +
    ggplot2::annotation_custom(bg, xmin = 0, xmax = W, ymin = 0, ymax = H) +
    ggplot2::coord_fixed(xlim = c(0, W), ylim = c(H, 0), expand = FALSE) +
    ggplot2::scale_y_continuous(trans = "reverse") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0, face = "bold"))
  
  # a nicer FOV title (named tool to flat tool to meta)
  fov_title <- ""
  if (is.list(tool_named) && !is.null(tool_named$paths$fov)) fov_title <- tool_named$paths$fov
  if (!nzchar(fov_title) && is.list(tool_flat) && !is.null(tool_flat$LoadIMCSegmented$paths$fov)) {
    fov_title <- tool_flat$LoadIMCSegmented$paths$fov
  }
  if (!nzchar(fov_title) && is.list(tool_flat) && !is.null(tool_flat$paths$fov)) {
    fov_title <- tool_flat$paths$fov
  }
  if (!nzchar(fov_title)) {
    md <- obj[[]]; if ("fov" %in% names(md)) fov_title <- as.character(md$fov[1])
  }
  
  title_suffix <- sprintf("assay=%s | layer=%s | combine=%s | boundary=%s | rot=%d",
                          assay, layer, combine, boundaries, as.integer(spatial_rotate_k) %% 4L)
  
  # composites
  composite_val <- NULL; rgb_cols <- NULL
  if (combine == "rgb") {
    if (length(features) != 3) stop("combine='rgb' requires exactly 3 features (R,G,B).")
    Vrgb <- V
    for (i in 1:3) {
      r <- range(Vrgb[i,], na.rm = TRUE)
      if (is.finite(r[2]) && r[2] > r[1]) {
        Vrgb[i,] <- (Vrgb[i,] - r[1]) / (r[2] - r[1])
        Vrgb[i,] <- pmax(0, pmin(1, Vrgb[i,]))
      } else Vrgb[i,] <- 0
    }
    rgb_cols <- grDevices::rgb(Vrgb[1,], Vrgb[2,], Vrgb[3,], alpha = alpha)
  } else if (combine == "max") {
    composite_val <- apply(V, 2, max, na.rm = TRUE)
  } else if (combine == "mean") {
    composite_val <- colMeans(V, na.rm = TRUE)
  } else if (combine == "sum") {
    composite_val <- colSums(V, na.rm = TRUE)
  }
  
  # draw
  if (boundaries == "segmentation") {
    seg_df <- .get_seg_df_for_obj(
      obj       = obj,
      mask_path = mask_path,
      rebuild   = rebuild_seg_if_missing,
      verbose   = verbose
    )
    if (is.null(seg_df)) stop("No segmentation boundary/coords available. Could not rebuild from mask.")
    
    rot <- .rotate_xy(seg_df$x, seg_df$y, k = spatial_rotate_k, cx = cx0, cy = cy0)
    seg_df$x <- rot$x; seg_df$y <- rot$y
    
    if (combine == "stack") {
      if (is.null(colors)) {
        colors <- if (requireNamespace("scales", quietly = TRUE)) scales::hue_pal()(length(features))
        else grDevices::rainbow(length(features))
      }
      col_map <- setNames(colors, features)
      long_list <- lapply(seq_along(features), function(i) {
        f <- features[i]
        v <- V[match(f, rownames(V)), ]
        data.frame(cell = colnames(obj), ftr = f, val = as.numeric(v), check.names = FALSE)
      })
      val_long <- do.call(rbind, long_list)
      seg_long <- merge(seg_df, val_long, by = "cell", all.x = TRUE, sort = FALSE)
      
      nv <- ave(seg_long$val, seg_long$ftr, FUN = function(x) {
        r <- range(x, na.rm = TRUE)
        if (!is.finite(r[2]) || r[2] <= r[1]) rep(0, length(x)) else (x - r[1]) / (r[2] - r[1])
      })
      ar <- if (is.null(alpha_range) || length(alpha_range) < 2 || any(!is.finite(alpha_range))) c(0, alpha) else alpha_range
      a_min <- ar[1]; a_max <- ar[2]
      if (isTRUE(stack_rescale_alpha_by_n)) a_max <- a_max / max(1, length(features))
      seg_long$alpha_val <- a_min + (a_max - a_min) * (nv ^ stack_alpha_gamma)
      
      p <- p +
        ggplot2::geom_polygon(
          ggplot2::aes(x = x, y = y, group = cell, fill = ftr, alpha = alpha_val),
          data = seg_long, color = NA, linewidth = 0, na.rm = TRUE
        ) +
        ggplot2::scale_fill_manual(values = col_map) +
        ggplot2::scale_alpha_identity(guide = "none") +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | features=%s",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          fill = "Feature"
        )
      
    } else if (combine == "rgb") {
      rgb_df <- data.frame(cell = colnames(obj), col = rgb_cols, stringsAsFactors = FALSE)
      seg_plot <- merge(seg_df, rgb_df, by = "cell", all.x = TRUE, sort = FALSE)
      p <- p +
        ggplot2::geom_polygon(
          ggplot2::aes(x = x, y = y, group = cell),
          data = seg_plot, fill = seg_plot$col, color = NA, linewidth = 0, alpha = 1, na.rm = TRUE
        ) +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | RGB = [%s]",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          fill = NULL
        )
      
    } else {
      val_df  <- data.frame(cell = colnames(obj), val = composite_val, stringsAsFactors = FALSE)
      seg_plot <- merge(seg_df, val_df, by = "cell", all.x = TRUE, sort = FALSE)
      p <- p +
        ggplot2::geom_polygon(
          ggplot2::aes(x = x, y = y, group = cell, fill = val),
          data = seg_plot, linewidth = 0, color = NA, alpha = alpha, na.rm = TRUE
        ) +
        (if (!is.null(limits)) ggplot2::scale_fill_continuous(limits = limits) else ggplot2::scale_fill_continuous()) +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | features=%s",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          fill = "Composite"
        )
    }
    
  } else { # centroids
    md <- obj[[]]
    df_pts <- data.frame(x = md$centroid_x, y = md$centroid_y, row.names = colnames(obj))
    rot <- .rotate_xy(df_pts$x, df_pts$y, k = spatial_rotate_k, cx = cx0, cy = cy0)
    df_pts$x <- rot$x; df_pts$y <- rot$y
    
    if (combine == "stack") {
      if (is.null(colors)) {
        colors <- if (requireNamespace("scales", quietly = TRUE)) scales::hue_pal()(length(features))
        else grDevices::rainbow(length(features))
      }
      long_list <- lapply(seq_along(features), function(i) {
        f <- features[i]; v <- V[match(f, rownames(V)), ]
        data.frame(cell = colnames(obj), ftr = f, val = as.numeric(v), check.names = FALSE)
      })
      val_long <- do.call(rbind, long_list)
      df_long <- merge(
        data.frame(cell = rownames(df_pts), x = df_pts$x, y = df_pts$y, stringsAsFactors = FALSE),
        val_long, by = "cell", all.x = TRUE, sort = FALSE
      )
      nv <- ave(df_long$val, df_long$ftr, FUN = function(x) {
        r <- range(x, na.rm = TRUE)
        if (!is.finite(r[2]) || r[2] <= r[1]) rep(0, length(x)) else (x - r[1]) / (r[2] - r[1])
      })
      ar <- if (is.null(alpha_range) || length(alpha_range) < 2 || any(!is.finite(alpha_range))) c(0, alpha) else alpha_range
      a_min <- ar[1]; a_max <- ar[2]
      if (isTRUE(stack_rescale_alpha_by_n)) a_max <- a_max / max(1, length(features))
      df_long$alpha_val <- a_min + (a_max - a_min) * (nv ^ stack_alpha_gamma)
      
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = x, y = y, color = ftr, alpha = alpha_val),
          data = df_long, size = pt_size, na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_alpha_identity(guide = "none") +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | features=%s",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          color = "Feature"
        )
      
    } else if (combine == "rgb") {
      df_pts$col <- rgb_cols
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = x, y = y),
          data = df_pts, size = pt_size, color = df_pts$col, alpha = 1, na.rm = TRUE
        ) +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | RGB = [%s]",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          color = NULL
        )
      
    } else {
      df_pts$val <- composite_val
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = x, y = y, color = val),
          data = df_pts, size = pt_size, alpha = alpha, na.rm = TRUE
        ) +
        (if (!is.null(limits)) ggplot2::scale_color_continuous(limits = limits) else ggplot2::scale_color_continuous()) +
        ggplot2::labs(
          title = sprintf("FOV: %s | %s | features=%s",
                          fov_title, title_suffix, paste(features, collapse = ", ")),
          color = "Composite"
        )
    }
  }
  
  if (!is.null(export_path)) {
    ggplot2::ggsave(export_path, p, width = W / dpi, height = H / dpi, dpi = dpi)
    if (isTRUE(verbose)) message("Saved overlay: ", normalizePath(export_path, winslash = "/"))
    return(invisible(p))
  }
  p
}
