
## Authors : Hanen Daayeb and Frederic Ouimet

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(LaplacesDemon)
  library(misc3d)
  library(parallel)
  library(patchwork)
  library(plot3D)
  library(plotly)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Paths (edit here)
# ============================================================

BASE_PATH <- "C:\\Users\\oufr6443\\Dropbox\\Dossier\\Prof\\Latex\\In_preparation\\Dirichlet_KDE_missing_data"

TABLE_DIR <- file.path(BASE_PATH, "tables")
FIG_DIR   <- file.path(BASE_PATH, "figures")

dir_create_safe <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
dir_create_safe(TABLE_DIR)
dir_create_safe(FIG_DIR)

# ============================================================
# Hyperparameters (edit here)
# ============================================================

SEED <- 1
EPS  <- 1e-12

# Global simulation mechanism parameters
RHO    <- 0.7
BETA_1 <- 1

# --- Simplex grids (3-part simplex -> 2D plots) ---
GRID_RES_VIS   <- 300
GRID_RES_TABLE <- 300
GRID_RES_COMP  <- 300

EPS_GRID_VIS   <- 0.01
EPS_GRID_TABLE <- 0.01
EPS_GRID_COMP  <- 0.01

# --- LSCV integration grid (deterministic interior grid) ---
LSCV_GRID_RES <- 40
LSCV_EPS_GRID <- 0.01

# --- Section toggles ---
RUN_SECTION_1_VISUALS     <- TRUE   # MAR + true/estimated densities (Model 1 & 2)
RUN_SECTION_2_TABLES      <- TRUE   # Tables 1–2 + Figure 4–6
RUN_SECTION_3_COMPARISON  <- TRUE   # Figure 7–9 (Dirichlet vs ALR vs ILR), Model1 vs Model2, feasible only
RUN_SECTION_4_APPLICATION <- TRUE   # Figure 10 (NHANES application)

# --- Section 1: visualization settings ---
N_VIS <- 2000
P_VIS <- 0.9  # observation probability

# Candidate Dirichlet bandwidth grid for Figures 2–3 (LSCV)
B_GRID_CONTOUR <- seq(0.05, 0.35, length.out = 7)

# --- Section 2: Monte Carlo tables settings ---
N_VALUES_TABLE <- c(100, 200, 400, 800)
P_VALUES_TABLE <- c(0.95, 0.90, 0.80, 0.60)
MC_REPS_TABLE  <- 1000

# Candidate Dirichlet bandwidth grid for Tables 1–2 (LSCV)
B_GRID_TABLE <- seq(0.01, 0.35, length.out = 35)

# Sample size variation (Figure 7)
N_VALUES_COMP_N  <- c(100, 200, 400, 800)
P_COMP_N_FIXED   <- 0.8
MC_REPS_COMP_N   <- 1000

# Missingness variation (Figure 8)
N_COMP_MISS_FIXED  <- 400
MISSING_LEVELS     <- c(0.05, 0.10, 0.20, 0.40)
MC_REPS_COMP_MISS  <- 1000

# Joint variation (Figure 9)
N_VALUES_COMP_GLOBAL      <- c(100, 200, 400, 800)
MISSING_LEVELS_GLOBAL     <- c(0.10, 0.20, 0.40)
MC_REPS_COMP_GLOBAL       <- 1000

# Bandwidth candidate grids for comparisons (LSCV)
B_GRID_COMP   <- seq(0.01, 0.35, length.out = 35)
H_GRID_FACTOR <- seq(0.05, 1.0, length.out = 35)

# --- Section 4: Application (NHANES) ---
NHANES_CYCLE_SUFFIX <- "J"  # 2017–2018
B_GRID_APP          <- seq(0.01, 0.35, length.out = 35)
GRID_RES_APP        <- 300
EPS_GRID_APP        <- 0.01

# ============================================================
# Figure label/text sizing toggles (EDIT HERE)
#   Organized by Figure number, as requested.
#   Each figure has its own independent set of font-size toggles.
# ============================================================

# ---------------------------
# Figure 1 — MAR visualization (plotly)
# ---------------------------
FIG1_PLOTLY_AXIS_TITLE_SIZE <- 26
FIG1_PLOTLY_TICK_SIZE       <- 24
FIG1_PLOTLY_LEGEND_SIZE     <- 26
FIG1_PLOTLY_TITLE_SIZE      <- 26   # (not used unless you add a title)
FIG1_PLOTLY_TEXT_SIZE       <- 26   # (for add_text, if used)

# Figure 1 legend symbol sizing (Option B: legend-only big markers)
FIG1_MARKER_SIZE_OBS        <- 6
FIG1_MARKER_SIZE_MISS       <- 6
FIG1_MARKER_SIZE_XSTAR      <- 6
FIG1_MARKER_SIZE_X          <- 6
FIG1_LEGEND_MARKER_MULT     <- 2.5      # double the legend marker size
FIG1_LEGEND_DUMMY_X         <- -999   # off-range dummy point
FIG1_LEGEND_DUMMY_Y         <- -999   # off-range dummy point

# ---------------------------
# Figure 2 — True density (ggplot contour)
# ---------------------------
FIG2_GG_BASE_SIZE         <- 16
FIG2_AXIS_TITLE_SIZE      <- 16
FIG2_AXIS_TEXT_SIZE       <- 14
FIG2_LEGEND_TITLE_SIZE    <- 14
FIG2_LEGEND_TEXT_SIZE     <- 14
FIG2_STRIP_TEXT_SIZE      <- 14
FIG2_PLOT_TITLE_SIZE      <- 14

# ---------------------------
# Figure 3 — Estimated density (ggplot contour)
# ---------------------------
FIG3_GG_BASE_SIZE         <- 16
FIG3_AXIS_TITLE_SIZE      <- 16
FIG3_AXIS_TEXT_SIZE       <- 14
FIG3_LEGEND_TITLE_SIZE    <- 14
FIG3_LEGEND_TEXT_SIZE     <- 14
FIG3_STRIP_TEXT_SIZE      <- 14
FIG3_PLOT_TITLE_SIZE      <- 14

# ---------------------------
# Figure 4 — Mean ISE plot (ggplot line + points)
# ---------------------------
FIG4_GG_BASE_SIZE         <- 16
FIG4_AXIS_TITLE_SIZE      <- 16
FIG4_AXIS_TEXT_SIZE       <- 14
FIG4_LEGEND_TITLE_SIZE    <- 14
FIG4_LEGEND_TEXT_SIZE     <- 14
FIG4_STRIP_TEXT_SIZE      <- 14
FIG4_PLOT_TITLE_SIZE      <- 14

# ---------------------------
# Figure 5 — ISE metrics (Model 1) (ggplot patchwork)
# ---------------------------
FIG5_GG_BASE_SIZE         <- 14
FIG5_AXIS_TITLE_SIZE      <- 14
FIG5_AXIS_TEXT_SIZE       <- 12
FIG5_LEGEND_TITLE_SIZE    <- 12
FIG5_LEGEND_TEXT_SIZE     <- 12
FIG5_STRIP_TEXT_SIZE      <- 12
FIG5_PLOT_TITLE_SIZE      <- 12

# ---------------------------
# Figure 6 — ISE metrics (Model 2) (ggplot patchwork)
# ---------------------------
FIG6_GG_BASE_SIZE         <- 14
FIG6_AXIS_TITLE_SIZE      <- 14
FIG6_AXIS_TEXT_SIZE       <- 12
FIG6_LEGEND_TITLE_SIZE    <- 12
FIG6_LEGEND_TEXT_SIZE     <- 12
FIG6_STRIP_TEXT_SIZE      <- 12
FIG6_PLOT_TITLE_SIZE      <- 12

# ---------------------------
# Figure 7 — Comparison: sample size variation (ggplot boxplot)
# ---------------------------
FIG7_GG_BASE_SIZE         <- 18
FIG7_AXIS_TITLE_SIZE      <- 18
FIG7_AXIS_TEXT_SIZE       <- 16
FIG7_LEGEND_TITLE_SIZE    <- 16   # (legend hidden in code, but kept as a toggle)
FIG7_LEGEND_TEXT_SIZE     <- 16
FIG7_STRIP_TEXT_SIZE      <- 16
FIG7_PLOT_TITLE_SIZE      <- 16

# ---------------------------
# Figure 8 — Comparison: missingness variation (ggplot boxplot)
# ---------------------------
FIG8_GG_BASE_SIZE         <- 18
FIG8_AXIS_TITLE_SIZE      <- 18
FIG8_AXIS_TEXT_SIZE       <- 16
FIG8_LEGEND_TITLE_SIZE    <- 14   # (legend hidden in code, but kept as a toggle)
FIG8_LEGEND_TEXT_SIZE     <- 14
FIG8_STRIP_TEXT_SIZE      <- 16
FIG8_PLOT_TITLE_SIZE      <- 16

# ---------------------------
# Figure 9 — Comparison: joint effect (n × missing) (ggplot facets)
# ---------------------------
FIG9_GG_BASE_SIZE         <- 22
FIG9_AXIS_TITLE_SIZE      <- 22
FIG9_AXIS_TEXT_SIZE       <- 20
FIG9_AXIS_TEXT_X_SIZE     <- 20  # special x tick label size (rotated)
FIG9_LEGEND_TITLE_SIZE    <- 18    # (legend hidden in code, but kept as a toggle)
FIG9_LEGEND_TEXT_SIZE     <- 18
FIG9_STRIP_TEXT_SIZE      <- 20
FIG9_PLOT_TITLE_SIZE      <- 20

# ---------------------------
# Figure 10a — NHANES 2D raster (ggplot) (NO TITLE)
# ---------------------------
FIG10A_GG_BASE_SIZE        <- 18
FIG10A_AXIS_TITLE_SIZE     <- 20
FIG10A_AXIS_TEXT_SIZE      <- 16
FIG10A_LEGEND_TITLE_SIZE   <- 16
FIG10A_LEGEND_TEXT_SIZE    <- 16
FIG10A_STRIP_TEXT_SIZE     <- 16  # (not used here, but kept as a toggle)
FIG10A_PLOT_TITLE_SIZE     <- 16  # (NO TITLE, but kept as a toggle)

# ---------------------------
# Figure 10b — NHANES 3D static plot (plot3D) (NO TITLE) — applies to all 4 views
# ---------------------------
FIG10B_CEX_LAB            <- 1.3
FIG10B_CEX_AXIS           <- 1.2
FIG10B_COLKEY_CEX_AXIS    <- 1.2
FIG10B_COLKEY_CEX_CLAB    <- 1.2
FIG10B_MODE_LABEL_CEX     <- 1.2    # text size for "Mode" label (text3D)

# ---------------------------
# Plotly export timing (not text size; kept as a tweakable control)
# ---------------------------
PLOTLY_EXPORT_DELAY_DEFAULT <- 0.2
PLOTLY_EXPORT_DELAY_WEBGL   <- 1.5

set.seed(SEED)

# ============================================================
# Timers (Sys.time) helpers
# ============================================================

timer_report <- function(label, start_time, end_time = Sys.time()) {
  elapsed_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  elapsed_min <- elapsed_sec / 60
  message(sprintf("[TIMER] %s: %.2f seconds (%.2f minutes)", label, elapsed_sec, elapsed_min))
  invisible(list(seconds = elapsed_sec, minutes = elapsed_min))
}

# ============================================================
# Theme helpers (GGPlot)
#   One helper, then per-figure wrappers (so each figure uses its own toggles).
# ============================================================

theme_gg_by_sizes <- function(base_size,
                              axis_title_size,
                              axis_text_size,
                              legend_title_size,
                              legend_text_size,
                              strip_text_size,
                              plot_title_size) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.title   = element_text(size = axis_title_size),
      axis.text    = element_text(size = axis_text_size),
      legend.title = element_text(size = legend_title_size),
      legend.text  = element_text(size = legend_text_size),
      strip.text   = element_text(size = strip_text_size),
      plot.title   = element_text(size = plot_title_size)
    )
}

theme_fig2_gg <- function() theme_gg_by_sizes(
  FIG2_GG_BASE_SIZE, FIG2_AXIS_TITLE_SIZE, FIG2_AXIS_TEXT_SIZE,
  FIG2_LEGEND_TITLE_SIZE, FIG2_LEGEND_TEXT_SIZE,
  FIG2_STRIP_TEXT_SIZE, FIG2_PLOT_TITLE_SIZE
)

theme_fig3_gg <- function() theme_gg_by_sizes(
  FIG3_GG_BASE_SIZE, FIG3_AXIS_TITLE_SIZE, FIG3_AXIS_TEXT_SIZE,
  FIG3_LEGEND_TITLE_SIZE, FIG3_LEGEND_TEXT_SIZE,
  FIG3_STRIP_TEXT_SIZE, FIG3_PLOT_TITLE_SIZE
)

theme_fig4_gg <- function() theme_gg_by_sizes(
  FIG4_GG_BASE_SIZE, FIG4_AXIS_TITLE_SIZE, FIG4_AXIS_TEXT_SIZE,
  FIG4_LEGEND_TITLE_SIZE, FIG4_LEGEND_TEXT_SIZE,
  FIG4_STRIP_TEXT_SIZE, FIG4_PLOT_TITLE_SIZE
)

theme_fig5_gg <- function() theme_gg_by_sizes(
  FIG5_GG_BASE_SIZE, FIG5_AXIS_TITLE_SIZE, FIG5_AXIS_TEXT_SIZE,
  FIG5_LEGEND_TITLE_SIZE, FIG5_LEGEND_TEXT_SIZE,
  FIG5_STRIP_TEXT_SIZE, FIG5_PLOT_TITLE_SIZE
)

theme_fig6_gg <- function() theme_gg_by_sizes(
  FIG6_GG_BASE_SIZE, FIG6_AXIS_TITLE_SIZE, FIG6_AXIS_TEXT_SIZE,
  FIG6_LEGEND_TITLE_SIZE, FIG6_LEGEND_TEXT_SIZE,
  FIG6_STRIP_TEXT_SIZE, FIG6_PLOT_TITLE_SIZE
)

theme_fig7_gg <- function() theme_gg_by_sizes(
  FIG7_GG_BASE_SIZE, FIG7_AXIS_TITLE_SIZE, FIG7_AXIS_TEXT_SIZE,
  FIG7_LEGEND_TITLE_SIZE, FIG7_LEGEND_TEXT_SIZE,
  FIG7_STRIP_TEXT_SIZE, FIG7_PLOT_TITLE_SIZE
)

theme_fig8_gg <- function() theme_gg_by_sizes(
  FIG8_GG_BASE_SIZE, FIG8_AXIS_TITLE_SIZE, FIG8_AXIS_TEXT_SIZE,
  FIG8_LEGEND_TITLE_SIZE, FIG8_LEGEND_TEXT_SIZE,
  FIG8_STRIP_TEXT_SIZE, FIG8_PLOT_TITLE_SIZE
)

theme_fig9_gg <- function() theme_gg_by_sizes(
  FIG9_GG_BASE_SIZE, FIG9_AXIS_TITLE_SIZE, FIG9_AXIS_TEXT_SIZE,
  FIG9_LEGEND_TITLE_SIZE, FIG9_LEGEND_TEXT_SIZE,
  FIG9_STRIP_TEXT_SIZE, FIG9_PLOT_TITLE_SIZE
)

theme_fig10a_gg <- function() theme_gg_by_sizes(
  FIG10A_GG_BASE_SIZE, FIG10A_AXIS_TITLE_SIZE, FIG10A_AXIS_TEXT_SIZE,
  FIG10A_LEGEND_TITLE_SIZE, FIG10A_LEGEND_TEXT_SIZE,
  FIG10A_STRIP_TEXT_SIZE, FIG10A_PLOT_TITLE_SIZE
)

# ============================================================
# Parallelization on cores (local laptop, cores = detectCores()-1)
# ============================================================

cores_detected <- parallel::detectCores()
if (is.na(cores_detected) || cores_detected < 2) cores_detected <- 2
NUM_CORES <- max(1, cores_detected - 1)

libraries_to_load <- c(
  "LaplacesDemon"
)

# Objects/functions that MC workers must have (exported once at cluster creation)
vars_to_export <- c(
  # constants
  "RHO", "BETA_1", "EPS",
  "B_GRID_TABLE", "B_GRID_COMP", "H_GRID_FACTOR", "LSCV_GRID",
  
  # core utilities
  "normalize_comp",
  "sample_mixture_dirichlet",
  "generate_data",
  "estimate_pi_hat",
  "true_density_fun",
  
  # Dirichlet KDE + LSCV
  "dirichlet_kernel_matrix",
  "dirichlet_kde_ipw",
  "lscv_dirichlet",
  
  # log-ratio KDE + LSCV (comparison study)
  "alr_transform",
  "ilr_transform",
  "jacobian_dzdy_alr",
  "jacobian_dzdy_ilr",
  "gaussian_kernel_matrix",
  "logratio_kde_ipw",
  "h_reference",
  "lscv_logratio",
  
  # replication-level wrappers
  "run_single_replication_table",
  "comparison_one_rep"
)

setup_parallel_cluster <- function(num_cores = NUM_CORES) {
  cl <- parallel::makeCluster(num_cores)
  
  # Export library list then load packages on each worker
  parallel::clusterExport(cl, varlist = "libraries_to_load", envir = .GlobalEnv)
  invisible(parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      lapply(libraries_to_load, library, character.only = TRUE)
    })
  }))
  
  # Export required functions/objects
  parallel::clusterExport(cl, varlist = vars_to_export, envir = .GlobalEnv)
  
  cl
}

# ============================================================
# Models (Dirichlet mixtures)
# ============================================================

MODELS <- list(
  Model1 = list(
    w1 = 0.4, w2 = 0.6,
    alpha1 = c(1.3, 1.6, 1.0),
    alpha2 = c(1.7, 1.2, 2.5)
  ),
  Model2 = list(
    w1 = 0.4, w2 = 0.6,
    alpha1 = c(4, 1, 2),
    alpha2 = c(1, 3, 2)
  )
)

# ============================================================
# Helper utilities
# ============================================================

normalize_comp <- function(Y, eps = EPS) {
  Y <- as.matrix(Y)
  Y <- pmax(Y, eps)
  Y / rowSums(Y)
}

generate_simplex_grid <- function(res = 80, eps_grid = 0.01) {
  H <- res - 1
  i_seq <- 1:(H - 2)
  j_count <- H - i_seq - 1
  
  i <- rep(i_seq, j_count)
  j <- sequence(j_count)
  k <- H - i - j
  
  scale <- 1 - 3 * eps_grid
  
  data.frame(
    s1 = eps_grid + scale * (i / H),
    s2 = eps_grid + scale * (j / H),
    s3 = eps_grid + scale * (k / H)
  )
}

LSCV_GRID <- generate_simplex_grid(res = LSCV_GRID_RES, eps_grid = LSCV_EPS_GRID)

true_density_fun <- function(grid, pars, eps = EPS) {
  S <- as.matrix(grid[, c("s1", "s2", "s3")])
  S <- normalize_comp(S, eps = eps)
  
  apply(S, 1, function(p) {
    pars$w1 * ddirichlet(p, pars$alpha1) +
      pars$w2 * ddirichlet(p, pars$alpha2)
  })
}

sample_mixture_dirichlet <- function(n, w1, alpha1, alpha2) {
  comps <- rbinom(n, 1, w1)
  t(vapply(comps, function(cmp) {
    if (cmp == 1) rdirichlet(1, alpha1) else rdirichlet(1, alpha2)
  }, FUN.VALUE = numeric(length(alpha1))))
}

generate_data <- function(n, pars, rho, beta_1, p_obs) {
  y <- sample_mixture_dirichlet(n, pars$w1, pars$alpha1, pars$alpha2)
  y <- normalize_comp(y)
  
  x_star <- matrix(rnorm(n * 3), n, 3)
  x <- rho * y + sqrt(1 - rho^2) * x_star
  
  beta_0 <- log(p_obs / (1 - p_obs)) - beta_1 * mean(x[, 1])
  pi <- 1 / (1 + exp(-(beta_0 + beta_1 * x[, 1])))
  
  delta <- rbinom(n, 1, pi)
  list(x = x, y = y, x_star = x_star, delta = delta)
}

estimate_pi_hat <- function(x, delta, eps = EPS) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  sigma <- mean(apply(x, 2, sd))
  h <- 1.06 * sigma * n^(-1 / (p + 4))
  
  D2 <- as.matrix(dist(x / h))^2
  W  <- exp(-0.5 * D2) # constants cancel in NW ratio
  
  numer <- as.numeric(crossprod(delta, W))
  denom <- colSums(W)
  
  pi_hat <- numer / pmax(denom, eps)
  pmin(pmax(pi_hat, eps), 1 - eps)
}

# Dirichlet kernel matrix K[j,i] = kappa_{S_i,b}(X_j)
dirichlet_kernel_matrix <- function(X, S, b, eps = EPS) {
  X <- normalize_comp(X, eps = eps)
  S <- normalize_comp(S, eps = eps)
  
  alpha <- S / b + 1
  logC  <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha))
  
  # sum_k (S_{i,k}/b) log(X_{j,k})
  Klog <- log(X) %*% t(S / b)
  
  exp(Klog + matrix(rep(logC, each = nrow(X)), nrow = nrow(X)))
}

# IPW Dirichlet KDE (HT scaling by 1/n_total)
dirichlet_kde_ipw <- function(y, delta, pi, b, eval_points, eps = EPS) {
  n_total <- nrow(y)
  idx <- which(delta == 1)
  if (length(idx) == 0) return(rep(NA_real_, nrow(eval_points)))
  
  y_obs <- y[idx, , drop = FALSE]
  w_obs <- 1 / pmax(pi[idx], eps)
  
  S <- as.matrix(eval_points[, c("s1", "s2", "s3")])
  K <- dirichlet_kernel_matrix(y_obs, S, b, eps = eps)
  
  as.numeric(crossprod(w_obs, K)) / n_total
}

# ============================================================
# Log-ratio transformations + Jacobians (checked)
# ============================================================

alr_transform <- function(Y, eps = EPS) {
  Y <- normalize_comp(Y, eps = eps)
  cbind(
    log(Y[, 1] / Y[, 3]),
    log(Y[, 2] / Y[, 3])
  )
}

ilr_transform <- function(Y, eps = EPS) {
  Y <- normalize_comp(Y, eps = eps)
  z1 <- sqrt(1 / 2) * log(Y[, 1] / Y[, 2])
  z2 <- sqrt(1 / 6) * log((Y[, 1] * Y[, 2]) / (Y[, 3]^2))
  cbind(z1, z2)
}

# |det(dz/d(y1,y2))|, with y3 = 1 - y1 - y2
jacobian_dzdy_alr <- function(S, eps = EPS) {
  S <- normalize_comp(S, eps = eps)
  1 / (S[, 1] * S[, 2] * S[, 3])
}

jacobian_dzdy_ilr <- function(S, eps = EPS) {
  S <- normalize_comp(S, eps = eps)
  1 / (sqrt(3) * S[, 1] * S[, 2] * S[, 3])
}

gaussian_kernel_matrix <- function(ZX, ZS, h) {
  ZX <- as.matrix(ZX)
  ZS <- as.matrix(ZS)
  d  <- ncol(ZX)
  
  n <- nrow(ZX)
  m <- nrow(ZS)
  
  d2 <- matrix(rowSums(ZX^2), n, m) +
    matrix(rowSums(ZS^2), n, m, byrow = TRUE) -
    2 * ZX %*% t(ZS)
  
  exp(-0.5 * d2 / (h^2)) / ((2 * pi)^(d / 2) * h^d)
}

logratio_kde_ipw <- function(Y, delta, pi, h, eval_points,
                             transform_fun, jacobian_fun, eps = EPS) {
  n_total <- nrow(Y)
  idx <- which(delta == 1)
  if (length(idx) == 0) return(rep(NA_real_, nrow(eval_points)))
  
  Y_obs <- Y[idx, , drop = FALSE]
  w_obs <- 1 / pmax(pi[idx], eps)
  
  ZY <- transform_fun(Y_obs, eps = eps)
  ZS <- transform_fun(eval_points, eps = eps)
  
  Kz <- gaussian_kernel_matrix(ZY, ZS, h)
  fz <- as.numeric(crossprod(w_obs, Kz)) / n_total
  
  f_simplex <- fz * jacobian_fun(eval_points, eps = eps)
  pmax(f_simplex, 0)
}

h_reference <- function(Y, transform_fun, eps = EPS) {
  Z <- transform_fun(Y, eps = eps)
  d <- ncol(Z)
  n <- nrow(Z)
  mean(apply(Z, 2, sd)) * n^(-1 / (d + 4))
}

# ============================================================
# LSCV bandwidth selection (IPW-adapted)
# ============================================================

lscv_dirichlet <- function(y, delta, pi, b, lscv_grid = LSCV_GRID, eps = EPS) {
  n_total <- nrow(y)
  idx <- which(delta == 1)
  if (length(idx) < 2) return(Inf)
  
  y_obs <- y[idx, , drop = FALSE]
  w_obs <- 1 / pmax(pi[idx], eps)
  
  # Term 1: approximate ∫ f_hat^2
  f_grid <- dirichlet_kde_ipw(y, delta, pi, b, lscv_grid, eps = eps)
  term1  <- mean(f_grid^2, na.rm = TRUE) * 0.5 # 0.5 is the area of S_2
  
  # Term 2: weighted leave-one-out at observed points
  K_yy   <- dirichlet_kernel_matrix(y_obs, y_obs, b, eps = eps)  # rows=j, cols=i
  sum_wK <- as.numeric(crossprod(w_obs, K_yy))                   # Σ_j w_j K_{ji}
  
  loo <- (sum_wK - w_obs * diag(K_yy)) / (n_total - 1)
  term2 <- (2 / n_total) * sum(w_obs * loo)
  
  term1 - term2
}

lscv_logratio <- function(Y, delta, pi, h, lscv_grid,
                          transform_fun, jacobian_fun, eps = EPS) {
  n_total <- nrow(Y)
  idx <- which(delta == 1)
  if (length(idx) < 2) return(Inf)
  
  Y_obs <- Y[idx, , drop = FALSE]
  w_obs <- 1 / pmax(pi[idx], eps)
  
  ZY <- transform_fun(Y_obs, eps = eps)
  
  # Term 1: approximate ∫ f_hat^2 using simplex grid
  f_grid <- logratio_kde_ipw(
    Y = Y, delta = delta, pi = pi, h = h,
    eval_points = lscv_grid,
    transform_fun = transform_fun,
    jacobian_fun  = jacobian_fun,
    eps = eps
  )
  term1 <- mean(f_grid^2, na.rm = TRUE) * 0.5 # 0.5 is the area of S_2
  
  # Term 2: leave-one-out at observed points
  K_zz   <- gaussian_kernel_matrix(ZY, ZY, h)
  sum_wK <- as.numeric(crossprod(w_obs, K_zz))
  loo_z  <- (sum_wK - w_obs * diag(K_zz)) / (n_total - 1)
  loo_s  <- loo_z * jacobian_fun(Y_obs, eps = eps)
  
  term2 <- (2 / n_total) * sum(w_obs * loo_s)
  term1 - term2
}

# ============================================================
# Output helpers (PDF figures + LaTeX tables)
# ============================================================

save_gg_pdf <- function(plot_obj, filename_pdf, width = 7, height = 5, dpi = 300) {
  if (!grepl("\\.pdf$", filename_pdf, ignore.case = TRUE)) {
    filename_pdf <- paste0(filename_pdf, ".pdf")
  }
  ggsave(filename = filename_pdf, plot = plot_obj, width = width, height = height, dpi = dpi)
}

save_plotly_pdf <- function(p, file_pdf, width = 900, height = 700,
                            prefer_webshot = FALSE,
                            webshot_delay = PLOTLY_EXPORT_DELAY_DEFAULT) {
  if (!grepl("\\.pdf$", file_pdf, ignore.case = TRUE)) {
    file_pdf <- paste0(file_pdf, ".pdf")
  }
  
  ok <- FALSE
  
  # 1) Try native export first (kaleido/orca), unless we force webshot (recommended for WebGL 3D)
  if (!prefer_webshot) {
    try({
      plotly::save_image(p, file_pdf, width = width, height = height)
      ok <- file.exists(file_pdf)
    }, silent = TRUE)
  }
  
  # 2) Fallback: save widget -> webshot2/pagedown/webshot to PDF
  if (!ok) {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("Failed to save plotly figure to PDF: please install 'htmlwidgets' and 'webshot2' (recommended), or configure plotly image export (kaleido/orca).")
    }
    
    tmp_html <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget(p, file = tmp_html, selfcontained = TRUE)
    
    if (requireNamespace("webshot2", quietly = TRUE)) {
      try({
        webshot2::webshot(
          url = tmp_html,
          file = file_pdf,
          vwidth = width,
          vheight = height,
          delay = webshot_delay
        )
        ok <- file.exists(file_pdf)
      }, silent = TRUE)
    }
    
    if (!ok && requireNamespace("pagedown", quietly = TRUE)) {
      try({
        pagedown::chrome_print(input = tmp_html, output = file_pdf)
        ok <- file.exists(file_pdf)
      }, silent = TRUE)
    }
    
    if (!ok && requireNamespace("webshot", quietly = TRUE)) {
      try({
        webshot::webshot(
          url = tmp_html,
          file = file_pdf,
          vwidth = width,
          vheight = height,
          delay = webshot_delay
        )
        ok <- file.exists(file_pdf)
      }, silent = TRUE)
    }
  }
  
  if (!ok) {
    stop("Failed to save plotly figure to PDF. Install one of: 'webshot2' (recommended), 'pagedown', or 'webshot', or configure plotly image export (kaleido/orca).")
  }
  
  invisible(ok)
}

write_latex_table <- function(summary_df, file_tex, caption, label, table_comment) {
  fmt_f <- function(x, digits) formatC(x, format = "f", digits = digits)
  
  header <- c("n", "Missing Rate (\\%)", "Mean ISE", "Median ISE", "SD ISE", "IQR ISE", "Mean $b^{\\star}$")
  
  lines <- c(
    paste0("%", table_comment),
    "\\begin{table}[H]",
    "\\renewcommand{\\arraystretch}{1.1}",
    "\\centering",
    "\\begin{tabular}{rrrrrrr}",
    "\\hline",
    paste0(paste(header, collapse = " & "), " \\\\"),
    "\\hline"
  )
  
  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    vals <- c(
      as.integer(r$n),
      fmt_f((1 - r$p) * 100, 0),
      fmt_f(r$mean_ISE, 4),
      fmt_f(r$median_ISE, 4),
      fmt_f(r$SD_ISE, 4),
      fmt_f(r$IQR_ISE, 4),
      fmt_f(r$mean_b, 4)
    )
    lines <- c(lines, paste0(paste(vals, collapse = " & "), " \\\\"))
  }
  
  lines <- c(
    lines,
    "\\hline",
    "\\end{tabular}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\end{table}",
    ""
  )
  
  writeLines(lines, con = file_tex)
}

triangle_outline_df <- function() {
  data.frame(s1 = c(0, 1, 0, 0), s2 = c(0, 0, 1, 0))
}

# ============================================================
# SECTION 1 — MAR + true/estimated densities (Model 1 & 2)
# ============================================================

make_mar_plot <- function(dat, title_text) {
  y     <- dat$y
  x     <- dat$x
  xstar <- dat$x_star
  delta <- dat$delta
  
  idx_obs  <- which(delta == 1)
  idx_miss <- which(delta == 0)
  
  triangle_2d <- data.frame(
    x = c(0, 1, 0, 0),
    y = c(0, 0, 1, 0)
  )
  
  # Legend-only big markers: keep real markers small, add dummy off-range traces for legend
  fig_2d <- plot_ly() %>%
    
    # Simplex outline (keep in legend as-is)
    add_trace(
      type = "scatter", mode = "lines",
      x = triangle_2d$x, y = triangle_2d$y,
      line = list(color = "black", width = 2),
      name = "Simplex",
      showlegend = TRUE
    ) %>%
    
    # -----------------------
  # Observed Y: real trace (no legend)
  add_trace(
    type = "scatter", mode = "markers",
    x = y[idx_obs, 1], y = y[idx_obs, 2],
    marker = list(size = FIG1_MARKER_SIZE_OBS, color = "#1f77b4"),
    name = "Observed Y",
    showlegend = FALSE,
    legendgroup = "obsY"
  ) %>%
    # Observed Y: legend-only big marker (dummy off-range point)
    add_trace(
      type = "scatter", mode = "markers",
      x = FIG1_LEGEND_DUMMY_X, y = FIG1_LEGEND_DUMMY_Y,
      marker = list(size = FIG1_MARKER_SIZE_OBS * FIG1_LEGEND_MARKER_MULT, color = "#1f77b4"),
      name = "Observed Y",
      showlegend = TRUE,
      legendgroup = "obsY",
      hoverinfo = "skip"
    ) %>%
    
    # -----------------------
  # Missing Y: real trace (no legend)
  add_trace(
    type = "scatter", mode = "markers",
    x = y[idx_miss, 1], y = y[idx_miss, 2],
    marker = list(size = FIG1_MARKER_SIZE_MISS, color = "red", symbol = "x"),
    name = "Missing Y",
    showlegend = FALSE,
    legendgroup = "missY"
  ) %>%
    # Missing Y: legend-only big marker
    add_trace(
      type = "scatter", mode = "markers",
      x = FIG1_LEGEND_DUMMY_X, y = FIG1_LEGEND_DUMMY_Y,
      marker = list(
        size = FIG1_MARKER_SIZE_MISS * FIG1_LEGEND_MARKER_MULT,
        color = "red",
        symbol = "x"
      ),
      name = "Missing Y",
      showlegend = TRUE,
      legendgroup = "missY",
      hoverinfo = "skip"
    ) %>%
    
    # -----------------------
  # Noise X*: real trace (no legend)
  add_trace(
    type = "scatter", mode = "markers",
    x = xstar[, 1], y = xstar[, 2],
    marker = list(size = FIG1_MARKER_SIZE_XSTAR, color = "#2ca02c", symbol = "cross"),
    name = "Noise X*",
    showlegend = FALSE,
    legendgroup = "xstar"
  ) %>%
    # Noise X*: legend-only big marker
    add_trace(
      type = "scatter", mode = "markers",
      x = FIG1_LEGEND_DUMMY_X, y = FIG1_LEGEND_DUMMY_Y,
      marker = list(
        size = FIG1_MARKER_SIZE_XSTAR * FIG1_LEGEND_MARKER_MULT,
        color = "#2ca02c",
        symbol = "cross"
      ),
      name = "Noise X*",
      showlegend = TRUE,
      legendgroup = "xstar",
      hoverinfo = "skip"
    ) %>%
    
    # -----------------------
  # Raw X: real trace (no legend)
  add_trace(
    type = "scatter", mode = "markers",
    x = x[, 1], y = x[, 2],
    marker = list(size = FIG1_MARKER_SIZE_X, color = "#ff7f0e", symbol = "diamond"),
    name = "Raw X",
    showlegend = FALSE,
    legendgroup = "rawX"
  ) %>%
    # Raw X: legend-only big marker
    add_trace(
      type = "scatter", mode = "markers",
      x = FIG1_LEGEND_DUMMY_X, y = FIG1_LEGEND_DUMMY_Y,
      marker = list(
        size = FIG1_MARKER_SIZE_X * FIG1_LEGEND_MARKER_MULT,
        color = "#ff7f0e",
        symbol = "diamond"
      ),
      name = "Raw X",
      showlegend = TRUE,
      legendgroup = "rawX",
      hoverinfo = "skip"
    ) %>%
    
    layout(
      margin = list(t = 40),
      xaxis = list(
        range = c(-0.5, 1.5),
        constrain = "domain",
        title = list(text = "s1", font = list(size = FIG1_PLOTLY_AXIS_TITLE_SIZE)),
        tickfont = list(size = FIG1_PLOTLY_TICK_SIZE)
      ),
      yaxis = list(
        range = c(-0.3, 1.3),
        constrain = "domain",
        title = list(text = "s2", font = list(size = FIG1_PLOTLY_AXIS_TITLE_SIZE)),
        tickfont = list(size = FIG1_PLOTLY_TICK_SIZE),
        scaleanchor = "x"
      ),
      legend = list(
        font = list(size = FIG1_PLOTLY_LEGEND_SIZE),
        itemsizing = "trace",
        groupclick = "togglegroup"
      )
    )
  
  fig_2d
}

# Plasma palette for contour plots (true vs estimated)
make_contour_plot <- function(grid, z, title_text, breaks = NULL, theme_fun) {
  df <- grid %>% mutate(z = z)
  
  p <- ggplot(df, aes(x = s1, y = s2, z = z)) +
    coord_fixed() +
    labs(x = "s1", y = "s2") +
    theme_fun() +
    theme(legend.position = "right")
  
  if (is.null(breaks)) {
    p <- p + geom_contour_filled(bins = 10, na.rm = TRUE)
  } else {
    p <- p + geom_contour_filled(breaks = breaks, na.rm = TRUE)
  }
  
  # Plasma color palette
  if (requireNamespace("viridis", quietly = TRUE)) {
    p <- p + viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", drop = FALSE)
  } else {
    p <- p + scale_fill_manual(values = viridisLite::plasma(20), drop = FALSE)
  }
  
  p
}

select_b_dirichlet_lscv <- function(y, x, delta, b_grid, lscv_grid = LSCV_GRID, eps = EPS) {
  pi_hat <- estimate_pi_hat(x, delta, eps = eps)
  vals <- vapply(b_grid, function(b) lscv_dirichlet(y, delta, pi_hat, b, lscv_grid, eps), numeric(1))
  b_opt <- b_grid[which.min(vals)]
  list(b_opt = b_opt, pi_hat = pi_hat, lscv_values = vals)
}

run_section_1_visuals <- function() {
  message("\n=== SECTION 1: MAR + density visualizations ===")
  
  grid_vis <- generate_simplex_grid(res = GRID_RES_VIS, eps_grid = EPS_GRID_VIS)
  
  for (model_name in names(MODELS)) {
    pars <- MODELS[[model_name]]
    
    dat_vis <- generate_data(N_VIS, pars, RHO, BETA_1, P_VIS)
    
    # Figure 1 — MAR visualization
    p_mar <- make_mar_plot(dat_vis, paste("Visualisation MAR -", model_name))
    save_plotly_pdf(
      p_mar,
      file.path(FIG_DIR, paste0("mar_visualization_", tolower(model_name), ".pdf")),
      webshot_delay = PLOTLY_EXPORT_DELAY_DEFAULT
    )
    
    # True density
    f_true <- true_density_fun(grid_vis, pars)
    
    # Estimated density (Dirichlet KDE + LSCV bandwidth)
    sel  <- select_b_dirichlet_lscv(dat_vis$y, dat_vis$x, dat_vis$delta, B_GRID_CONTOUR, LSCV_GRID, EPS)
    bopt <- sel$b_opt
    pih  <- sel$pi_hat
    
    f_est <- dirichlet_kde_ipw(dat_vis$y, dat_vis$delta, pih, bopt, grid_vis)
    
    # Common breaks for true vs estimated (same legend, per model)
    rng <- range(c(f_true, f_est), finite = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
      rng <- c(0, 1)
    }
    breaks_common <- seq(rng[1], rng[2], length.out = 11)
    
    # Figure 2 — True density (theme controlled by FIG2_* toggles)
    p_true <- make_contour_plot(grid_vis, f_true, "", breaks = breaks_common, theme_fun = theme_fig2_gg)
    save_gg_pdf(
      p_true,
      file.path(FIG_DIR, paste0("true_density_", tolower(model_name), ".pdf")),
      width = 7, height = 6
    )
    
    # Figure 3 — Estimated density (theme controlled by FIG3_* toggles)
    p_est <- make_contour_plot(grid_vis, f_est, "", breaks = breaks_common, theme_fun = theme_fig3_gg)
    save_gg_pdf(
      p_est,
      file.path(FIG_DIR, paste0("estimated_density_", tolower(model_name), ".pdf")),
      width = 7, height = 6
    )
    
    message(sprintf("  %s done (b_LSCV=%.4f).", model_name, bopt))
  }
}

# ============================================================
# SECTION 2 — Tables 1–2 + Figure 4–6 (Dirichlet KDE only)
# ============================================================

run_single_replication_table <- function(n, p_obs, pars, grid_eval) {
  dat <- generate_data(n, pars, RHO, BETA_1, p_obs)
  
  pi_hat <- estimate_pi_hat(dat$x, dat$delta, eps = EPS)
  
  lscv_vals <- vapply(B_GRID_TABLE, function(b) {
    lscv_dirichlet(dat$y, dat$delta, pi_hat, b, lscv_grid = LSCV_GRID, eps = EPS)
  }, numeric(1))
  b_opt <- B_GRID_TABLE[which.min(lscv_vals)]
  
  dens_est  <- dirichlet_kde_ipw(dat$y, dat$delta, pi_hat, b_opt, grid_eval, eps = EPS)
  true_grid <- true_density_fun(grid_eval, pars, eps = EPS)
  
  list(
    ise          = mean((dens_est - true_grid)^2, na.rm = TRUE) * 0.5, # 0.5 is the area of S_2
    percent_miss = mean(dat$delta == 0) * 100,
    b_opt        = b_opt
  )
}

run_simulation_table_for_model <- function(cl, pars, model_name) {
  message("\n--- Running Monte Carlo table for: ", model_name, " ---")
  message("Using ", NUM_CORES, " cores (detectCores()-1).")
  
  grid_eval <- generate_simplex_grid(res = GRID_RES_TABLE, eps_grid = EPS_GRID_TABLE)
  
  # Export large objects once (avoid serializing them at every replication)
  parallel::clusterExport(cl, varlist = c("pars", "grid_eval"), envir = environment())
  
  raw_all    <- data.frame()
  summary_df <- data.frame()
  
  for (n in N_VALUES_TABLE) {
    for (p_obs in P_VALUES_TABLE) {
      message(sprintf("  n = %d, Missing Rate = %.0f%%", n, (1 - p_obs) * 100))
      
      # Deterministic, unique seeds per (model,n,p,rep)
      off_model <- ifelse(model_name == "Model1", 10000000L, 20000000L)
      off_n     <- 100000L * which(N_VALUES_TABLE == n)
      off_p     <- 1000L   * which(P_VALUES_TABLE == p_obs)
      seeds     <- SEED + off_model + off_n + off_p + seq_len(MC_REPS_TABLE)
      
      rep_res <- parallel::parLapply(
        cl,
        X = seeds,
        fun = function(seed, n, p_obs) {
          set.seed(seed)
          run_single_replication_table(n, p_obs, pars, grid_eval)
        },
        n = n,
        p_obs = p_obs
      )
      
      ise_vals  <- vapply(rep_res, `[[`, numeric(1), "ise")
      miss_vals <- vapply(rep_res, `[[`, numeric(1), "percent_miss")
      b_vals    <- vapply(rep_res, `[[`, numeric(1), "b_opt")
      
      raw_all <- rbind(
        raw_all,
        data.frame(
          model        = model_name,
          n            = n,
          p            = p_obs,
          replication  = seq_len(MC_REPS_TABLE),
          ISE          = ise_vals,
          percent_miss = miss_vals,
          b_opt        = b_vals
        )
      )
      
      summary_df <- rbind(
        summary_df,
        data.frame(
          n               = n,
          p               = p_obs,
          percent_missing = round(mean(miss_vals), 2),
          mean_ISE        = round(mean(ise_vals), 4),
          median_ISE      = round(median(ise_vals), 4),
          SD_ISE          = round(sd(ise_vals), 4),
          IQR_ISE         = round(IQR(ise_vals), 4),
          mean_b          = round(mean(b_vals), 4)
        )
      )
    }
  }
  
  # ---- Save outputs (RAW CSV -> SUMMARY CSV -> LaTeX) ----
  raw_csv     <- file.path(TABLE_DIR, paste0("raw_ise_", tolower(model_name), ".csv"))
  summary_csv <- file.path(TABLE_DIR, paste0("summary_ise_", tolower(model_name), ".csv"))
  tex_file    <- file.path(TABLE_DIR, paste0("table_lscv_", tolower(model_name), ".tex"))
  
  write.csv(raw_all, raw_csv, row.names = FALSE)
  write.csv(summary_df, summary_csv, row.names = FALSE)
  
  if (model_name == "Model1") {
    caption <- "The mean, median, standard deviation, and interquartile range of 1000 ISEs in Model I for the IPW Dirichlet KDE, with sample sizes $n\\in \\{100, 200, 400, 800\\}$ and missing rates of $5\\%$, $10\\%$, $20\\%$, and $40\\%$."
    label   <- "tab:results_LSCV_Model1"
    comment <- "TABLE MODEL 1"
  } else {
    caption <- "The mean, median, standard deviation, and interquartile range of 1000 ISEs in Model II for the IPW Dirichlet KDE, with sample sizes $n\\in \\{100, 200, 400, 800\\}$ and missing rates of $5\\%$, $10\\%$, $20\\%$, and $40\\%$."
    label   <- "tab:results_LSCV_Model2"
    comment <- "TABLE MODEL 2"
  }
  
  write_latex_table(
    summary_df    = summary_df,
    file_tex      = tex_file,
    caption       = caption,
    label         = label,
    table_comment = comment
  )
  
  list(raw = raw_all, summary = summary_df)
}

make_mean_ise_plot <- function(tab, model_name) {
  # Calculate expected missing percentage: (1 - p_obs) * 100
  tab <- tab %>% mutate(missing_level = (1 - p) * 100)
  
  ggplot(tab, aes(x = missing_level, y = mean_ISE, color = factor(n), group = n)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      x = "Missing Rate (%)",
      y = "Mean ISE",
      color = "      n"
    ) +
    theme_fig4_gg()
}

make_ise_metrics_plot <- function(tab, model_name, theme_fun) {
  # Calculate expected missing percentage: (1 - p_obs) * 100
  tab <- tab %>% mutate(missing_level = (1 - p) * 100)
  
  p1 <- ggplot(tab, aes(missing_level, mean_ISE, color = factor(n), group = n)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    labs(x = "Missing Rate (%)", y = "Mean ISE", color = "         n") +
    theme_fun()
  
  p2 <- ggplot(tab, aes(missing_level, median_ISE, color = factor(n), group = n)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    labs(x = "Missing Rate (%)", y = "Median ISE", color = "         n") +
    theme_fun()
  
  p3 <- ggplot(tab, aes(missing_level, SD_ISE, color = factor(n), group = n)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    labs(x = "Missing Rate (%)", y = "SD ISE", color = "         n") +
    theme_fun()
  
  p4 <- ggplot(tab, aes(missing_level, IQR_ISE, color = factor(n), group = n)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    labs(x = "Missing Rate (%)", y = "IQR ISE", color = "         n") +
    theme_fun()
  
  (p1 + p2) / (p3 + p4)
}

run_section_2_tables_and_plots <- function(cl) {
  message("\n=== SECTION 2: Tables 1–2 + Mean/metrics plots ===")
  
  res1 <- run_simulation_table_for_model(cl, MODELS$Model1, "Model1")
  res2 <- run_simulation_table_for_model(cl, MODELS$Model2, "Model2")
  
  table_model1 <- res1$summary
  table_model2 <- res2$summary
  
  # Figure 4: mean ISE lines
  p_mean_1 <- make_mean_ise_plot(table_model1, "Model 1")
  p_mean_2 <- make_mean_ise_plot(table_model2, "Model 2")
  
  save_gg_pdf(p_mean_1, file.path(FIG_DIR, "mean_ise_model1.pdf"), width = 7.5, height = 5.5)
  save_gg_pdf(p_mean_2, file.path(FIG_DIR, "mean_ise_model2.pdf"), width = 7.5, height = 5.5)
  
  # Figure 5–6: ISE metrics (Model 1 -> Figure 5 toggles; Model 2 -> Figure 6 toggles)
  p_metrics_1 <- make_ise_metrics_plot(table_model1, "Model 1", theme_fun = theme_fig5_gg)
  p_metrics_2 <- make_ise_metrics_plot(table_model2, "Model 2", theme_fun = theme_fig6_gg)
  
  save_gg_pdf(p_metrics_1, file.path(FIG_DIR, "ise_metrics_model1.pdf"), width = 10, height = 8)
  save_gg_pdf(p_metrics_2, file.path(FIG_DIR, "ise_metrics_model2.pdf"), width = 10, height = 8)
  
  invisible(list(Model1 = res1, Model2 = res2))
}

# ============================================================
# SECTION 3 — Dirichlet vs ALR vs ILR (LSCV bandwidths)
#   -> Model1 vs Model2, feasible only
#   -> Figure 7 baseline deleted; renumbered figures:
#        Figure 7: sample size variation
#        Figure 8: missingness variation
#        Figure 9: joint effect (n × missing)
# ============================================================

comparison_one_rep <- function(n, p_obs, pars, s_grid, f_true_grid) {
  dat <- generate_data(n, pars, RHO, BETA_1, p_obs)
  
  Y     <- dat$y
  X     <- dat$x
  delta <- dat$delta
  
  # Estimated pi for feasible weights only
  pi_hat <- estimate_pi_hat(X, delta, eps = EPS)
  
  # --- Dirichlet b via LSCV (feasible weights) ---
  lscv_b <- vapply(B_GRID_COMP, function(b) {
    lscv_dirichlet(Y, delta, pi_hat, b, lscv_grid = LSCV_GRID, eps = EPS)
  }, numeric(1))
  b_opt <- B_GRID_COMP[which.min(lscv_b)]
  
  # --- ALR h via LSCV ---
  idx_obs <- which(delta == 1)
  h0_alr <- h_reference(Y[idx_obs, , drop = FALSE], alr_transform, eps = EPS)
  h_grid_alr <- h0_alr * H_GRID_FACTOR
  lscv_h_alr <- vapply(h_grid_alr, function(h) {
    lscv_logratio(Y, delta, pi_hat, h, LSCV_GRID, alr_transform, jacobian_dzdy_alr, eps = EPS)
  }, numeric(1))
  h_opt_alr <- h_grid_alr[which.min(lscv_h_alr)]
  
  # --- ILR h via LSCV ---
  h0_ilr <- h_reference(Y[idx_obs, , drop = FALSE], ilr_transform, eps = EPS)
  h_grid_ilr <- h0_ilr * H_GRID_FACTOR
  lscv_h_ilr <- vapply(h_grid_ilr, function(h) {
    lscv_logratio(Y, delta, pi_hat, h, LSCV_GRID, ilr_transform, jacobian_dzdy_ilr, eps = EPS)
  }, numeric(1))
  h_opt_ilr <- h_grid_ilr[which.min(lscv_h_ilr)]
  
  # --- Density estimates on s_grid (feasible only) ---
  f_dir <- dirichlet_kde_ipw(Y, delta, pi_hat,  b_opt,     s_grid, eps = EPS)
  f_alr <- logratio_kde_ipw( Y, delta, pi_hat,  h_opt_alr, s_grid, alr_transform, jacobian_dzdy_alr, eps = EPS)
  f_ilr <- logratio_kde_ipw( Y, delta, pi_hat,  h_opt_ilr, s_grid, ilr_transform, jacobian_dzdy_ilr, eps = EPS)
  
  ise <- colMeans(
    (cbind(f_dir, f_alr, f_ilr) - f_true_grid)^2 * 0.5, # 0.5 is the area of S_2
    na.rm = TRUE
  )
  
  names(ise) <- c("Dirichlet", "ALR", "ILR")
  ise
}

run_section_3_comparisons <- function(cl) {
  message("\n=== SECTION 3: Dirichlet vs ALR vs ILR comparisons (LSCV) ===")
  message("Using ", NUM_CORES, " cores (detectCores()-1).")
  
  for (model_name in names(MODELS)) {
    message("\n--- Comparative study for: ", model_name, " (feasible only) ---")
    
    pars <- MODELS[[model_name]]
    s_grid <- generate_simplex_grid(res = GRID_RES_COMP, eps_grid = EPS_GRID_COMP)
    f_true_grid <- true_density_fun(s_grid, pars)
    
    # Export once for this model (avoid re-serializing large objects each rep)
    parallel::clusterExport(cl, varlist = c("pars", "s_grid", "f_true_grid"), envir = environment())
    
    off_model <- if (model_name == "Model1") 0L else 5000000L
    
    # ----------------------------------------------------------
    # Figure 7: sample size variation — PARALLEL over reps (per n)
    # ----------------------------------------------------------
    message("  Figure 7: impact of sample size (parallel MC reps)")
    
    ISE_all_n <- data.frame()
    
    for (n in N_VALUES_COMP_N) {
      off_n <- 100000L * which(N_VALUES_COMP_N == n)
      seeds <- SEED + off_model + 31000000L + off_n + seq_len(MC_REPS_COMP_N)
      
      res_list_n <- parallel::parLapply(
        cl,
        X = seeds,
        fun = function(seed, n, p_obs) {
          set.seed(seed)
          comparison_one_rep(n, p_obs, pars, s_grid, f_true_grid)
        },
        n = n,
        p_obs = P_COMP_N_FIXED
      )
      
      M  <- do.call(rbind, res_list_n)
      MC <- nrow(M)
      
      ISE_all_n <- rbind(
        ISE_all_n,
        data.frame(
          ISE = c(M[, "Dirichlet"], M[, "ALR"], M[, "ILR"]),
          Estimator = rep(c("Dirichlet", "ALR", "ILR"), each = MC),
          SampleSize = paste0("n=", n)
        )
      )
      
      message("    Finished n = ", n)
    }
    
    ISE_all_n$Estimator  <- factor(ISE_all_n$Estimator,  levels = c("ALR", "ILR", "Dirichlet"))
    ISE_all_n$SampleSize <- factor(ISE_all_n$SampleSize, levels = paste0("n=", N_VALUES_COMP_N))
    
    lims_n <- ISE_all_n %>%
      group_by(SampleSize, Estimator) %>%
      summarize(
        lower = boxplot.stats(ISE)$stats[1],
        upper = boxplot.stats(ISE)$stats[5],
        .groups = "drop"
      ) %>%
      summarize(min_y = min(lower), max_y = max(upper))
    
    p_n <- ggplot(ISE_all_n, aes(x = SampleSize, y = ISE, fill = Estimator)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      coord_cartesian(ylim = c(lims_n$min_y, lims_n$max_y)) +
      labs(x = "Sample size", y = "ISE") +
      theme_fig7_gg() +
      theme(
        legend.position = c(0.99, 0.99),       
        legend.justification = c("right", "top"), 
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA) 
      )
    
    save_gg_pdf(
      p_n,
      file.path(FIG_DIR, paste0(tolower(model_name), "_comparison_n.pdf")),
      width = 8.5, height = 5.5
    )
    
    # ----------------------------------------------------------
    # Figure 8: missingness variation — PARALLEL over reps (per miss)
    # ----------------------------------------------------------
    message("  Figure 8: impact of missing rate (parallel MC reps)")
    
    ISE_all_missing <- data.frame()
    
    for (miss in MISSING_LEVELS) {
      p_obs <- 1 - miss
      off_miss <- 100000L * which(MISSING_LEVELS == miss)
      seeds <- SEED + off_model + 32000000L + off_miss + seq_len(MC_REPS_COMP_MISS)
      
      res_list_m <- parallel::parLapply(
        cl,
        X = seeds,
        fun = function(seed, n, p_obs) {
          set.seed(seed)
          comparison_one_rep(n, p_obs, pars, s_grid, f_true_grid)
        },
        n = N_COMP_MISS_FIXED,
        p_obs = p_obs
      )
      
      M  <- do.call(rbind, res_list_m)
      MC <- nrow(M)
      
      ISE_all_missing <- rbind(
        ISE_all_missing,
        data.frame(
          ISE = c(M[, "Dirichlet"], M[, "ALR"], M[, "ILR"]),
          Estimator = rep(c("Dirichlet", "ALR", "ILR"), each = MC),
          Missing = paste0(miss * 100, "%")
        )
      )
      
      message("    Finished missing = ", miss * 100, "%")
    }
    
    ISE_all_missing$Estimator <- factor(ISE_all_missing$Estimator, levels = c("ALR", "ILR", "Dirichlet"))
    ISE_all_missing$Missing   <- factor(ISE_all_missing$Missing,   levels = paste0(MISSING_LEVELS * 100, "%"))
    
    lims_m <- ISE_all_missing %>%
      group_by(Missing, Estimator) %>%
      summarize(
        lower = boxplot.stats(ISE)$stats[1],
        upper = boxplot.stats(ISE)$stats[5],
        .groups = "drop"
      ) %>%
      summarize(min_y = min(lower), max_y = max(upper))
    
    p_m <- ggplot(ISE_all_missing, aes(x = Missing, y = ISE, fill = Estimator)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      coord_cartesian(ylim = c(lims_m$min_y, lims_m$max_y)) +
      labs(x = "Missing rate", y = "ISE") +
      theme_fig8_gg() +
      theme(
        legend.position = c(0.99, 0.99),       
        legend.justification = c("right", "top"), 
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA) 
      )
    
    save_gg_pdf(
      p_m,
      file.path(FIG_DIR, paste0(tolower(model_name), "_comparison_missing.pdf")),
      width = 8, height = 5.5
    )
    
    # ----------------------------------------------------------
    # Figure 9: joint effect (n × missing) — PARALLEL over reps (per cell)
    # ----------------------------------------------------------
    message("  Figure 9: joint effect (n × missing) (parallel MC reps per cell)")
    
    ISE_global <- data.frame()
    
    for (n in N_VALUES_COMP_GLOBAL) {
      for (miss in MISSING_LEVELS_GLOBAL) {
        p_obs <- 1 - miss
        
        off_n    <- 1000000L * which(N_VALUES_COMP_GLOBAL == n)
        off_miss <- 10000L   * which(MISSING_LEVELS_GLOBAL == miss)
        seeds <- SEED + off_model + 33000000L + off_n + off_miss + seq_len(MC_REPS_COMP_GLOBAL)
        
        res_list_g <- parallel::parLapply(
          cl,
          X = seeds,
          fun = function(seed, n, p_obs) {
            set.seed(seed)
            comparison_one_rep(n, p_obs, pars, s_grid, f_true_grid)
          },
          n = n,
          p_obs = p_obs
        )
        
        M  <- do.call(rbind, res_list_g)
        MC <- nrow(M)
        
        ISE_global <- rbind(
          ISE_global,
          data.frame(
            ISE = c(M[, "Dirichlet"], M[, "ALR"], M[, "ILR"]),
            Estimator = rep(c("Dirichlet", "ALR", "ILR"), each = MC),
            SampleSize = paste0("n=", n),
            Missing = paste0(miss * 100, "%")
          )
        )
        
        message(sprintf("    Finished n=%d, missing=%s", n, paste0(miss * 100, "%")))
      }
    }
    
    ISE_global$Estimator  <- factor(ISE_global$Estimator,  levels = c("ALR", "ILR", "Dirichlet"))
    ISE_global$SampleSize <- factor(ISE_global$SampleSize, levels = paste0("n=", N_VALUES_COMP_GLOBAL))
    ISE_global$Missing    <- factor(ISE_global$Missing,    levels = paste0(MISSING_LEVELS_GLOBAL * 100, "%"))
    
    lims_g <- ISE_global %>%
      group_by(SampleSize, Missing, Estimator) %>%
      summarize(
        lower = boxplot.stats(ISE)$stats[1],
        upper = boxplot.stats(ISE)$stats[5],
        .groups = "drop"
      ) %>%
      summarize(min_y = min(lower), max_y = max(upper))
    
    p_g <- ggplot(ISE_global, aes(x = SampleSize, y = ISE, fill = Estimator)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      facet_wrap(~Missing) +
      coord_cartesian(ylim = c(lims_g$min_y, lims_g$max_y)) +
      labs(x = "Sample size", y = "ISE") +
      theme_fig9_gg() +
      theme(
        axis.text.x = element_text(size = FIG9_AXIS_TEXT_X_SIZE, angle = 45, hjust = 1),
        legend.position = c(0.99, 0.99),       
        legend.justification = c("right", "top"), 
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
      )
    
    save_gg_pdf(
      p_g,
      file.path(FIG_DIR, paste0(tolower(model_name), "_comparison_missing_n.pdf")),
      width = 10, height = 6
    )
  }
  
  invisible(TRUE)
}

# ============================================================
# SECTION 4 — NHANES application (Dirichlet KDE + LSCV)
#   Figure 10a: 2D raster (viridis) — PDF-safe (NO TITLE)
#   Figure 10b: 3D density (STATIC PDF-safe, no plotly/WebGL)
#              -> 4 different viewing angles (4 PDFs)
#              -> NO TITLE (legend/color key label is fine)
# ============================================================

run_section_4_application <- function() {
  message("\n=== SECTION 4: NHANES application ===")
  
  # --- Optional: crop PDFs to reduce white margins (requires TeX pdfcrop in PATH) ---
  maybe_pdfcrop_local <- function(pdf_file) {
    if (!file.exists(pdf_file)) return(invisible(FALSE))
    pdfcrop_bin <- Sys.which("pdfcrop")
    if (!nzchar(pdfcrop_bin)) return(invisible(FALSE))
    
    tmp_out <- tempfile(fileext = ".pdf")
    ok <- TRUE
    tryCatch({
      system2(pdfcrop_bin, c(shQuote(pdf_file), shQuote(tmp_out)), stdout = TRUE, stderr = TRUE)
    }, error = function(e) ok <<- FALSE)
    
    if (ok && file.exists(tmp_out)) {
      file.copy(tmp_out, pdf_file, overwrite = TRUE)
      unlink(tmp_out)
      return(invisible(TRUE))
    }
    invisible(FALSE)
  }
  
  # --- NHANES download ---
  if (!requireNamespace("nhanesA", quietly = TRUE)) {
    message("Package 'nhanesA' not found. Install with install.packages('nhanesA') or adapt to local XPT files.")
    return(invisible(FALSE))
  }
  suppressPackageStartupMessages(library(nhanesA))
  
  cbc_name <- paste0("CBC_", NHANES_CYCLE_SUFFIX)
  bmx_name <- paste0("BMX_", NHANES_CYCLE_SUFFIX)  # NEW: Body measures (BMI)
  
  cbc <- nhanesA::nhanes(cbc_name)
  bmx <- nhanesA::nhanes(bmx_name)
  
  cbc_keep <- cbc %>%
    dplyr::select(SEQN, LBXNEPCT, LBXLYPCT, LBXMOPCT, LBXEOPCT, LBXBAPCT)
  
  bmx_keep <- bmx %>%
    dplyr::select(SEQN, BMXBMI) %>%
    # IMPORTANT: ensure X has no NA (estimate_pi_hat uses dist(), which breaks with NA)
    dplyr::filter(!is.na(BMXBMI))
  
  # Merge by respondent sequence number (SEQN).
  # Use left_join so participants with BMI but missing CBC differential are kept with delta = 0.
  dat <- bmx_keep %>%
    dplyr::left_join(cbc_keep, by = "SEQN") %>%
    dplyr::mutate(
      OthersPCT = LBXMOPCT + LBXEOPCT + LBXBAPCT
    )
  
  n_total <- nrow(dat)
  if (n_total == 0) {
    message("  NHANES: no records after filtering for non-missing BMI; aborting Section 4.")
    return(invisible(FALSE))
  }
  
  # --- Missingness indicator (unit-level for CBC differential) ---
  # delta = 1 iff *all five* CBC differential percentage components are observed.
  delta <- as.integer(stats::complete.cases(
    dat$LBXNEPCT, dat$LBXLYPCT, dat$LBXMOPCT, dat$LBXEOPCT, dat$LBXBAPCT
  ))
  
  n_obs  <- sum(delta == 1)
  n_miss <- n_total - n_obs
  message(sprintf("  NHANES: n = %d (BMI observed); delta=1: %d, delta=0: %d",
                  n_total, n_obs, n_miss))
  
  if (n_obs == 0) {
    message("  NHANES: no observed CBC differential compositions (delta=1 count is 0); skipping Section 4.")
    return(invisible(FALSE))
  }
  
  # --- Build Y (compositions) and X (covariates) ---
  # NHANES variables are PERCENTAGES (0–100); convert to proportions (/100).
  Y_raw <- dat %>%
    dplyr::transmute(
      s1 = LBXNEPCT / 100,       # Neutrophils proportion
      s2 = LBXLYPCT / 100,       # Lymphocytes proportion
      s3 = OthersPCT / 100       # Others = MO + EO + BA, as proportion
    )
  
  # Optional robustness: mark obviously invalid observed rows as missing
  # (should be rare/nonexistent in NHANES, but protects downstream Dirichlet code)
  ok_range <- with(Y_raw,
                   is.finite(s1) & is.finite(s2) & is.finite(s3) &
                     (s1 >= 0) & (s2 >= 0) & (s3 >= 0) &
                     (s1 <= 1) & (s2 <= 1) & (s3 <= 1))
  delta <- as.integer((delta == 1) & ok_range)
  n_obs2 <- sum(delta == 1)
  if (n_obs2 != n_obs) {
    message(sprintf("  NHANES: dropped %d invalid observed CBC rows after range checks.",
                    n_obs - n_obs2))
    n_obs <- n_obs2
    n_miss <- n_total - n_obs
  }
  
  # p = 1 continuous auxiliary covariate
  X <- dat %>%
    dplyr::transmute(
      bmi = as.numeric(BMXBMI)
    ) %>%
    as.matrix()
  
  Y <- as.matrix(Y_raw)
  
  # Normalize observed compositions only (enforce interior simplex for Dirichlet KDE)
  if (any(delta == 1)) {
    Y[delta == 1, ] <- normalize_comp(Y[delta == 1, , drop = FALSE], eps = EPS)
  }
  
  # --- Estimate pi(x) ---
  # NOTE: This is an unweighted NHANES analysis (treating data as i.i.d.).
  pi_hat <- estimate_pi_hat(X, delta, eps = EPS)
  
  # --- LSCV for Dirichlet bandwidth b ---
  grid_app <- generate_simplex_grid(res = GRID_RES_APP, eps_grid = EPS_GRID_APP)
  
  lscv_vals <- vapply(B_GRID_APP, function(b) {
    lscv_dirichlet(Y, delta, pi_hat, b, lscv_grid = LSCV_GRID, eps = EPS)
  }, numeric(1))
  
  b_opt <- B_GRID_APP[which.min(lscv_vals)]
  message(sprintf("  NHANES: b_opt (LSCV) = %.4f", b_opt))
  
  # --- Density estimate on grid ---
  f_app <- dirichlet_kde_ipw(Y, delta, pi_hat, b_opt, grid_app, eps = EPS)
  
  df_plot <- data.frame(
    neutro  = grid_app$s1,
    lympho  = grid_app$s2,
    density = f_app
  )
  
  # ============================================================
  # Figure 10a — 2D raster plot (tight margins, viridis) — NO TITLE
  # ============================================================
  
  p2d <- ggplot(df_plot, aes(x = neutro, y = lympho, fill = density)) +
    geom_raster(interpolate = TRUE, na.rm = TRUE) +
    scale_fill_viridis_c(option = "viridis") +
    coord_fixed(expand = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0))) +
    labs(
      x = "Neutrophils proportion",
      y = "Lymphocytes proportion",
      fill = "Density"
    ) +
    theme_fig10a_gg() +
    theme(
      # reduce white space without truncating
      plot.margin       = ggplot2::margin(6, 6, 6, 6, unit = "pt"),
      legend.margin     = ggplot2::margin(4, 4, 4, 4, unit = "pt"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
    )
  
  out2d_pdf <- file.path(FIG_DIR, "leukocyte_kde_2d.pdf")
  save_gg_pdf(p2d, out2d_pdf, width = 9.0, height = 7.5)
  # maybe_pdfcrop_local(out2d_pdf)
  
  # ============================================================
  # Figure 10b — 3D density (STATIC, PDF-safe) — 4 viewpoints
  #   Primary: surf3D via akima interpolation
  #   Fallback: scatter3D
  #   NO TITLE
  # ============================================================
  
  # Clean non-finite densities
  df3 <- df_plot %>%
    dplyr::filter(is.finite(neutro), is.finite(lympho), is.finite(density))
  
  if (nrow(df3) == 0) {
    message("  NHANES: no finite density values on grid; skipping Figure 10b.")
    return(invisible(TRUE))
  }
  
  idx_max <- which.max(replace(df3$density, !is.finite(df3$density), -Inf))
  mode_neutro  <- df3$neutro[idx_max]
  mode_lympho  <- df3$lympho[idx_max]
  mode_density <- df3$density[idx_max]
  mode_others  <- max(0, 1 - mode_neutro - mode_lympho)
  
  message(sprintf(
    "  Approximate mode on grid: (%.2f, %.2f, %.2f), density=%.4f",
    mode_neutro, mode_lympho, mode_others, mode_density
  ))
  
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    message("Package 'plot3D' not installed. Install with install.packages('plot3D'). Skipping Figure 10b.")
    return(invisible(TRUE))
  }
  
  # Palette
  n_col <- 200
  cols <- if (requireNamespace("viridisLite", quietly = TRUE)) {
    viridisLite::viridis(n_col)
  } else {
    grDevices::colorRampPalette(c("navy", "blue", "green", "yellow"))(n_col)
  }
  
  # ---- Precompute smooth surface interpolation ONCE (if possible) ----
  did_surface_possible <- FALSE
  xmat_surf <- ymat_surf <- zmat_surf <- NULL
  
  if (requireNamespace("akima", quietly = TRUE)) {
    xo <- seq(0, 1, length.out = 300)
    yo <- seq(0, 1, length.out = 300)
    
    interp <- suppressWarnings(akima::interp(
      x = df3$neutro,
      y = df3$lympho,
      z = df3$density,
      xo = xo,
      yo = yo,
      duplicate = "mean",
      linear = FALSE    # smooth interpolation
    ))
    
    if (!is.null(interp) && !is.null(interp$z)) {
      zmat <- interp$z
      
      # Build x/y as MATRICES (required by plot3D::surf3D)
      xmat <- matrix(rep(interp$x, times = length(interp$y)),
                     nrow = length(interp$x), ncol = length(interp$y))
      ymat <- matrix(rep(interp$y, each = length(interp$x)),
                     nrow = length(interp$x), ncol = length(interp$y))
      
      # Mask outside simplex (neutro + lympho <= 1)
      mask <- (xmat + ymat) <= 1
      zmat[!mask] <- NA_real_
      
      xmat_surf <- xmat
      ymat_surf <- ymat
      zmat_surf <- zmat
      did_surface_possible <- TRUE
    }
  }
  
  # ---- 4 distinct viewpoints ----
  views <- list(
    theta045_phi20 = list(theta = 45,  phi = 20),
    theta135_phi20 = list(theta = 135, phi = 20),
    theta225_phi20 = list(theta = 225, phi = 20),
    theta315_phi20 = list(theta = 315, phi = 20)
  )
  
  for (vname in names(views)) {
    theta_v <- views[[vname]]$theta
    phi_v   <- views[[vname]]$phi
    
    out3d_pdf <- file.path(FIG_DIR, paste0("leukocyte_kde_3d_mode_", vname, ".pdf"))
    
    grDevices::pdf(out3d_pdf, width = 9.0, height = 7.5, onefile = TRUE, useDingbats = FALSE)
    
    op <- par(no.readonly = TRUE)
    par(
      mar = c(4.4, 4.4, 3.5, 5.2),
      cex.lab  = FIG10B_CEX_LAB,
      cex.axis = FIG10B_CEX_AXIS,
      font.lab  = 1,
      font.main = 1
    )
    
    if (did_surface_possible) {
      plot3D::surf3D(
        x = xmat_surf,
        y = ymat_surf,
        z = zmat_surf,
        colvar = zmat_surf,
        col = cols,
        border = NA,
        theta = theta_v,
        phi = phi_v,
        xlab = "Neutrophils proportion",
        ylab = "Lymphocytes proportion",
        zlab = "Density",
        clab = "Density",
        ticktype = "detailed",
        bty = "b2",
        colkey = list(
          side = 4,
          length = 0.6,
          width  = 0.6,
          cex.axis = FIG10B_COLKEY_CEX_AXIS,
          cex.clab = FIG10B_COLKEY_CEX_CLAB
        )
      )
    } else {
      plot3D::scatter3D(
        x = df3$neutro,
        y = df3$lympho,
        z = df3$density,
        colvar = df3$density,
        col = cols,
        pch = 16,
        cex = 0.55,
        alpha = 0.9,
        theta = theta_v,
        phi = phi_v,
        xlab = "Neutrophils proportion",
        ylab = "Lymphocytes proportion",
        zlab = "Density",
        clab = "Density",
        ticktype = "detailed",
        bty = "b2",
        colkey = list(
          side = 4,
          length = 0.6,
          width  = 0.6,
          cex.axis = FIG10B_COLKEY_CEX_AXIS,
          cex.clab = FIG10B_COLKEY_CEX_CLAB
        )
      )
    }
    
    # Small vertical offset so point/label float above surface
    z_offset <- mode_density * 0.0035
    
    plot3D::points3D(
      x = mode_neutro,
      y = mode_lympho,
      z = mode_density + z_offset,
      add = TRUE,
      col = "red",
      pch = 19,
      cex = 2.0
    )
    
    plot3D::text3D(
      x = mode_neutro,
      y = mode_lympho,
      z = mode_density + z_offset + (mode_density * 0.05),
      add = TRUE,
      labels = sprintf("Mode\n(%.2f, %.2f)", mode_neutro, mode_lympho),
      col = "red",
      cex = FIG10B_MODE_LABEL_CEX
    )
    
    par(op)
    grDevices::dev.off()
    
    # maybe_pdfcrop_local(out3d_pdf)
    message(sprintf("  Saved 3D view: %s", basename(out3d_pdf)))
  }
  
  invisible(TRUE)
}

# ============================================================
# MAIN — Run sections in the same order as the paper
# ============================================================

main <- function() {
  global_start <- Sys.time()
  message("\n=== GLOBAL TIMER STARTED ===")
  
  cl <- NULL
  needs_parallel <- RUN_SECTION_2_TABLES || RUN_SECTION_3_COMPARISON
  
  # Safety stop on exit (in case of error)
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
  
  # ----------------------------
  # Parallel cluster init timer
  # ----------------------------
  if (needs_parallel) {
    t_init <- Sys.time()
    message("\n=== Initializing parallel cluster ===")
    message("Detected cores: ", cores_detected, " | Using cores: ", NUM_CORES)
    
    cl <- setup_parallel_cluster(num_cores = NUM_CORES)
    
    timer_report("Parallel cluster initialization", t_init)
  }
  
  # ----------------------------
  # Section 1 timer
  # ----------------------------
  if (RUN_SECTION_1_VISUALS) {
    t1 <- Sys.time()
    run_section_1_visuals()
    timer_report("SECTION 1 (Visuals)", t1)
  }
  
  # ----------------------------
  # Section 2 timer
  # ----------------------------
  if (RUN_SECTION_2_TABLES) {
    if (is.null(cl)) stop("Parallel cluster not initialized but RUN_SECTION_2_TABLES=TRUE.")
    t2 <- Sys.time()
    run_section_2_tables_and_plots(cl)
    timer_report("SECTION 2 (Tables + Figures 4–6)", t2)
  }
  
  # ----------------------------
  # Section 3 timer
  # ----------------------------
  if (RUN_SECTION_3_COMPARISON) {
    if (is.null(cl)) stop("Parallel cluster not initialized but RUN_SECTION_3_COMPARISON=TRUE.")
    t3 <- Sys.time()
    run_section_3_comparisons(cl)
    timer_report("SECTION 3 (Comparisons Figures 7–9)", t3)
  }
  
  # ----------------------------
  # Section 4 timer
  # ----------------------------
  if (RUN_SECTION_4_APPLICATION) {
    t4 <- Sys.time()
    run_section_4_application()
    timer_report("SECTION 4 (NHANES Application Figure 10)", t4)
  }
  
  # ----------------------------
  # Parallel cluster stop timer
  # ----------------------------
  if (!is.null(cl)) {
    t_stop <- Sys.time()
    try(parallel::stopCluster(cl), silent = TRUE)
    cl <- NULL
    timer_report("Parallel cluster shutdown", t_stop)
  }
  
  timer_report("TOTAL runtime", global_start)
  message("=== GLOBAL TIMER ENDED ===\n")
  
  invisible(TRUE)
}

main()

