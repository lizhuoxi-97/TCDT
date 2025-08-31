pkg_all <- c(
  "ggplot2", "gtable", "dplyr", "forcats", "ggthemes", "glue", "openxlsx",
  "purrr", "reshape2", "scales", "ggh4x", "rlang", "gridExtra"
)
for (pkg in pkg_all) {
  suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
}

source("numerical/utils.R")

# Parameterized setting ----
args <- commandArgs(TRUE)
if (length(args) > 0) {
  type_args <- "shell"
} else if (length(args) == 0) {
  type_args <- "Rscript"
}

if (type_args == "shell") {
  example <- args[1]
} else if (type_args == "Rscript") {
  example <- "local__univariate"
  # example <- "local__multivariate"
  # example <- "global__covariate_shift__univariate"
  # example <- "global__covariate_shift__multivariate"
  # example <- "global__label_shift"
  # example <- "realdata__airfoil"
}

if (example == "local__univariate") {
  # Simulations ----
  ## Local - Univariate ----
  setting_dgp <- "different_mean"
  n1 <- n2 <- 150
  p <- 1
  methods <- c("CED", "CMMD")
  x0_all <- round(seq(-1.0, 1.0, 0.1), digits = 1)
  dist_x <- "normt"
  settings_m <- c("x")
  shift_x <- "true"
  num_simu <- 1000
  mode_local <- "parallel"

  power_all <- array(
    dim = c(length(x0_all), length(methods), length(settings_m)),
    dimnames = list(format(x0_all, trim = TRUE, nsmall = 1), methods, settings_m)
  )
  df_plot <- df_signal <- data.frame()
  for (setting_m in settings_m) {
    for (x0 in x0_all) {
      obj_result <- env()
      load(
        file.path(
          "output", "results", "simulations", "local", "univariate", "different_mean",
          glue(
            "m_{setting_m}",
            "n_{n1}_{n2}",
            "p_{p}",
            "dist_x_{dist_x}",
            "shift_x_{shift_x}",
            "mode_local_{mode_local}",
            "x0_{x0}",
            "num_simu_{num_simu}.RData",
            .sep = "__"
          )
        ),
        obj_result
      )
      reject <- obj_result$reject[, methods, , drop = FALSE]
      power_all[format(x0, trim = TRUE, nsmall = 1), , setting_m] <- apply(reject, c(2, 3), mean)
    }

    form_m1 <- obj_result$form_m1
    form_m2 <- obj_result$form_m2
    diff_m <- function(x) abs(form_m1(x) - form_m2(x))

    df_plot <- rbind(
      df_plot,
      melt(
        power_all[, , setting_m, drop = FALSE], varnames = c("x0", "method", "setting"), value.name = "power"
      ) %>% mutate(setting_m = setting_m, signal_orig = diff_m(x0), signal = rescale(signal_orig))
    )

    x0_grid <- seq(min(x0_all), max(x0_all), 0.001)
    df_signal <- rbind(
      df_signal,
      data.frame(
        x0 = x0_grid,
        sample1 = form_m1(x0_grid),
        sample2 = form_m2(x0_grid),
        signal = diff_m(x0_grid),
        setting_m = setting_m
      )
    )
  }


  ### power table ----
  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  write.csv(
    power_all,
    file = file.path(dir, glue("power__{example}.csv"))
  )


  ### power plot ----
  fig <- df_plot %>%
    ggplot() +
    geom_line(aes(x = x0, y = power, linetype = method, color = method)) +
    geom_line(aes(x = x0, y = signal, linetype = "signal", color = "signal")) +
    scale_linetype_manual(
      "", values = c("CMMD" = 2, "CED" = 1, "signal" = 3),
      breaks = c("signal", "CED", "CMMD")
    ) +
    scale_color_manual(
      "", values = c("CMMD" = "#00AAFF", "CED" = "#F8766D", "signal" = "black"),
      breaks = c("signal", "CED", "CMMD")
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1)
    ) +
    xlab("Local x") +
    ylab("Power") +
    my_theme_base() +
    # theme(legend.position = "bottom")
    theme(legend.position = "right")

  dir <- file.path("output", "figures")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  ggsave(
    file.path(dir, glue("power__{example}.pdf")),
    fig,
    # width = 8, height = 4
    width = 8, height = 2.5
  )
} else if (example == "local__multivariate") {
  ## Local - Multivariate ----

  n1 <- n2 <- 500
  signal <- 1
  xmarg <- "TN_U"
  grid_size <- 21

  name_file <- glue(
    "xmarg_{xmarg}",
    "n_{n1}",
    "signal_{signal}",
    "grid_{grid_size}",
    .sep = "__"
  )
  dir_result <- file.path("output", "results", "simulations", "local", "multivariate")
  load(file.path(dir_result, glue("{name_file}.RData")))

  figs <- vector("list", length(methods))
  names(figs) <- methods
  for (method in methods) {
    figs[[method]] <- data.frame(Xgrid, reject = reject[, method]) %>%
      mutate(reject = factor(reject)) %>%
      ggplot() +
      geom_tile(
        aes(x = X1, y = X2, fill = reject),
        color = "lightgrey", lwd = 0.01, linetype = 1
      ) +
      geom_hline(yintercept = 0, linetype = 1, linewidth = 1) +
      geom_vline(xintercept = 0, linetype = 1, linewidth = 1) +
      scale_fill_discrete(guide = "none") +
      coord_fixed() +
      my_theme_base() +
      ggtitle(method) +
      theme(plot.title = element_text(hjust = 0.5, size = 20))
  }
  fig <- grid.arrange(figs$CED, figs$CMMD, nrow = 1)

  dir_figure <- file.path("output", "figures")
  if (!dir.exists(dir_figure)) {
    dir.create(dir_figure, recursive = TRUE)
  }
  ggsave(
    file.path(dir_figure, glue("reject__local__multivariate.pdf")),
    plot = fig, width = 10, height = 5
  )
} else if (example == "global__covariate_shift__univariate") {
  ## Global - Covariate Shift - Univariate ----
  method_colors <- c("CED" = "#F8766D", "CMMD" = "#00AAFF", "CONF" = "#7CAE00", "ECF" = "#00BFC4")
  method_linetypes <- c("CED" = "solid", "CMMD" = "dashed", "CONF" = "dotted", "ECF" = "longdash")

  ### different mean functions ----
  params_fixed <- list(
    problem = "global",
    setting_shift = "covariate_shift",
    setting_dgp = "different_mean",
    methods = c("CED", "CMMD", "CONF", "ECF"),
    p = 1,
    shift_x = "true",
    dist_x = "normt",
    dist_eps = "norm",
    num_simu = 1000
  )

  params_varying <- list(
    signal_m = seq(0, 2, 0.4),
    setting_m = c("x2+x(x+2)(x-2)", "exp+sin"),
    n = c(50, 100)
  )

  power_array_diffm <- extract_result(
    "power_array", params_varying, params_fixed
  )
  dimnames(power_array_diffm)[["method"]] <- c("CED", "CMMD", "CONF", "ECF")
  dimnames(power_array_diffm)[["setting_m"]] <- c("Setting 3.1", "Setting 3.2")
  dimnames(power_array_diffm)[["n"]] <- c("n1 = n2 = 50", "n1 = n2 = 100")


  #### power table ----
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Power", gridLines = TRUE)

  row_cur <- 1
  for (n in dimnames(power_array_diffm)[["n"]]) {
    col_cur <- 1
    writeData(wb, 1, x = paste0("n = ", n), startCol = col_cur, startRow = row_cur)
    row_cur <- row_cur + 1

    for (setting_m in dimnames(power_array_diffm)[["setting_m"]]) {
      col_cur <- 2
      writeData(wb, 1, x = paste0("setting_m = ", setting_m), startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1

      power_block <- power_array_diffm[, , setting_m, n]
      writeData(wb, 1, x = power_block, rowNames = TRUE, startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1 + nrow(power_block)
    }
  }

  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveWorkbook(
    wb,
    file.path(dir, glue("power__{example}__diffm.xlsx")),
    overwrite = TRUE
  )

  #### power plot ----
  fig_diffm <- melt(
    power_array_diffm,
    varnames = c("c", "Methods", "Setting", "n"),
    value.name = "Power"
  ) %>%
    mutate(n = factor(n), type = "Difference in conditional means") %>%
    ggplot(
      aes(
        x = c, y = Power, color = Methods,
        linetype = Methods, shape = Methods
      )
    ) +
    geom_line() +
    geom_point() +
    facet_nested(
      n ~ type + Setting,
      nest_line = element_line(linetype = 1)
    ) +
    scale_shape_discrete("Method", solid = FALSE) +
    scale_linetype_manual("Method", values = method_linetypes) +
    scale_color_manual("Method", values = method_colors) +
    scale_x_continuous(
      breaks = seq(0, 2, by = 0.4),
      minor_breaks = seq(0, 2, by = 0.4)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1)
    ) +
    my_theme_base()


  ### different variance functions ----
  params_varying <- list(
    signal_v = seq(0, 1, 0.2),
    setting_v = c("homo", "hetero"),
    n = c(100, 200)
  )

  params_fixed <- list(
    problem = "global",
    setting_shift = "covariate_shift",
    setting_dgp = "different_var",
    methods = c("CED", "CMMD", "CONF", "ECF"),
    p = 1,
    shift_x = "true",
    # dist_x = "norm",
    dist_x = "normt",
    dist_eps = "norm",
    setting_m = "quad",
    num_simu = 1000
  )

  power_array_diffv <- extract_result(
    "power_array", params_varying, params_fixed
  )
  dimnames(power_array_diffv)[["method"]] <- c("CED", "CMMD", "CONF", "ECF")
  dimnames(power_array_diffv)[["setting_v"]] <- c("Setting 3.3", "Setting 3.4")
  dimnames(power_array_diffv)[["n"]] <- c("n1 = n2 = 100", "n1 = n2 = 200")

  #### power table ----
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Power", gridLines = TRUE)

  row_cur <- 1
  for (n in dimnames(power_array_diffv)[["n"]]) {
    col_cur <- 1
    writeData(wb, 1, x = paste0("n = ", n), startCol = col_cur, startRow = row_cur)
    row_cur <- row_cur + 1

    for (setting_v in dimnames(power_array_diffv)[["setting_v"]]) {
      col_cur <- 2
      writeData(wb, 1, x = paste0("setting_v = ", setting_v), startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1

      power_block <- power_array_diffv[, , setting_v, n]
      writeData(wb, 1, x = power_block, rowNames = TRUE, startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1 + nrow(power_block)
    }
  }

  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveWorkbook(
    wb,
    file.path(dir, glue("power__{example}__diffv.xlsx")),
    overwrite = TRUE
  )

  #### power plot ----
  fig_diffv <- melt(
    power_array_diffv,
    varnames = c("c", "Methods", "Setting", "n"),
    value.name = "Power"
  ) %>%
    mutate(n = factor(n), type = "Difference in conditional variances") %>%
    ggplot(
      aes(
        x = c, y = Power, color = Methods,
        linetype = Methods, shape = Methods)
    ) +
    geom_line() +
    geom_point() +
    facet_nested(
      n ~ type + Setting,
      nest_line = element_line(linetype = 1)
    ) +
    scale_shape_discrete("Method", solid = FALSE) +
    scale_linetype_manual("Method", values = method_linetypes) +
    scale_color_manual("Method", values = method_colors) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.2)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1)
    ) +
    my_theme_base ()

  ### Combine two plots ----
  dir <- file.path("output", "figures")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  g_legend <- function(p) {
    g <- ggplotGrob(p)
    guide_grob <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(guide_grob)
  }
  legend <- g_legend(fig_diffm)
  fig_diffm <- fig_diffm + theme(legend.position = "none")
  fig_diffv <- fig_diffv + theme(legend.position = "none")

  combined <- grid.arrange(
    arrangeGrob(fig_diffm, fig_diffv, ncol = 2),
    legend,
    ncol = 1,
    heights = c(10, 1) # Adjust heights to control space for the legend
  )

  ggsave(
    file.path(dir, glue("power__{example}.pdf")),
    combined,
    width = 11, height = 4.5
  )

} else if (example == "global__covariate_shift__multivariate") {
  ## Global - Covariate Shift - Multivariate  ----
  params_fixed <- list(
    problem = "global",
    setting_shift = "covariate_shift",
    setting_dgp = "different_mean",
    methods = c("CED", "CMMD", "CONF"),
    p = 4,
    shift_x = "true",
    dist_x = "norm",
    dist_eps = "t",
    num_simu = 1000
  )

  params_varying <- list(
    signal_m = seq(0, 2, 0.4),
    setting_m = c("sindex_quad"),
    n = c(100, 200)
  )

  power_array <- extract_result(
    "power_array", params_varying, params_fixed
  )
  dimnames(power_array)[["method"]] <- c("CED", "CMMD", "CONF")
  dimnames(power_array)[["setting_m"]] <- c("Setting 4.1")
  dimnames(power_array)[["n"]] <- c("n1 = n2 = 100", "n1 = n2 = 200")

  #### power table ----
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Power", gridLines = TRUE)

  row_cur <- 1
  for (n in dimnames(power_array)[["n"]]) {
    col_cur <- 1
    writeData(wb, 1, x = paste0("n = ", n), startCol = col_cur, startRow = row_cur)
    row_cur <- row_cur + 1

    for (setting_m in dimnames(power_array)[["setting_m"]]) {
      col_cur <- 2
      writeData(wb, 1, x = paste0("setting_m = ", setting_m), startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1

      power_block <- power_array[, , setting_m, n]
      writeData(wb, 1, x = power_block, rowNames = TRUE, startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1 + nrow(power_block)
    }
  }

  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveWorkbook(
    wb,
    file.path(dir, glue("power__{example}.xlsx")),
    overwrite = TRUE
  )

  ### power plot ----
  method_colors <- c("CED" = "#F8766D", "CMMD" = "#00AAFF", "CONF" = "#7CAE00")

  fig <- melt(
    power_array,
    varnames = c("c", "Methods", "Setting", "n"),
    value.name = "Power"
  ) %>%
    mutate(n = factor(n)) %>%
    ggplot(
      aes(
        x = c, y = Power, color = Methods,
        linetype = Methods, shape = Methods
      )
    ) +
    geom_line() +
    geom_point() +
    facet_wrap(~ n) +
    scale_shape_discrete("Method", solid = FALSE) +
    scale_linetype_discrete("Method") +
    scale_color_manual("Method", values = method_colors) +
    scale_x_continuous(
      breaks = seq(0, 2, by = 0.4),
      minor_breaks = seq(0, 2, 0.4)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1)
    ) +
    my_theme_base()

  dir <- file.path("output", "figures")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  ggsave(
    file.path(dir, glue("power__{example}.pdf")),
    # width = 6, height = 4
    width = 8, height = 4
  )

} else if (example == "global__label_shift") {
  ## Global - Label Shift ----
  params_fixed <- list(
    problem = "global",
    setting_shift = "label_shift",
    setting_dgp = "different_mean",
    methods = c("CED", "CMMD", "CONF"),
    shift_y = "true",
    num_simu = 1000
  )

  params_varying <- list(
    signal_m = seq(0, 1, 0.2),
    setting_m = c("nonlinear"),
    n = c(100, 200),
    p = c(5, 20)
  )

  power_array <- extract_result(
    "power_array", params_varying, params_fixed
  )
  dimnames(power_array)[["method"]] <- c("CED", "CMMD", "CONF")
  dimnames(power_array)[["setting_m"]] <- c("Setting (5.1)")
  dimnames(power_array)[["n"]] <- c("n1 = n2 = 100", "n1 = n2 = 200")
  dimnames(power_array)[["p"]] <- c("p = 5", "p = 20")

  ### power table ----
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Power", gridLines = TRUE)

  row_cur <- 1
  for (p in dimnames(power_array)[["p"]]) {
    col_cur <- 1
    writeData(wb, 1, x = paste0("p = ", p), startCol = col_cur, startRow = row_cur)
    row_cur <- row_cur + 1

    for (n in dimnames(power_array)[["n"]]) {
      col_cur <- 2
      writeData(wb, 1, x = paste0("n = ", n), startCol = col_cur, startRow = row_cur)
      row_cur <- row_cur + 1

      for (setting_m in dimnames(power_array)[["setting_m"]]) {
        col_cur <- 3
        writeData(wb, 1, x = paste0("setting_m = ", setting_m), startCol = col_cur, startRow = row_cur)
        row_cur <- row_cur + 1

        power_block <- power_array[, , setting_m, n, p]
        writeData(wb, 1, x = power_block, rowNames = TRUE, startCol = col_cur, startRow = row_cur)
        row_cur <- row_cur + 1 + nrow(power_block)
      }
    }
  }

  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveWorkbook(
    wb,
    file.path(dir, glue("power__{example}.xlsx")),
    overwrite = TRUE
  )

  ### power plot ----
  method_colors <- c("CED" = "#F8766D", "CMMD" = "#00AAFF", "CONF" = "#7CAE00")

  fig <- melt(
    power_array,
    varnames = c("c", "Methods", "Setting", "n", "p"),
    value.name = "Power"
  ) %>%
    mutate(n = factor(n)) %>%
    ggplot(
      aes(
        x = c, y = Power, color = Methods,
        linetype = Methods, shape = Methods
      )
    ) +
    geom_line() +
    geom_point() +
    facet_grid(p ~ n) +
    scale_shape_discrete("Method", solid = FALSE) +
    scale_linetype_discrete("Method") +
    scale_color_manual("Method", values = method_colors) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, 0.2)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      minor_breaks = seq(0, 1, by = 0.1)
    ) +
    my_theme_base() +
    theme(legend.position = "right")

  dir <- file.path("output", "figures")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  ggsave(
    file.path(dir, glue("power__{example}.pdf")),
    width = 8, height = 4
  )
} else if (example == "realdata__airfoil") {
  # Realdata ----
  ## Airfoil (global) ----
  methods <- c("CED", "CMMD")
  settings_null <- c("cov_shift_null", "prior_shift_null")
  settings_alt <- c("cov_shift_alt", "prior_shift_alt")

  results_reject_rate <- matrix(
    nrow = 2, ncol = 2,
    dimnames = list(settings_null, methods)
  )
  results_pvalue <- matrix(
    nrow = 2, ncol = 2,
    dimnames = list(settings_alt, methods)
  )

  for (setting in settings_null) {
    obj_result <- env()
    load(
      file.path(
        "output", "results", "realdata",
        glue("realdata__airfoil__{setting}.RData")
      ),
      obj_result
    )
    results_reject_rate[setting, ] <- colMeans(obj_result$reject)[methods]
  }
  rownames(results_reject_rate) <- paste("Setting", 1:2)

  for (setting in settings_alt) {
    obj_result <- env()
    load(
      file.path(
        "output", "results", "realdata",
        glue("realdata__airfoil__{setting}.RData")
      ),
      obj_result
    )
    results_pvalue[setting, ] <- obj_result$pvalue[1, methods]
  }
  rownames(results_pvalue) <- paste("Setting", 3:4)


  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Power", gridLines = TRUE)

  col_cur <- 1
  writeData(wb, 1, x = "Rejection rate", startCol = col_cur, startRow = 1)
  writeData(wb, 1, x = results_reject_rate, rowNames = TRUE, startCol = col_cur + 1, startRow = 2)

  col_cur <- col_cur + ncol(results_reject_rate) + 3
  writeData(wb, 1, x = "p-value", startCol = col_cur, startRow = 1)
  writeData(wb, 1, x = results_pvalue, rowNames = TRUE, startCol = col_cur + 1, startRow = 2)

  dir <- file.path("output", "tables")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveWorkbook(
    wb,
    file.path(dir, glue("{example}.xlsx")),
    overwrite = TRUE
  )

  ## Airfoil (local) ----
  path <- file.path("output", "results", "realdata")
  load(file.path(path, glue("realdata__airfoil__prior_shift_alt__local.RData")))

  df_pvalue <- rbind(
    data.frame(
      x = x0_all,
      pvalue = pvalue[, "CMMD"],
      method = "CMMD"
    ),
    data.frame(
      x = x0_all,
      pvalue = pvalue[, "CED"],
      method = "CED"
    )
  )

  data %>%
    ggplot() +
    geom_line(aes(x = x, y = pvalue, color = method), data = df_pvalue) +
    my_theme_base() +
    scale_color_manual("Method", values = c("CED" = "#F8766D", "CMMD" = "#00AAFF"), breaks = c("CED", "CMMD")) +
    xlab("Local y (scaled sound pressure level)") +
    ylab("p-value") +
    ylim(c(0, 0.15)) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
    ) +
    geom_hline(yintercept = 0.05, linetype = 2)


  dir_fig <- file.path("output", "figures")
  if (!dir.exists(dir_fig)) {
    dir.create(dir_fig, recursive = TRUE)
  }
  ggsave(
    file.path(
      dir_fig,
      glue("airfoil_data_local.pdf")
    ),
    width = 8, height = 3
  )

}
