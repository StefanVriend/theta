plot_theta_mult <- function(filename,
                            model_list,
                            pars = c("sigma_e2", "K", "mu_r1", "theta", "gamma"),
                            true = c(sigma_e2, K, mu_r1, theta, gamma)) {
  # filename: character. Name of plot file without extension.
  # model_list: named list of models.
  # pars: character. Parameters to plot.
  # true: numeric. True values of {pars}

  # Create temporary folder to save individual pdfs in
  # NB: ggsave + cairo_pdf do not allow multi-page pdfs
  dir.create("temp-figures")

  # Table of possible models
  model_markers <- tibble::tibble(
    model_name = c("Approx", "MLE", "RE-r", "RE-s", "RE-r-hyp", "RE-s-hyp"),
    model_col = c("#2c91a0", "#ff9f1c", "#B54057", "#99b12f", "#007a5e", "#40476d")
  )

  # Select names and colors
  models <- model_markers %>%
    dplyr::filter(model_name %in% names(model_list))

  model_names <- models %>% dplyr::pull(model_name)
  cl <- models %>% dplyr::pull(model_col)
  trace_cl <- purrr::map(.x = cl,
                         .f = ~{c(colorspace::darken(.x, amount = 0.5), .x,
                                  colorspace::lighten(.x, amount = 0.5))}
  ) %>%
    unlist()

  # Create pdfs
  purrr::walk(.x = seq_len(pops),
              .f = ~{

                cat(paste0("Population: ", .x, "\n"))

                data <- tibble::tibble(
                  N = obs_N[.x, ],
                  t = seq_len(tmax)
                ) %>%
                  dplyr::mutate(r = log(lead(N)) - log(N))

                # Plot theme
                theme_set(theme_classic(base_family = "Roboto"))
                theme_update(
                  axis.text = element_text(size = 11, color = "black"),
                  axis.title = element_text(size = 12, color = "black"),
                  plot.margin = margin(l = 5, b = 5, t = 10, r = 15),
                  legend.position = "none"
                )

                # Plot 1: N versus time
                p1 <- ggplot(data, aes(x = t, y = N)) +
                  geom_point(size = 2, shape = 21) +
                  geom_line() +
                  scale_y_continuous(breaks = scales::pretty_breaks()) +
                  labs(x = "Time", y = "N")

                # Plot 2: r versus N
                # Retrieve posterior samples for N, and r or s
                # post_N <- as.matrix(mod)[,grepl(paste0("N\\[", .x, ",.*"),
                #                                 colnames(as.matrix(mod)))][,first[.x]:(last[.x]-1)]
                #
                # if(any(grepl("r\\[", colnames(as.matrix(mod))))) {
                #
                #   post_growth <- as.matrix(mod)[,grepl(paste0("r\\[", .x, ",.*"),
                #                                        colnames(as.matrix(mod)))][,first[.x]:(last[.x]-1)]
                #
                # } else {
                #
                #   post_growth <- as.matrix(mod)[,grepl(paste0("s\\[", .x, ",.*"),
                #                                        colnames(as.matrix(mod)))][,first[.x]:(last[.x]-1)]
                #
                # }
                #
                # post <- tibble::tibble(
                #   N_low = matrixStats::colQuantiles(post_N, probs = 0.025),
                #   N_med = matrixStats::colQuantiles(post_N, probs = 0.5),
                #   N_high = matrixStats::colQuantiles(post_N, probs = 0.975),
                #   N_mode = apply(post_N, 2, max_den),
                #   gr_low = matrixStats::colQuantiles(post_growth, probs = 0.025),
                #   gr_med = matrixStats::colQuantiles(post_growth, probs = 0.5),
                #   gr_high = matrixStats::colQuantiles(post_growth, probs = 0.975),
                #   gr_mode = apply(post_growth, 2, max_den)
                # )

                p2 <- ggplot() +
                  geom_point(mapping = aes(x = N, y = r), data = data, size = 2, shape = 21) +
                  #geom_point(mapping = aes(x = N_mode, y = gr_mode), data = post, size = 2, shape = 21, col = "red") +
                  scale_x_continuous(breaks = scales::pretty_breaks()) +
                  scale_y_continuous(breaks = scales::pretty_breaks()) +
                  labs(x = "N", y = bquote(log(N[t+1]~"/"~N[t])))

                # Plots 3: Posterior distributions of parameters
                plot_posterior <- function(par, true) {

                  post_data <- purrr::map2_dfr(.x = model_list,
                                               .y = model_names,
                                               .f = ~{

                                                 tibble::tibble(
                                                   y = as.matrix(.x)[, par],
                                                   group = as.character(.y)
                                                 )

                                               }) %>%
                    dplyr::mutate(group = forcats::fct_relevel(group, rev(model_names)))

                  mean_data <- post_data %>%
                    dplyr::group_by(group) %>%
                    dplyr::summarise(median = quantile(y, probs = 0.5, na.rm = TRUE),
                                     low = quantile(y, probs = 0.025, na.rm = TRUE),
                                     high = quantile(y, probs = 0.975, na.rm = TRUE),
                                     .groups = "drop")

                  if(stringr::str_detect(par, "theta")) {
                    # MLE has no theta bootstrap info

                    mle <- which(model_names == "MLE")
                    model_names <- model_names[-mle]
                    cl <- cl[-mle]

                    post_data <- post_data %>%
                      dplyr::filter(group != "MLE") %>%
                      dplyr::mutate(group = forcats::fct_drop(group))

                    mean_data <- mean_data %>%
                      dplyr::filter(group != "MLE") %>%
                      dplyr::mutate(group = forcats::fct_drop(group))

                  }

                  ggplot(data = post_data, aes(y = group, x = y)) +
                    geom_vline(xintercept = true, linetype = "dashed") +
                    ggridges::geom_density_ridges(aes(fill = group, color = group), scale = 0.85, alpha = 0.5) +
                    geom_segment(aes(y = group, yend = group, x = low, xend = high, color = group), data = mean_data, size = 1) +
                    geom_point(aes(y = group, x = median, color = group), data = mean_data, size = 3) +
                    labs(x = par, y = "") +
                    scale_color_manual(values = rev(cl)) +
                    scale_fill_manual(values = colorspace::lighten(rev(cl), amount = 0.35))

                }

                p3 <- purrr::map2(.x = paste0(pars, "[", .x, "]"),
                                  .y = true,
                                  .f = ~plot_posterior(.x, .y))

                # Plots 4: Trace plots of parameters
                plot_trace <- function(par) {

                  trace_data <- purrr::map2_dfr(.x = model_list,
                                                .y = model_names,
                                                .f = ~{

                                                  if(.y != "MLE") {

                                                    tibble::tibble(
                                                      y = c(.x$chain1[, par], .x$chain2[, par], .x$chain3[, par]),
                                                      sample = rep(1:nrow(.x$chain1), 3),
                                                      chain = as.character(rep(c(1, 2, 3), each = nrow(.x$chain1))),
                                                      group = as.character(.y),
                                                      gr_ch = paste(group, chain, sep = "-")
                                                    )

                                                  } else {

                                                    tibble::tibble(
                                                      y = .x %>% dplyr::pull(par),
                                                      sample = rep(1:(nrow(.x) / 3), 3),
                                                      chain = as.character(rep(c(1, 2, 3), each = nrow(.x) / 3)),
                                                      group = as.character(.y),
                                                      gr_ch = paste(group, chain, sep = "-")
                                                    )

                                                  }

                                                }) %>%
                    dplyr::mutate(group = forcats::fct_relevel(group, model_names),
                                  gr_ch = forcats::fct_relevel(gr_ch, purrr::map(model_names, ~paste0(.x, "-", 1:3)) %>% unlist()))


                  ggplot(data = trace_data, aes(x = sample, y = y, color = gr_ch)) +
                    geom_line(size = 1/3, alpha = 0.8) +
                    facet_wrap(~group, ncol = length(model_names)) +
                    labs(x = "", y = par) +
                    scale_color_manual(values = trace_cl) +
                    theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

                }

                p4 <- purrr::map(.x = paste0(pars, "[", .x, "]"),
                                 .f = ~plot_trace(.x))

                # Title
                title <- cowplot::ggdraw() +
                  cowplot::draw_label(paste0("Population: ", .x), fontface = "bold")

                # Combine plots
                cowplot::plot_grid(
                  title,
                  cowplot::plot_grid(
                    cowplot::plot_grid(p1, p2, NULL, NULL, nrow = 4, ncol = 1, align = "v"),
                    cowplot::plot_grid(plotlist = p3, nrow = length(pars), ncol = 1, align = "v"),
                    cowplot::plot_grid(plotlist = p4, nrow = length(pars), ncol = 1, align = "v"),
                    nrow = 1, ncol = 3, rel_widths = c(1.2, 1.2, 3)
                  ),
                  nrow = 2, rel_heights = c(0.075, 1)
                )

                ggsave(filename = paste0("temp-figures/temp", .x, ".pdf"),
                       device = cairo_pdf,
                       plot = ggplot2::last_plot(),
                       width = 16, height = 12)

              })

  # Combine individual pdfs
  pdftools::pdf_combine(paste0("temp-figures/temp", seq_len(pops), ".pdf"),
                        output = paste0(filename, ".pdf"))

  # Remove temporary folder
  unlink("temp-figures", recursive = TRUE)

}
