proliferation_location = function() {
    "/oak/stanford/groups/pritch/users/jweinstk/resources/proliferation_legut/table_s4_legut_t_cell_proliferation.csv"
}

read_proliferation = function() {
    vroom::vroom(proliferation_location(), col_select = c("gene_symbol", "stat", "padj", "log2FoldChange")) %>%
        dplyr::mutate(
            proliferation = padj < 5e-2
        ) %>%
        dplyr::select(
            gene_name = gene_symbol,
            proliferation,
            proliferation_log2FC = log2FoldChange,
            zscore = stat
        )
}

plot_proliferation = function(meta, proliferation) {
    fc_counts = meta %>%
        dplyr::distinct(gene_name = KO, gene_FC_counts) %>%
        tidyr::drop_na(.)

    dfm = proliferation %>%
        dplyr::inner_join(fc_counts, by = "gene_name")

    plot = ggplot(
            data = dfm, 
            aes(x = proliferation_log2FC, y = gene_FC_counts, color = proliferation, label = gene_name)) +
        geom_point() +
        # geom_smooth() + 
        cowplot::theme_cowplot() +
        ggrepel::geom_label_repel() +
        labs(x = "Proliferative effect", y = "Normalized cell counts", color = "Proliferative\neffect pvalue < .05")

    ggsave(
        file.path(figure_dir(), "proliferation_scatter.pdf"), plot,
        width = 6, height = 4, units = "in"
    )
}

plot_proliferation_ko_specific = function(ko_specific_effects, proliferation, meta, tag = "diffeq") {

    counts = purrr::map_dbl(ko_specific_effects, ~{
            sig_proliferation = proliferation$gene_name[proliferation$proliferation]
            sum(.x %in% sig_proliferation) 
        }) %>%
        tibble::enframe(.) %>%
        dplyr::arrange(desc(value)) %>%
        dplyr::inner_join(meta %>% dplyr::distinct(KO, gene_group), by = c("name" = "KO"))

    print(summary(MASS::glm.nb(value ~ gene_group, data = counts)))

    top_counts = counts %>%
        dplyr::slice(1:20) 

    normalized_counts = purrr::map_dbl(ko_specific_effects, ~{
            sig_proliferation = proliferation$gene_name[proliferation$proliferation]
            sum(.x %in% sig_proliferation) / length(.x)
        }) %>%
        tibble::enframe(.) %>%
        dplyr::arrange(desc(value)) %>%
        dplyr::inner_join(meta %>% dplyr::distinct(KO, gene_group), by = c("name" = "KO"))

    top_normalized_counts = normalized_counts %>%
        dplyr::slice(1:20) 

    print(summary(MASS::glm.nb(value ~ gene_group, data = normalized_counts)))

    fnames = c(
        "top" = file.path(figure_dir(), glue::glue("top_ranked_effects_on_proliferation_{tag}.pdf")),
        "top_normalized" = file.path(figure_dir(), glue::glue("top_normalized_ranked_effects_on_proliferation_{tag}.pdf")),
        "distribution" = file.path(figure_dir(), glue::glue("ridgeline_effects_on_proliferation_{tag}.pdf")),
        "distribution_normalized" = file.path(figure_dir(), glue::glue("normalized_ridgeline_effects_on_proliferation_{tag}.pdf"))
    )

    plot = ggplot(data = top_counts, aes(x = forcats::fct_reorder(name, value, .desc = TRUE), y = value, fill = gene_group)) +
            geom_col() +
            labs(y = "Number of perturbation specific effects on proliferation", fill = "Gene group") +
            scale_y_continuous(expand = expand_scale(mult = c(0.0, 0.0))) +
            add_gene_group_fill_ggplot2() +
            cowplot::theme_cowplot(font_size = 12) +
            theme(
                legend.position = "top",
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle=45, vjust=1, hjust=1)
            ) 
            # guides(fill = "none")

    ggsave(fnames["top"], plot, width = 8, height = 6, units = "in")

    plot = ggplot(data = top_normalized_counts, aes(x = forcats::fct_reorder(name, value, .desc = TRUE), y = value, fill = gene_group)) +
            geom_col() +
            labs(y = "Number of perturbation specific effects on proliferation", fill = "Gene group") +
            scale_y_continuous(expand = expand_scale(mult = c(0.0, 0.0))) +
            add_gene_group_fill_ggplot2() +
            cowplot::theme_cowplot(font_size = 12) +
            theme(
                legend.position = "top",
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle=45, vjust=1, hjust=1)
            ) 
            # guides(fill = "none")

    ggsave(fnames["top_normalized"], plot, width = 8, height = 6, units = "in")


    # plot = ggplot(data = counts, aes(y = gene_group, x = value, fill = factor(stat(quantile)))) +
    plot = ggplot(data = counts, aes(y = gene_group, x = value)) +
                ggridges::stat_density_ridges(
                        geom = "density_ridges_gradient", calc_ecdf = TRUE
                        # quantiles = 4, quantile_lines = TRUE
                ) +
                # scale_fill_viridis_d(name = "Quartiles") +
                scale_x_continuous(trans = "log1p") + 
                cowplot::theme_cowplot(font_size = 12) +
                cowplot::panel_border() +
                theme(axis.title.y = element_blank(), axis.title.x = element_blank())


    ggsave(fnames["distribution"], plot, width = 8, height = 6, units = "in")

    # plot = ggplot(data = normalized_counts, aes(y = gene_group, x = value, fill = factor(stat(quantile)))) +
    plot = ggplot(data = normalized_counts, aes(y = gene_group, x = value)) +
                ggridges::stat_density_ridges(
                        geom = "density_ridges_gradient", calc_ecdf = TRUE
                        # quantiles = 4, quantile_lines = TRUE
                ) +
                # scale_fill_viridis_d(name = "Quartiles") +
                # scale_x_continuous(trans = "log1p") + 
                cowplot::theme_cowplot(font_size = 12) +
                cowplot::panel_border() +
                theme(axis.title.y = element_blank(), axis.title.x = element_blank())


    ggsave(fnames["distribution_normalized"], plot, width = 8, height = 6, units = "in")
}


plot_cytokine_proliferation_effect = function(proliferation, cytokine_hits, mashr, meta, lfsr_threshold = set_lfsr_threshold()) {

    pm = ashr::get_pm(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "beta", -to) 

    lfsr = ashr::get_lfsr(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "lfsr", -to) %>%
        dplyr::inner_join(pm) %>%
        dplyr::inner_join(
            meta %>% 
                add_gene_group_colors %>%
                dplyr::distinct(KO, gene_group, color), 
            by = c("from" = "KO")
        )

    print(lfsr)

    edges = lfsr %>%
        dplyr::filter(lfsr < lfsr_threshold) %>%
        dplyr::inner_join(
            proliferation %>%
                dplyr::select(
                    gene_name, 
                    proliferation_zscore = zscore
                ),
            by = c("to" = "gene_name")
        ) %>%
        dplyr::inner_join(
            cytokine_hits %>%
                dplyr::filter(Cytokine == "IFNG") %>%
                dplyr::select(
                    gene_name,
                    IFNg_zscore = zscore
                ),
            by = c("to" = "gene_name")
        )

    summarized_edges = edges %>%
        dplyr::group_by(from, gene_group, color) %>%
        dplyr::summarize(
            proliferation = sum(proliferation_zscore * sign(beta), na.rm = T),
            IFNg = sum(IFNg_zscore * sign(beta), na.rm = T),
        )

    plot = ggplot(summarized_edges, aes(x = proliferation, y = IFNg, color = gene_group, label = from)) +
        geom_point() +
        add_gene_group_color_ggplot2() +
        ggrepel::geom_label_repel() +
        cowplot::theme_cowplot(font_size = 11) +
        labs(
            x = "Aggregate downstream effect on proliferation",
            y = "Aggregate downstream effect on IFNg secretion"
        )

    ggsave(
        file.path(figure_dir(), "joint_downstream_proliferation_IFNg_scatter.pdf"),
        plot,
        width = 7,
        height = 5,
        units = "in"
    )

    return(edges)
}
