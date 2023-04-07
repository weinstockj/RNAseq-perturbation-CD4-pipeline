group_colors_ = function() {
    colors = c(
        "Control" = "#d95f02",
        "IEI Target" = "#1b9e77",
        "IL2RA Regulators" = "#7570b3"
    )

    return(colors)
}

group_colors_extra_ = function() {
    colors = c(group_colors_(), "not-KO'd" = "black")
    return(colors)
}

get_gene_group_color = function(group) {

    colors = group_colors_()
    return(colors[group])
}

get_gene_group_color_extra = function(group) {

    colors = group_colors_extra_()
    return(colors[group])
}

add_gene_group_colors = function(meta) {
    meta %>%
        dplyr::mutate(
            color = get_gene_group_color(gene_group)
        ) 
}

add_gene_group_colors_extra = function(meta) {
    meta %>%
        dplyr::mutate(
            color = get_gene_group_color_extra(gene_group)
        ) 
}

add_gene_group_color_ggplot2 = function() {
    
    colors = group_colors_()

    return(scale_color_manual(values = colors))
}

add_gene_group_color_extra_ggplot2 = function() {
    
    colors = group_colors_extra_()

    return(scale_color_manual(values = colors))
}

add_gene_group_fill_ggplot2 = function() {
    
    colors = group_colors_()

    return(scale_fill_manual(values = colors))
}

add_gene_group_fill_extra_ggplot2 = function() {
    
    colors = group_colors_extra_()

    return(scale_fill_manual(values = colors))
}

create_tidy_graph = function(parsed_chain, meta) {
    
    # targets = meta %>%
    #     dplyr::filter(!is_control) %>%
    #     dplyr::pull(KO) %>%
    #     as.character %>%
    #     unique

    targets = unique(c(parsed_chain$row, parsed_chain$col))

    status = meta %>% 
        dplyr::mutate(KO = as.character(KO)) %>%
        dplyr::distinct(KO, gene_group) %>%
        add_gene_group_colors %>%
        dplyr::rename(name = KO)

    graph = tidygraph::tbl_graph(
        nodes = tibble::tibble(name = targets) %>% dplyr::inner_join(status),
        edges = parsed_chain %>% 
            dplyr::select(from = row, to = col, signed_weight = estimate) %>%
            dplyr::mutate(weight = abs(signed_weight))
    )

    return(graph)
}

create_total_tidy_graph = function(
    parsed_chain, mashr, expression,
    meta,
    constraint,
    iei_genes = iei_genes(),
    gwas_genes = pics_ai_genes(),
    trans_egenes,
    threshold = set_causal_threshold(), lfsr_threshold = set_lfsr_threshold()) {

    upstream_edges = parsed_chain %>% 
        dplyr::select(from = row, to = col, signed_weight = estimate) %>%
        dplyr::mutate(weight = abs(signed_weight)) %>%
        dplyr::filter(abs(signed_weight) > threshold)

    pm = ashr::get_pm(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "signed_weight", -to) %>%
        dplyr::mutate(weight = abs(signed_weight))

    lfsr = ashr::get_lfsr(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "lfsr", -to) %>%
        dplyr::inner_join(pm) %>%
        dplyr::inner_join(
            meta %>% dplyr::distinct(KO, gene_group), 
            by = c("from" = "KO")
        )

    edges = lfsr %>%
        dplyr::filter(lfsr < lfsr_threshold) %>%
        dplyr::select(from, to, signed_weight, weight) %>%
        dplyr::bind_rows(upstream_edges)

    targets = unique(c(edges$from, edges$to))

    status = meta %>% 
                dplyr::mutate(KO = as.character(KO)) %>%
                dplyr::distinct(KO, gene_group) %>%
                dplyr::bind_rows(
                    tibble::tibble(
                        KO = lfsr$to,
                        gene_group = rep("not-KO'd", nrow(lfsr))
                    ) %>%
                    dplyr::distinct(.)
                ) %>%
                add_gene_group_colors_extra %>%
                dplyr::rename(name = KO) %>%
                dplyr::left_join(constraint, by = c("name" = "gene_name")) %>%
                dplyr::left_join(expression, by = c("name" = "gene_name")) %>%
                dplyr::mutate(
                    expression = log1p(expression),
                    is_iei = name %in% iei_genes,
                    is_gwas = name %in% gwas_genes,
                    is_trans_egene = name %in% trans_egenes
                )

    graph = tidygraph::tbl_graph(
        nodes = tibble::tibble(name = targets) %>% 
            dplyr::inner_join(status, by = "name"),

        edges = edges %>% 
            dplyr::select(from, to, signed_weight, weight)  
    )

    return(graph)
    
}

plot_network = function(parsed_chain, meta, threshold = set_causal_threshold(), results = NULL, txdb = NULL, output_dir = figure_dir(), tag = "", edges = "straight", layout = "stress", scale_node_size_by_degree = TRUE, width = 13.5, height = 10, label_size = 6.3, ...) {


    stopifnot(xor(is.null(threshold), is.null(results)))
    
    if(!is.null(threshold)) {

        reduced_graph = parsed_chain %>% 
                dplyr::filter(abs(estimate) > threshold)

    } else {
        reduced_graph = parsed_chain %>%
            filter_on_sig_deg(results, txdb) 

        threshold = "sig_dge"

    }

    logger::log_info(glue::glue("plotting {nrow(reduced_graph)} edges"))
    
    graph = create_tidy_graph(
        reduced_graph,
        meta
    ) %>%
        tidygraph::mutate(degree = tidygraph::centrality_degree())

    # plot = ggraph(graph, "stress", bbox = 10) + 
    # plot = ggraph(graph, "auto", circular = TRUE) + 
    plot = ggraph(graph, layout, ...) 
    # plot = ggraph(graph, "centrality", cen = tidygraph::graph.strength(graph)) + 
                # geom_edge_arc2(

    if(edges == "straight") {
        plot = plot + 
                geom_edge_link(
                    aes(colour = signed_weight),
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        alpha = .25,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) 

    } else {
        plot = plot + 
                geom_edge_arc(
                    aes(colour = signed_weight),
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        strength = .2,
                        alpha = .5,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) 

    }

    plot = plot + 
            # geom_node_point(aes(colour = gene_group)) +
            scale_edge_colour_gradient2(low = "blue", mid = "gray", high = "red")  
            # geom_node_point(aes(size = 2.5 + 0.08 * log1p(degree), color = I(color)), alpha = .7) +

    if(scale_node_size_by_degree) {

        plot = plot +
            geom_node_point(aes(size = 3.0 + 0.08 * sqrt(1L + degree), color = I(color)), alpha = .7) 
    } else {
        plot = plot + 
            geom_node_point(aes(color = I(color)), alpha = .7) 
    }

    plot = plot + 
            geom_node_label(
                aes(label = name, color = I(color)),
                size = label_size,
                # colour = "black", 
                # family = "serif",
                check_overlap = TRUE,
                repel = TRUE
            ) + 
            labs(colour = "Control TF") +
            theme_graph(base_family = "Helvetica")

    fname = file.path(output_dir, glue::glue("ggplot_network_direct_effects_{threshold}_{tag}.pdf"))
    ggsave(
        fname,
        plot,
        # width = 13.5,
        width = width,
        # height = 10,
        height = height,
        units = "in"
    )

    return(fname)
}

filter_on_sig_deg = function(dfm, results, txdb) {
    stopifnot(all(c("row", "col") %in% names(dfm)))

    sig_results = results %>%
        dplyr::inner_join(txdb, by = c("gene" = "gene_id")) %>%
        dplyr::rename(row = KO, col = gene_name) %>%
        dplyr::filter(padj < .05)

    logger::log_info(glue::glue("Identified {nrow(sig_results)} significant total effects"))

    dfm %>%
        dplyr::inner_join(sig_results, by = c("row", "col"))
}

plot_indegree_outdegree_scatter = function(centrality) {

    plot = ggplot(data = centrality, aes(x = out_degree, y = in_degree, label = name, color = gene_group)) +
        geom_point() +
        add_gene_group_color_ggplot2() +
        ggrepel::geom_label_repel(max.overlaps = 5) +
        labs(x = "Outgoing connections", y = "Incoming connections", color = "Group") +
        cowplot::theme_cowplot(font_size = 12)


    ggsave(
        file.path(figure_dir(), "indegree_outdegree_scatter.pdf"),
        plot,
        width = 6,
        height = 4
    )
}

plot_direct_indirect = function(causal_network, total_effects, results, txdb, cyclic_genes, output_dir = figure_dir()) {

    fnames = c(
        "direct_indirect" = file.path(output_dir, "direct_indirect_scatter.pdf"),
        "direct_total" = file.path(output_dir, "direct_total_scatter.pdf")
    )

    dfm = total_effects %>%
        dplyr::inner_join(causal_network) %>%
        dplyr::filter(abs(estimate) > set_causal_threshold()) %>%
        mutate(
            label = glue::glue("{row}%->%{col}"),
            is_cyclic = (row %in% cyclic_genes) | (col %in% cyclic_genes),
            parsed_label = map_chr(label, ~as.character(as.expression(.x)))
        ) 
        

    plot = dfm %>%
        mutate(
            label = glue::glue("{row}%->%{col}"),
            is_cyclic = (row %in% cyclic_genes) | (col %in% cyclic_genes),
            indirect = total_effect - estimate,
            # parsed_label = map_chr(label, ~parse(text = as.expression(.x)))
            parsed_label = map_chr(label, ~as.character(as.expression(.x)))
        ) %>%
        ggplot(data = ., aes(x = estimate, y = indirect, label = parsed_label)) +
            cowplot::theme_cowplot(font_size = 12) +
            geom_smooth(method = "lm") +
            geom_point(aes(x = estimate, y = indirect, color = is_cyclic), alpha = .5) + 
            ggrepel::geom_text_repel(max.overlaps = 7, parse = TRUE) +
            labs(x = "Direct effect", y = "Indirect effect", color = "Participates in\nlength 2 cycle") 

    ggsave(fnames["direct_indirect"], plot, width = 7, height = 5, units = "in")

    plot = dfm %>%
        ggplot(data = ., aes(x = estimate, y = total_effect, label = parsed_label)) +
            cowplot::theme_cowplot(font_size = 12) +
            geom_smooth(method = "lm") +
            # geom_point(aes(x = estimate, y = total_effect, color = is_cyclic), alpha = .5) + 
            geom_point(aes(x = estimate, y = total_effect), alpha = .5) + 
            ggrepel::geom_text_repel(max.overlaps = 5, parse = TRUE) +
            # labs(x = "Direct effect", y = "Total effect", color = "Participates in\nlength 2 cycle") 
            labs(x = "Direct effect", y = "Total effect") 

    ggsave(fnames["direct_total"], plot, width = 7, height = 5, units = "in")


    model = lm(total_effect ~ estimate, data = dfm)
    print(broom::tidy(model)) 
    print(broom::glance(model)) 

    # return(fnames)
    return(dfm)
}

find_cycles_igraph <- function(graph, k) {
    ring <- igraph::graph.ring(k, TRUE)

    igraph::subgraph_isomorphisms(ring, graph)
}

find_cycles = function(causal_network, meta, k = 2, threshold = set_causal_threshold()) {

    graph = create_tidy_graph(
        causal_network %>%
            dplyr::filter(abs(estimate) > threshold),
        meta
    ) 

    cycles = find_cycles_igraph(graph, k)

    return(unique(names(unlist(cycles))))
}

compute_centrality = function(causal_network, meta, threshold = set_causal_threshold()) {

    graph = create_tidy_graph(
        causal_network %>% 
            dplyr::filter(abs(estimate) > threshold),
        meta
    )

    centrality = graph %>%
        tidygraph::activate(nodes) %>%
        tidygraph::mutate(
            page_rank = tidygraph::centrality_pagerank(),
            degree = tidygraph::centrality_degree(mode = "all"),
            out_degree = tidygraph::centrality_degree(mode = "out"),
            in_degree = tidygraph::centrality_degree(mode = "in"),
            betweenness = tidygraph::centrality_betweenness()
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::arrange(desc(page_rank))

    print(centrality)

    centrality_long =  centrality %>%
        tidyr::pivot_longer(names_to = "centrality_measure", values_to = "value", cols = c(page_rank:betweenness)) 

    print(centrality_long)

    plot = centrality_long %>%
        ggplot(data = ., aes(y = gene_group, x = value, fill = factor(stat(quantile)))) +
                ggridges::stat_density_ridges(
                        geom = "density_ridges_gradient", calc_ecdf = TRUE,
                        quantiles = 2, quantile_lines = TRUE
                ) +
                scale_fill_viridis_d(name = "quantile") +
                facet_wrap(~centrality_measure, nrow = 1, scales = "free_x") +
                # scale_x_continuous(trans = "log2") + 
                cowplot::theme_cowplot(font_size = 12) +
                cowplot::panel_border() +
                theme(axis.title.y = element_blank(), axis.title.x = element_blank())
            
    ggsave(file.path(figure_dir(), "centrality_measures_stratify.pdf"), plot, 
        width = 7, height = 5, units = "in")

    readr::write_tsv(centrality, file.path(txt_dir(), glue::glue("node_centrality_{threshold}.tsv")))

    return(centrality)
}


constraint_upstream_scatter = function(total_graph) {

    nodes = total_graph %>%
        tidygraph::activate(nodes) %>%
        mutate(
            out_degree = tidygraph::centrality_degree(mode = "out"),
            in_degree = tidygraph::centrality_degree(mode = "in")
        ) %>%
        as_tibble %>%
        mutate(
            upstream = out_degree - in_degree,
            upstream_bin = cut(upstream, breaks = c(-100, 0, 100, 400, 1300)),
            pLI_bin = cut(pLI, breaks = c(0, 0.1, .2, .3, .4, .5, .6, .8, .9, 1.0))
        )
    # pos = position_jitter(width = 0.3, seed = 1)

    print(nodes %>% dplyr::filter(is.na(upstream)))


    graph = ggplot(data = nodes, aes(x = pLI, y = upstream + 84, color = gene_group, label = name)) +
        # geom_point(position = pos, alpha = .5) +
        geom_point(alpha = .5) +
        cowplot::theme_cowplot(font_size = 12) +
        labs(
            x = "Constraint (pLI)",
            y = "log10(Outgoing connections - incoming connections + 84)",
            color = "Gene group"
        ) +
        add_gene_group_color_extra_ggplot2() +
        geom_hline(yintercept = 84, linetype = "dashed", alpha = .5) +
        scale_y_continuous(trans = "log10") +
        ggrepel::geom_label_repel(
            # position = pos,
            max.time = 20,
            seed = 100,
            max.iter = 8e5,
            size = 2.1,
            max.overlaps = 3
        )

    ggsave(file.path(figure_dir(), "upstream_constraint_scatter.pdf"), graph, 
        width = 7, height = 5, units = "in")

    graph = ggplot(data = nodes %>% dplyr::filter(!is.na(pLI) & !is.na(upstream_bin)),
        # aes(y = upstream_bin, x = pLI, fill = factor(stat(quantile)))) +
        aes(y = upstream_bin, x = pLI)) +
        # geom_point(position = pos, alpha = .5) +
        xlim(0, 1) +
        # ggridges::stat_density_ridges(
        #         geom = "density_ridges_gradient", calc_ecdf = TRUE,
        #         quantiles = 3, quantile_lines = TRUE,
        #         jittered_points = TRUE,
        #         position = ggridges::position_points_jitter(width = 0.001, height = 0),
        #         point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7
        # ) +
        # scale_fill_viridis_d(name = "Tertiles") +
        ggridges::geom_density_ridges(alpha = .5) +
        cowplot::theme_cowplot(font_size = 12) +
        labs(
            x = "Constraint (pLI)",
            y = "Outgoing connections - incoming connections"
            # color = "Quantile"
        ) 

    ggsave(file.path(figure_dir(), "ridgline_upstream_constraint.pdf"), graph, 
        width = 7, height = 5, units = "in")
        
}

plot_total_graph = function(total_graph, color, color_label, tag) {

    graph = total_graph %>%
        tidygraph::activate(nodes) %>%
        tidygraph::mutate(
            component_id = tidygraph::group_components()
        ) %>%
        tidygraph::filter(component_id == 1)


    plot = ggraph(graph, "stress") + 
                # geom_edge_arc2(
                geom_edge_link(
                    aes(colour = signed_weight),
                    # colour = gray,
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        alpha = .1,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) +
                scale_edge_colour_gradient2(low = "blue", mid = "gray", high = "red") + 
                geom_node_point(aes_string(color = color), alpha = .4) +
                labs(color = color_label) +
                theme_graph(base_family = "Helvetica")

    fname = file.path(figure_dir(), glue::glue("ggplot_network_total_{tag}.pdf"))
    ggsave(
        fname,
        plot,
        width = 14,
        height = 10,
        units = "in"
    )
}

create_pathfinder_graph = function(causal_network, pathfindr_results, cluster_membership, meta, threshold = set_causal_threshold()) {

    select_pathways = c(
        "T cell receptor signaling pathway",
        "TNF signaling pathway",
        "NF-kappa B signaling pathway",
        "Ribosome",
        "Th1 and Th2 cell differentiation"
    )
    edges = causal_network %>%
        dplyr::filter(abs(estimate) > threshold) %>%
        dplyr::inner_join(
            cluster_membership %>% dplyr::select(from_cluster = cluster, gene_name),
            by = c("row" = "gene_name")
        ) %>%
        dplyr::inner_join(
            cluster_membership %>% dplyr::select(to_cluster = cluster, gene_name),
            by = c("col" = "gene_name")
        ) 

    pathfindr_results = pathfindr_results %>%
        dplyr::filter(Term_Description %in% .env[["select_pathways"]]) %>%
        dplyr::inner_join(
            cluster_membership %>% dplyr::select(from_cluster = cluster, gene_name),
            by = c("KO" = "gene_name")
        )

    print(pathfindr_results)

    pathfinder_edges = pathfindr_results %>%
        # dplyr::rename(from = KO) %>%
        dplyr::rename(from = from_cluster) %>%
        dplyr::filter(Fold_Enrichment >= 4) %>%
        dplyr::select(from, to = Term_Description) %>%
        dplyr::distinct(.)

    print(pathfinder_edges)

    edges = edges %>%
        # dplyr::select(from = row, to = col) %>%
        dplyr::select(from = from_cluster, to = to_cluster) %>%
        dplyr::distinct(.) %>%
        dplyr::bind_rows(pathfinder_edges)

    active_nodes = unique(c(edges$from, edges$to))

    pathfinder_nodes = pathfinder_edges %>%
        # dplyr::filter(from %in% .env[["active_nodes"]]) %>%
        dplyr::select(name = to) %>%
        dplyr::mutate(cluster = "downstream") %>%
        dplyr::distinct(.)

    # status = meta %>%
    #     dplyr::filter(!is_control) %>%
    #     dplyr::distinct(name = KO, gene_group) %>%
    #     dplyr::inner_join(cluster_membership, by = c("name" = "gene_name")) %>%
    status = cluster_membership %>% 
        dplyr::select(name = cluster) %>%
        dplyr::bind_rows(pathfinder_nodes) %>%
        dplyr::distinct(.)

    print(pathfinder_nodes)
    print(status)

    logger::log_info('here')

    # browser()
    graph = tidygraph::tbl_graph(
        status,
        edges
    )
   

    plot = ggraph(graph, "auto") + 
                geom_edge_arc2(
                # geom_edge_link(
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        alpha = .3,
                        strength = 0.3,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) +
                scale_edge_colour_gradient2(low = "blue", mid = "gray", high = "red") + 
                geom_node_point(aes(color = name)) +
                geom_node_text(
                    aes(label = name, color = name),
                    size = 4.5,
                    # colour = "black", 
                    # family = "serif",
                    check_overlap = TRUE,
                    repel = TRUE
                ) + 
                theme_graph(base_family = "Helvetica")

    fname = file.path(figure_dir(), glue::glue("pathfinder_cluster_network.pdf"))
    ggsave(
        fname,
        plot,
        width = 11,
        height = 6,
        units = "in"
    )
}

create_module_gwas_graph = function(causal_network, cluster_membership, downstream_immune_GWAS, meta, threshold = set_causal_threshold()) {

    select_modules = c(
                        "4A"
    )

    cluster_genes = cluster_membership %>%
        dplyr::filter(cluster == .env[["select_modules"]]) %>%
        dplyr::pull(gene_name)

    edges = causal_network %>%
        dplyr::filter(abs(estimate) > threshold) %>%
        dplyr::filter(col %in% cluster_genes & row %in% cluster_genes) %>%
        dplyr::rename(from = row, to = col)

    downstream_edges = downstream_immune_GWAS[[2]] %>%
        dplyr::filter(cluster == .env[["select_modules"]]) %>%
        dplyr::select(
            from,
            to,
            estimate = beta,
            `Multiple sclerosis`,
            `Type 1 diabetes`,
            `Lupus`,
            `Rheumatoid arthritis`
        )

    print(downstream_edges %>% dplyr::distinct(to))

    downstream_edges %>%
        dplyr::filter(
            `Multiple sclerosis` |
            `Type 1 diabetes` | 
            `Lupus` |
            `Rheumatoid arthritis`
        ) %>%
        dplyr::filter(!stringr::str_detect(to, "HLA-")) %>%
        dplyr::distinct(to) %>%
        print

    downstream_edges = downstream_edges %>%
        dplyr::filter(
            `Multiple sclerosis` 
            # `Type 1 diabetes` 
        ) %>%
        dplyr::filter(!stringr::str_detect(to, "HLA-")) 

    edges = edges %>%
        dplyr::bind_rows(downstream_edges) 

    status = tibble::tibble(
                name = unique(c(edges$from, edges$to))
            )

    gwas_genes = downstream_immune_GWAS[[3]] %>%
        dplyr::filter(
            `Multiple sclerosis` 
            # `Type 1 diabetes` 
        ) %>%
        dplyr::pull(gene_name)

    status = status %>%
        dplyr::mutate(is_gwas = name %in% gwas_genes)

    graph = tidygraph::tbl_graph(
        status,
        edges
    )

    manual_layout = tibble::tribble(
        ~name, ~x, ~y,
        "IRF4", -5, 8,
        "JAK3", 5, 8,
        "KMT2A", 0, 16,
        "STAT5A", 2, 8,
        "STAT5B", -2, 8,
        "IL2RA", 10, 8,
        "NCR3", 14, 0,
        "IER3", -10, 0,
        "PSMB9", -14, 0,
        "UQCC2", 16, 0,
        "BAK1", -8, 0,
        "MDC1", -6, 0,
        "TCF19", -4, 0,
        "BTN3A3", -12, 0,
        "ABCF1", -1, 0,
        "MRPS18B", 2, 0,
        "GABBR1", 5, 0,
        "GNL1", 7, 0,
        "AGPAT1", 9, 0,
        "AIF1", 12, 0
    )

    graph = graph %>%
        tidygraph::activate(nodes) %>%
        tidygraph::inner_join(manual_layout, by = "name")

    # plot = ggraph(graph, "tree") + 
    plot = ggraph(graph, x = x, y = y) + 
                geom_edge_arc2(
                # geom_edge_link(
                        aes(colour = estimate),
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        alpha = .6,
                        strength = 0.1,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) +
                scale_edge_colour_gradient2(low = "blue", mid = "gray", high = "red") + 
                geom_node_point() +
                geom_node_text(
                    # aes(label = name, color = is_gwas),
                    aes(label = name),
                    size = 4.0,
                    # colour = "black", 
                    # family = "serif",
                    check_overlap = TRUE,
                    repel = TRUE
                ) + 
                theme_graph(base_family = "Helvetica") +
                theme(legend.position = "top")

    fname = file.path(figure_dir(), glue::glue("cluster_4_network_manual_layout.pdf"))
    ggsave(
        fname,
        plot,
        width = 9,
        height = 6,
        units = "in"
    )

    return(graph)
}

create_module_pathway_graph = function(causal_network, mashr, cluster_membership, downstream_pathway = "04630", meta, threshold = set_causal_threshold(), lfsr_threshold = set_lfsr_threshold()) {

    select_modules = c(
                        "4A"
    )

    pathway_entrez_ids = as.list(org.Hs.eg.db::org.Hs.egPATH2EG)[[downstream_pathway]] %>%
        as.character

    print("here 1")
    pathway_gene_names = AnnotationDbi::mapIds(
                            org.Hs.eg.db::org.Hs.eg.db,
                            keys = pathway_entrez_ids,
                            column = "SYMBOL",
                            keytype = "ENTREZID"
                        )

    cluster_genes = cluster_membership %>%
        dplyr::filter(cluster == .env[["select_modules"]]) %>%
        dplyr::pull(gene_name)

    print("here 1.5")

    edges = causal_network %>%
        dplyr::filter(abs(estimate) > threshold) %>%
        dplyr::filter(col %in% cluster_genes & row %in% cluster_genes) %>%
        dplyr::rename(from = row, to = col, signed_weight = estimate)

    pm = ashr::get_pm(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "signed_weight", -to) %>%
        dplyr::mutate(weight = abs(signed_weight))

    lfsr = ashr::get_lfsr(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "lfsr", -to) %>%
        dplyr::inner_join(pm) %>%
        dplyr::inner_join(
            meta %>% dplyr::distinct(KO, gene_group), 
            by = c("from" = "KO")
        )

    print("here 2")
    downstream_edges = lfsr %>%
        dplyr::filter(lfsr < lfsr_threshold) %>%
        dplyr::filter(to %in% .env[["pathway_gene_names"]]) %>%
        dplyr::filter(from %in% cluster_genes) %>%
        dplyr::select(from, to, signed_weight, weight)

    print(downstream_edges %>% dplyr::distinct(to))


    edges = edges %>%
        dplyr::bind_rows(downstream_edges)
    
    print(edges)

    status = tibble::tibble(
                name = unique(c(edges$from, edges$to))
            )

    graph = tidygraph::tbl_graph(
        status,
        edges
    )

    level_0 = 16
    level_1 = 8
    level_2 = 0
    level_3 = -4
    level_4 = -2
    manual_layout = tibble::tribble(
        ~name, ~x, ~y,
        "IRF4", -5, level_1,
        "JAK3", 5, level_1,
        "KMT2A", 0, level_0,
        "STAT5A", 2, level_1,
        "STAT5B", -2, level_1,
        "IL2RA", 10, level_1,
        "JAK1", 5, level_2,
        "SOS1", 16, level_2,
        "CBLB", -20, level_2,
        "IL7R", -18, level_3,
        "MYC", -16, level_2,
        "JAK2", -14, level_2,
        "IL10RA", 15, level_3,
        "PTPN6", -12, level_2,
        "IL26", -10, level_3,
        "PIAS1", -9, level_2,
        "OSM", 2, level_2,
        "IL12RB2", 17, level_3,
        "IFNG", -2, level_3,
        "IL6R", 0, level_3,
        "IL21", -1, level_3,
        "PIK3CB", 15, level_2,
        "IFNGR2", 1, level_3,
        "CCND3", -7, level_2,
        "SOCS3", -4, level_2,
        "CCND2", 4, level_2,
        "CISH", -5, level_2,
        "IL2RB", 3, level_3,
        "IL22", 9, level_3,
        "IL23R", 7, level_3,
        "IL6ST", 18, level_3,
        "LIF", 20, level_2,
        "SOCS2", 7, level_2,
        "STAT4", 4.5, level_2,
        "IL4R", 22, level_3,
        "SOCS1", 2.5, level_2,
        "IL2RG", 23, level_3,
        "IL15", 3.5, level_3,
        "IL13", 21, level_3,
        "IFNGR1", 8, level_3,
        "SPRED2", 9, level_2,
        "IL15RA", 24, level_3,
        "IFNAR2", 25, level_3,
        "IL9", 19, level_3,
        "IL21R", -11, level_3,
        "AKT3", 26, level_2,
        "IL12RB1", 14, level_3,
        "IL5", 6, level_3
    )

    manual_layout %>%
        dplyr::add_count(x, y) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::filter(n > 1) %>%
        as.data.frame %>%
        print

    graph = graph %>%
        tidygraph::activate(nodes) %>%
        tidygraph::inner_join(manual_layout, by = "name")

    # plot = ggraph(graph, "tree") + 
    plot = ggraph(graph, x = x, y = y) + 
                geom_edge_arc2(
                # geom_edge_link(
                        aes(colour = signed_weight),
                        arrow = arrow(
                                angle = 15,
                                length = unit(0.13, "inches"),
                                # ends = "last",
                                type = "closed"
                        ),
                        alpha = .4,
                        strength = 0.1,
                        start_cap = circle(2.6, 'mm'),
                        end_cap = circle(2.6, 'mm')
                ) +
                scale_edge_colour_gradient2(low = "blue", mid = "gray", high = "red") + 
                geom_node_point() +
                geom_node_text(
                    # aes(label = name, color = is_gwas),
                    aes(label = name),
                    size = 3.0,
                    # colour = "black", 
                    # family = "serif",
                    check_overlap = TRUE,
                    repel = TRUE
                ) + 
                theme_graph(base_family = "Helvetica") +
                theme(legend.position = "top", legend.key.width = unit(1, "cm"))

    fname = file.path(figure_dir(), glue::glue("cluster_4_network_pathway_manual_layout.pdf"))
    ggsave(
        fname,
        plot,
        width = 12,
        height = 8,
        units = "in"
    )

    return(graph)
}

compute_modules_with_downstream = function(mashr, cluster_membership, lfsr_threshold = set_lfsr_threshold()) {

    cluster_labels = unique(cluster_membership$cluster)

    cluster_genes = cluster_membership %>%
        dplyr::group_by(cluster) %>%
        dplyr::group_map(~.x$gene_name) %>%
        purrr::set_names(cluster_labels)

    pm = ashr::get_pm(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "signed_weight", -to) %>%
        dplyr::mutate(weight = abs(signed_weight))

    lfsr = ashr::get_lfsr(mashr) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(to = mashr$readout_gene) %>%
        tidyr::pivot_longer(names_to = "from", values_to = "lfsr", -to) %>%
        dplyr::inner_join(pm)

    downstream_edges = lfsr %>%
        dplyr::filter(lfsr < lfsr_threshold) %>%
        dplyr::select(from, to, signed_weight, weight) %>%
        dplyr::inner_join(
            cluster_membership %>% dplyr::select(from_cluster = cluster, gene_name),
            by = c("from" = "gene_name")
        )

    cluster_genes = purrr::map(cluster_labels, ~{
            upstream = cluster_genes[[.x]]
            downstream = downstream_edges %>%
                dplyr::filter(from_cluster == .env[[".x"]]) %>%
                dplyr::pull(to) %>%
                unique

            return(c(upstream, downstream))
        }) %>%
        purrr::set_names(cluster_labels)

    return(cluster_genes)
}
