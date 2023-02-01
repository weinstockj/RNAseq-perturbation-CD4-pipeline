# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c(
    "tibble",
    "dplyr",
    "purrr", 
    "DESeq2",
    "stringr",
    "rtracklayer",
    "vroom",
    "arrow",
    "ggplot2",
    "patchwork",
    "mashr",
    "logger",
    "enrichR",
    "ggridges",
    "ggraph"
  ), # packages that your targets need to run
  format = "qs" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "slurm")
# options(clustermq.template = "clustermq.tmpl")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Load the R scripts with your custom functions:
source("R/functions.R")
source("R/plot_network.R")
source("R/enrichment_analyses.R")
source("R/mashr_analysis.R")
source("R/proliferation.R")
source("R/joint_downstream_model.R")
source("R/validate_network.R")
source("R/ABC_enrichment.R")
source("R/pathway_plot.R")

# Replace the target list below with your own:
list(
    tar_target(
        name = rna,
        command =  read_rna(rna_seq_location()),
        format = "parquet" 
    ),
    tar_target(
        name = txdb,
        command = get_txdb(gencode_location()),
        format = "parquet"
    ),
    tar_target(
        name = meta,
        command = parse_metadata(rna, vroom::vroom(master_metadata())),
        format = "parquet"
    ),
    tar_target(
        name = counts,
        command = get_counts(rna)
    ),
    tar_target(
        name = dds,
        command = make_dds(counts, meta)
    ),
    tar_target(
        name = filtered_dds,
        command = filter_low_counts(dds)
    ),
    tar_target(
        name = pca,
        command = apply_pca(SummarizedExperiment::assay(DESeq2::vst(filtered_dds, blind = TRUE)))
    ),
    tar_target(
        name = dds_with_pcs,
        command = modify_dds(filtered_dds, pca)
    ),
    tar_target(
        name = diffeq_no_pcs,
        command = run_diffeq(filtered_dds)
    ),
    tar_target(
        name = diffeq_with_pcs,
        command = run_diffeq(dds_with_pcs)
    ),
    tar_target(
        name = gene_lookup,
        command = get_readout_gene_lookup(filtered_dds, txdb),
        format = "parquet"
    ),
    tar_target(
        name = plot_pca_,
        command = plot_pca(pca, meta)
    ),
    tar_target(
        name = results_with_pcs,
        command = parse_results(diffeq_with_pcs)
    ),
    tar_target(
        name = results_no_pcs,
        command = parse_results(diffeq_no_pcs)
    ),
    tar_target(
        name = compare_dge_,
        command = compare_dge(results_with_pcs, results_no_pcs, txdb),
        format = "parquet"
    ),
    tar_target(
        name = plot_compare_dge_,
        command = plot_compare_dge(compare_dge_),
        format = "file"
    ),
    tar_target(
        name = mashr,
        command = apply_mashr(results_with_pcs, txdb)
    ),
    tar_target(
        name = mashr_no_pcs,
        command = apply_mashr(results_no_pcs, txdb)
    ),
    tar_target(
        name = joint_downstream_model,
        command = read_joint_downstream_model(txdb)
    ),
    tar_target(
        name = pairwise,
        command = extract_pairwise_effects(mashr)
    ),
    tar_target(
        name = pairwise_joint_downstream,
        command = extract_pairwise_effects(joint_downstream_model)
    ),
    tar_target(
        name = plot_pairwise,
        command = plot_pairwise_effects(pairwise, meta)
    ),
    tar_target(
        name = plot_pairwise_joint_downstream,
        command = plot_pairwise_effects(pairwise_joint_downstream, meta, tag = "joint_downstream")
    ),
    tar_target(
        name = plot_covariance_pis_ ,
        command = plot_covariance_pis(mashr)
    ),
    tar_target(
        name = causal_network,
        command = read_direct_effects(direct_effects_location())
    ),
    tar_target(
        name = write_out_total_covar,
        command = compute_covariance_reorder(causal_network, txdb)
    ),
    tar_target(
        name = length_2_cycle_genes,
        command = find_cycles(causal_network, meta, 2)
    ),
    tar_target(
        name = length_3_cycle_genes,
        command = find_cycles(causal_network, meta, 3)
    ),
    tar_target(
        name = compute_centrality_,
        command = compute_centrality(causal_network, meta)
    ),
    tar_target(
        name = plot_causal_network,
        command = plot_network(causal_network, meta, threshold = 0.025),
        format = "file"
    ),
    tar_target(
        name = plot_causal_network_dge,
        command = plot_network(causal_network, meta, threshold = NULL, results = results_with_pcs, txdb),
        format = "file"
    ),
    tar_target(
        name = network_degree,
        command = compute_degree(causal_network)
    ),
    tar_target(
        name = stratify_network_degree_,
        command = stratify_degree_by_gene_group(network_degree, meta)
    ),
    tar_target(
        name = total_effects,
        command = read_total_effects()
    ),
    tar_target(
        name = compare_total_effects_,
        command = compare_total_effects(total_effects, results_with_pcs, txdb)
    ),
    tar_target(
        name = plot_direct_indirect_,
        command = plot_direct_indirect(causal_network, total_effects, results_with_pcs, txdb, length_2_cycle_genes)
    ),
    tar_target(
        name = plot_indegree_outdegree_scatter_,
        command = plot_indegree_outdegree_scatter(compute_centrality_)
    ),
    tar_target(
        name = plot_causal_scatter,
        command = compare_pairwise_to_direct_effects(pairwise, causal_network)
    ),
    tar_target(
        name = plot_causal_scatter_joint_downstream,
        command = compare_pairwise_to_direct_effects(pairwise_joint_downstream, causal_network, tag = "joint_downstream")
    ),
    tar_target(
        name = plot_causal_scatter_distance,
        command = compare_pairwise_to_distance(pairwise, causal_network, meta)
    ),
    tar_target(
        name = plot_causal_scatter_distance_joint_downstream,
        command = compare_pairwise_to_distance(pairwise_joint_downstream, causal_network, meta, tag = "joint_downstream")
    ),
    tar_target(
        name = plot_degree_scatter,
        command = compare_degree_estimates(results_with_pcs, causal_network)
    ),
    tar_target(
        name = ko_specific_effects,
        command = extract_ko_specific_effects(mashr, txdb)
    ),
    tar_target(
        name = ko_specific_effects_no_pcs,
        command = extract_ko_specific_effects(mashr_no_pcs, txdb)
    ),
    tar_target(
        name = ko_specific_effects_joint_downstream,
        command = extract_ko_specific_effects(joint_downstream_model, txdb)
    ),
    tar_target(
        name = downstream_indegree,
        command = compute_downstream_indegree(mashr, txdb)
    ),
    tar_target(
        name = downstream_indegree_by_group,
        command = compute_downstream_indegree_by_group(mashr, txdb, meta, type = "lfsr", pm_threshold = 5e-3)
    ),
    tar_target(
        name = downstream_indegree_by_group_joint_downstream,
        command = compute_downstream_indegree_by_group(joint_downstream_model, txdb, meta, type = "lfsr", lfsr_threshold = 5e-3)
    ),
    tar_target(
        name = constraint,
        command = vroom::vroom(constraint_metrics_location(), col_select = c(gene, pLI)) %>%
            dplyr::rename(gene_name = gene)
    ),
    tar_target(
        name = proliferation,
        command = read_proliferation()
    ),
    tar_target(
        name = plot_proliferation_scatter_,
        command = plot_proliferation(meta, proliferation)
    ),
    tar_target(
        name = enrichment_,
        command = enrichment(ko_specific_effects)
    ),
    tar_target(
        name = enrichment_joint_downstream_,
        command = enrichment(ko_specific_effects_joint_downstream)
    ),
    tar_target(
        name = plot_ko_specific_effects_,
        command = plot_ko_specific_effects(ko_specific_effects, meta, tag = "regressed_pcs")
    ),
    tar_target(
        name = plot_ko_specific_effects_no_pcs_,
        command = plot_ko_specific_effects(ko_specific_effects_no_pcs, meta, tag = "no_pcs")
    ),
    tar_target(
        name = plot_ko_specific_effects_joint_downstream_,
        command = plot_ko_specific_effects(ko_specific_effects_joint_downstream, meta, tag = "joint_downstream")
    ),
    tar_target(
        name = plot_ko_specific_effects_proliferation_,
        command = plot_proliferation_ko_specific(ko_specific_effects, proliferation, meta)
    ),
    tar_target(
        name = plot_ko_specific_effects_proliferation_joint_downstream,
        command = plot_proliferation_ko_specific(ko_specific_effects_joint_downstream, proliferation, meta, tag = "joint_downstream")
    ),
    tar_target(
        name = plot_enrichments_,
        command = plot_enrichments(enrichment_)
    ),
    tar_target(
        name = plot_enrichments_joint_downstream,
        command = plot_enrichments(enrichment_joint_downstream_, tag = "joint_downstream")
    ),
    tar_target(
        name = downstream_gene_enrichment,
        command = disease_enrichment_of_coregulation(downstream_indegree, diffeq_no_pcs, gene_lookup)
    ),
    tar_target(
        name = plot_strata_,
        command = plot_strata(downstream_gene_enrichment)
    ),
    tar_target(
        name = plot_downstream_cluster_,
        command = downstream_cluster(
            mashr, 
            txdb, 
            row_cluster_k = 8,
            col_cluster_k = 8,
            min_overlap_prop = 0,
            min_overlap = 0,
            enrichment_function = "enrichr", 
            tag = "diffeq"
        )
    ),
    tar_target(
        name = plot_downstream_cluster_joint_downstream_,
        command = downstream_cluster(
            joint_downstream_model, 
            txdb, 
            row_cluster_k = 7,
            col_cluster_k = 4,
            lower_limit_n_downstream = 100,
            min_overlap_prop = 0,
            min_overlap = 0,
            enrichment_function = "enrichr", 
            tag = "joint_downstream"
        )
    ),
    tar_target(
        name = plot_downstream_cluster_predefined_,
        command = downstream_cluster(
            mashr,
            txdb,
            enrichment_function = "predefined",
            gene_set = iei_gene_list(),
            gene_names = iei_genes(),
            row_cluster_k = 5,
            min_overlap = 1,
            min_overlap_prop = .03,
            gene_group_label = "predefined",
            tag = "diffeq"
        )
    ),
    tar_target(
        name = plot_downstream_cluster_predefined_joint_downstream_,
        command = downstream_cluster(
            joint_downstream_model,
            txdb,
            enrichment_function = "predefined",
            gene_set = iei_gene_list(),
            gene_names = iei_genes(),
            row_cluster_k = 5,
            min_overlap = 1,
            min_overlap_prop = .03,
            gene_group_label = "predefined",
            tag = "joint_downstream"
        )
    ),
    tar_target(
        name = plot_downstream_cluster_iei_,
        command = downstream_cluster(
                    mashr,
                    txdb,
                    gene_names = iei_genes(),
                    gene_group_label = "iei",
                    lower_limit_n_downstream = 10,
                    row_cluster_k = 8,
                    col_cluster_k = 8,
                    # min_overlap = 5,
                    # min_overlap_prop = .18,
                    enrichment_function = "enrichr",
                    tag = "diffeq",
                    database = "GO_Biological_Process_2021"
                )
    ),
    tar_target(
        name = plot_downstream_cluster_iei_joint_downstream_,
        command = downstream_cluster(
                    joint_downstream_model,
                    txdb,
                    gene_names = iei_genes(),
                    gene_group_label = "iei",
                    lower_limit_n_downstream = 20,
                    row_cluster_k = 6,
                    col_cluster_k = 4,
                    # min_overlap = 5,
                    # min_overlap_prop = .18,
                    enrichment_function = "enrichr",
                    tag = "joint_downstream",
                    database = "GO_Biological_Process_2021"
                )
    ),
    tar_target(
        name = plot_downstream_cluster_pics_,
        command = downstream_cluster(
                    mashr,
                    txdb,
                    gene_names = pics_ai_genes(),
                    gene_group_label = "pics",
                    row_cluster_k = 8,
                    col_cluster_k = 8,
                    # min_overlap = 5,
                    # min_overlap_prop = .18,
                    enrichment_function = "enrichr",
                    tag = "diffeq",
                    "GO_Biological_Process_2021"
                )
    ),
    tar_target(
        name = plot_downstream_cluster_pics_joint_downstream_,
        command = downstream_cluster(
                    joint_downstream_model,
                    txdb,
                    gene_names = pics_ai_genes(),
                    gene_group_label = "pics",
                    row_cluster_k = 7,
                    col_cluster_k = 7,
                    # min_overlap = 5,
                    # min_overlap_prop = .18,
                    enrichment_function = "enrichr",
                    tag = "joint_downstream",
                    "GO_Biological_Process_2021"
                )
    ),
    tar_target(
        name = gene_mask,
        command = create_gene_mask(gene_lookup, meta),
        format = "feather"
    ),
    tar_target(
        name = write_gene_mask_,
        command = write_gene_mask(gene_mask),
        format = "file"
    ),
    tar_target(
        name = read_yazar_,
        command = read_yazar_trans_eqtls()
    ),
    tar_target(
        name = read_dice_,
        command = read_dice(),
        format = "feather"
    ),
    tar_target(
        name = blood_specific_genes,
        command = identify_blood_specific_genes()
    ),
    tar_target(
        name = tcell_specific_data,
        command = identify_tcell_specific_genes(blood_specific_genes, read_dice_),
        format = "feather"
    ),
    tar_target(
        name = tcell_specific_genes,
        command = tcell_specific_data %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = THSTAR_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(THSTAR_non_th_ratio > 1.3) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = CD4_STIM_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(CD4_STIM_THSTAR_ratio > 1.3) %>%
            # dplyr::filter(CD4_STIM_non_th_ratio > 1.2) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = CD4_NAIVE_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(CD4_NAIVE_THSTAR_ratio > 1.3) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = TH1_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(TH1_THSTAR_ratio > 1.3) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = TH2_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(TH2_THSTAR_ratio > 1.3) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    tar_target(
        name = TH17_enriched_genes,
        command = tcell_specific_data %>%
            dplyr::filter(TH17_THSTAR_ratio > 1.3) %>%
            dplyr::pull(gene_id) %>%
            convert_ensembl_to_symbol(remove_ensembl_version_txdb(txdb))
    ),
    # tar_target(
    #     name = annotate_t_cell_expression,
    #     command = dice_max_cell_type(read_dice_, blood_specific_genes),
    #     format = "feather"
    # ),
    tar_target(
        name = read_eqtlgen_,
        command = read_eqtlgen_trans_eqtls()
    ),
    tar_target(
        name = expression,
        command = extract_expression_by_gene(results_with_pcs, txdb),
    ),
    tar_target(
        name = forest_plot_downstream_indegree_,
        command = forest_plot_downstream_indegree(
            downstream_indegree_by_group,
            expression,
            meta,
            constraint,
            iei_genes = iei_genes(),
            gwas_genes = pics_ai_genes(),
            trans_egenes = unique(read_eqtlgen_)$gene_name,
            blood_specific_genes = convert_ensembl_to_symbol(blood_specific_genes, remove_ensembl_version_txdb(txdb)),
            tcell_specific_genes = tcell_specific_genes,
            tag = "diffeq_lfsr"
        )
    ),
    tar_target(
        name = forest_plot_downstream_indegree_joint_downstream_,
        command = forest_plot_downstream_indegree(
            downstream_indegree_by_group_joint_downstream,
            expression,
            meta,
            constraint,
            iei_genes = iei_genes(),
            gwas_genes = pics_ai_genes(),
            trans_egenes = unique(read_eqtlgen_)$gene_name,
            blood_specific_genes = convert_ensembl_to_symbol(blood_specific_genes, remove_ensembl_version_txdb(txdb)),
            tcell_specific_genes = tcell_specific_genes,
            tag = "joint_downstream_lfsr"
        )
    ),
    tar_target(
        name = forest_plot_downstream_indegree_cell_type_specific_joint_downstream_,
        command = forest_plot_downstream_indegree_cell_type_specific(
            downstream_indegree_by_group_joint_downstream,
            expression,
            meta,
            TH1_genes = TH1_enriched_genes,
            TH2_genes = TH2_enriched_genes,
            TH17_genes = TH17_enriched_genes,
            CD4_stimulated_genes = CD4_STIM_enriched_genes,
            CD4_naive_genes = CD4_NAIVE_enriched_genes,
            tag = "joint_downstream_lfsr"
        )
    ),
    tar_target(
        name = forest_plot_causal_centrality_,
        command = forest_plot_causal_centrality(
            compute_centrality_,
            expression,
            meta,
            constraint,
            iei_genes = iei_genes(),
            gwas_genes = pics_ai_genes(),
            trans_egenes = unique(read_eqtlgen_)$gene_name
        )
    ),
    # total graph analyses
    tar_target(
        name = total_graph,
        command = create_total_tidy_graph(
            causal_network,
            joint_downstream_model,
            expression,
            meta,
            constraint,
            iei_genes = iei_genes(),
            gwas_genes = pics_ai_genes(),
            trans_egenes = unique(read_eqtlgen_)$gene_name
        )
    ),
    tar_target(
        name = scatter_upstream_constraint_,
        command = constraint_upstream_scatter(total_graph)
    ),
    tar_target(
        name = plot_total_graph_gwas,
        command = plot_total_graph(total_graph, color = "is_gwas", color_label = "AI GWAS gene", tag = "gwas")
    ),
    tar_target(
        name = plot_total_graph_iei,
        command = plot_total_graph(total_graph, color = "is_iei", color_label = "IEI gene", tag = "iei")
    ),
    tar_target(
        name = plot_total_graph_eqlten,
        command = plot_total_graph(total_graph, color = "is_trans_egene", color_label = "eQTLgen trans-eGene", tag = "eqtlgen")
    ),
    tar_target(
        name = hbase,
        command = read_hbase(meta)
    ),
    tar_target(
        name = stringdb,
        command = read_stringdb(meta)
    ),
    tar_target(
        name = determine_enrichment_hbase,
        command = iterate_hbase_thresholding(meta, hbase, causal_network)
    ),
    tar_target(
        name = ABC_GRN,
        command = read_ABC_GRN(meta)
    ),
    tar_target(
        name = DAC,
        command = read_DAC()
    ),
    tar_target(
        name = DAC_ABC_GRN,
        command = merge_DAC_ABC_GRN(ABC_GRN, DAC)
    ),
    tar_target(
        name = determine_enrichment_DAC_ABC_GRN,
        command = iterate_ABC_GRN_thresholding(meta, DAC_ABC_GRN, causal_network)
    ),
    tar_target(
        name = matrix_factorization,
        command = read_matrix_factorization(txdb)
    ),
    tar_target(
        name = matrix_factorization_enrichments,
        command = enrich_loadings(matrix_factorization, txdb)
    ),
    tar_target(
        name = program_labels,
        command = process_loading_enrichments(matrix_factorization_enrichments)
    ),
    tar_target(
        name= plot_loadings_,
        command = plot_loadings(matrix_factorization, program_labels, txdb, meta)
    ),
    tar_target(
        name = pathway_plot_,
        command = plot_pathway(joint_downstream_model, txdb, "DR1", "R-HSA-202424")
    ),
    tar_target(
        name = write_out_experiment_,
        command = write_out_experiment(meta, diffeq_with_pcs, pca, mashr, results_with_pcs, results_no_pcs, txdb),
        format = "file"
    )
)
