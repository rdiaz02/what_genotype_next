load("transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15.RData")
transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["CBN"]]
## All 0
transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["CBN"]][6, ]

transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["CBN_td"]][6, ]

transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["MHN"]][8, ]

transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["MHN_td"]][8, ]


transition_matrices_num1_dataname_GBM_coo_cores_20_cores_cbn_version_max_genes_15[["OT"]]




load("transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15.RData")

## all to a 1
transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["OT"]][14, ]


G_dna_Damage_control, G_JNK, G_Small_GTPase_signaling, G_Wnt_Notch_signaling


transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["CBN"]]


load("transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15.RData")


which(transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15[["CBN"]][10, ] != 0)

transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15[["CBN"]]["G_5q_n_, G_Xp_n_", ]

transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15[["OT"]]["G_5q_n_, G_Xp_n_", ]


uu <- transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15[["CBN"]]["G_5q_n_, G_Xp_n_", ]


vv <- transition_matrices_num6_dataname_Ov_CNV_cores_20_cores_cbn_version_max_genes_15[["MHN"]]["G_5q_n_, G_Xp_n_", ]
vv[vv!= 0]
uu[uu!= 0]

## The names are given in different order

uu2 <- uu[uu!= 0]
vv2 <- vv[vv!= 0]

uu2 <- c(uu2, 0, 0)
vv2 <- vv2[c(5, 1, 3,  2, 4)]

jensen_shannon(uu2, vv2)


transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["OT"]][14,]

transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["OT"]]["G_dna_Damage_control, G_JNK, G_Small_GTPase_signaling, G_Wnt_Notch_signaling",]


transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["CBN"]]["G_dna_Damage_control, G_JNK, G_Small_GTPase_signaling, G_Wnt_Notch_signaling",]
















yy <- transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["OT"]]["G_dna_Damage_control, G_JNK, G_Small_GTPase_signaling, G_Wnt_Notch_signaling",]
xx <- transition_matrices_num5_dataname_Pan_pa_cores_20_cores_cbn_version_max_genes_15[["MHN"]]["G_dna_Damage_control, G_JNK, G_Small_GTPase_signaling, G_Wnt_Notch_signaling",]

yy[yy != 0]
xx[xx != 0]

jensen_shannon(c(1, 0, 0), xx[xx != 0][c(3, 1, 2)])
