library(RColorBrewer)
algorithm_colors <- c(
    "MOBSTER" = "darkorange2",
    "DECODE" = "magenta4"
)
cluster_shapes <- c(
    "Cluster_1" = 21,
    "Cluster_2" = 22,
    "Cluster_3" = 23,
    "Cluster_4" = 24,
    "Cluster_5" = 25
)
cluster_labels <- c(
    "Cluster_1" = "1",
    "Cluster_2" = "2",
    "Cluster_3" = "3",
    "Cluster_4" = "4",
    "Cluster_5" = "5"
)
cluster_colors <- c(
    "1" = "#D55E00",
    "2" = "#0072B2",
    "3" = "#009E73",
    "4" = "#E69F00",
    "5" = "#56B4E9"
)
ICGC_cohorts <- c(
    "CNS-PiloAstro", "Liver-HCC", "Kidney-RCC", "Prost-AdenoCA", "Lymph-BNHL",
    "Panc-Endocrine", "Panc-AdenoCA", "CNS-Medullo", "Breast-AdenoCA", "Stomach-AdenoCA",
    "Skin-Melanoma", "Lymph-CLL", "Eso-AdenoCA", "Myeloid-AML", "Head-SCC",
    "Billary-AdenoCA", "Ovary-AdenoCA", "Billary-AdenoCA", "Bone-Benign", "Bone-Epith",
    "Bone-Osteosarc", "Breast-DCIS", "Breast-LobularCA", "Myeloid-MDS", "Myeloid-MPN"
)
ICGC_cohort_colors <- setNames(
    colorRampPalette(brewer.pal(11, "Spectral"))(length(ICGC_cohorts)),
    ICGC_cohorts
)
POG_cohorts <- c(
    "COAD", "IDC", "BCC", "ATRT", "MCL", "NPC", "LUAD",
    "ANSC", "CHOL", "APAD", "PAAD", "SACA", "GRCT", "ILC",
    "LUSC", "LUACC", "EMBC", "MBC", "PRAD", "PTCLNOS", "HGSOC",
    "UEC", "AODG", "WDLS", "AECA", "PEMESO", "ACC", "SOC",
    "MFH", "ODG", "DDLS", "LGSOC", "DA", "MEL", "ES",
    "CACC", "CEAIS", "OS", "FL", "UM", "UCEC", "BLAD",
    "LMS", "ODGC", "VMM", "VSC", "HGNEC", "MLAV", "PANET",
    "GEJ", "POCA", "THME", "OUTT", "CCOV", "PLLS", "CSCLC",
    "STAS", "ESCA", "MRLS", "STAD", "IHCH", "THYC", "ESCC",
    "RCSNOS", "NSGCT", "CHDM", "MUP", "HCC", "PLEMESO", "ECAD",
    "PUO", "UCS", "BRCA", "GBM", "ULMS", "IMMC", "LUNE",
    "CCS", "RCC", "MPE", "CDRCC", "SCCO", "MPNST", "SFTCNS",
    "ALUCA", "ACYC", "MDLC", "SDCA", "PLMESO", "NSCLCPD", "MACR",
    "DCIS", "APE", "NSCLC", "CLL", "OCS", "SCLC", "THYM",
    "CHOS", "CCHM", "OPHSC", "RAML", "PGNG", "GBC", "MYEC",
    "GIST", "PB", "ARMS", "SARCL", "MNG", "ACRM", "COM",
    "SBWDNET", "BCL", "HNMUCM", "SNA", "PAST", "CSCC", "SYNS",
    "SKCM", "CHS", "THFO", "GNG", "CESC", "DLBCL", "EPM",
    "ASPS", "FIBS", "PAMPCA", "FIOS", "PNET", "PRCC", "BA",
    "ACBC", "LNET", "MXOV", "SARCNOS"
)
POG_cohort_colors <- setNames(
    colorRampPalette(brewer.pal(11, "Spectral"))(length(POG_cohorts)),
    POG_cohorts
)
