install.packages("BiocManager")
install.packages("plotrix")
install.packages("gridExtra")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db") # Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
install.packages("ggpubr")
install.packages("mutSignatures")
install.packages("cowplot")
install.packages("ggrepel")
BiocManager::install("fgsea")
install.packages("ggplotify")
BiocManager::install("Gviz")

