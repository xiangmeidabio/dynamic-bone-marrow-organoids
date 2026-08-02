#!/usr/bin/env Rscript

# CellChat comparison of static and dynamic BMO conditions
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. split the annotated BMO object into static and dynamic day-25 groups;
#   2. construct one CellChat object per group using the human database;
#   3. estimate and aggregate ligand-receptor communication probabilities; and
#   4. compare interaction counts, strengths, centrality, and pathways.

project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

library(CellChat)
library(Seurat)

BMO_combined_path <- file.path(project_root, "data", "processed", "BMO_combined.rds")
if (!file.exists(BMO_combined_path)) {
  stop("Missing annotated BMO object: ", BMO_combined_path)
}
BMO_combined <- readRDS(BMO_combined_path)
unique(BMO_combined$group)

Static.25d <- subset(BMO_combined, group=='Static.25d')

Static.25d.input <- LayerData(Static.25d, layer = 'data')
Static.25d.meta <- Static.25d@meta.data[,c("group","celltype")]
colnames(Static.25d.meta) <-  c("group","labels")
identical(colnames(Static.25d.input),rownames(Static.25d.meta))

Static.25d.cellchat <- createCellChat(object = Static.25d.input, meta = Static.25d.meta, group.by = "labels")

levels(Static.25d.cellchat@idents)
groupSize <- as.numeric(table(Static.25d.cellchat@idents))

CellChatDB <- CellChatDB.human
Static.25d.cellchat@DB <- CellChatDB

Static.25d.cellchat <- subsetData(Static.25d.cellchat)
Static.25d.cellchat <- identifyOverExpressedGenes(Static.25d.cellchat)
Static.25d.cellchat <- identifyOverExpressedInteractions(Static.25d.cellchat)

Static.25d.cellchat <- computeCommunProb(Static.25d.cellchat, type = "triMean")
Static.25d.cellchat <- filterCommunication(Static.25d.cellchat, min.cells = 10)
Static.25d.cellchat <- computeCommunProbPathway(Static.25d.cellchat)
Static.25d.cellchat <- aggregateNet(Static.25d.cellchat)

saveRDS(Static.25d.cellchat, file = file.path(project_root, "data", "processed", "Static.25d.cellchat.rds"))
Static.25d.cellchat <- readRDS(file.path(project_root, "data", "processed", "Static.25d.cellchat.rds"))

Dynamic.25d <- subset(BMO_combined, group=='Dynamic.25d')

Dynamic.25d.input <- LayerData(Dynamic.25d, layer = 'data')
Dynamic.25d.meta <- Dynamic.25d@meta.data[,c("group","celltype")]
colnames(Dynamic.25d.meta) <-  c("group","labels")
identical(colnames(Dynamic.25d.input),rownames(Dynamic.25d.meta))

Dynamic.25d.cellchat <- createCellChat(object = Dynamic.25d.input, meta = Dynamic.25d.meta, group.by = "labels")

levels(Dynamic.25d.cellchat@idents)
groupSize <- as.numeric(table(Dynamic.25d.cellchat@idents))

CellChatDB <- CellChatDB.human
Dynamic.25d.cellchat@DB <- CellChatDB

Dynamic.25d.cellchat <- subsetData(Dynamic.25d.cellchat)
Dynamic.25d.cellchat <- identifyOverExpressedGenes(Dynamic.25d.cellchat)
Dynamic.25d.cellchat <- identifyOverExpressedInteractions(Dynamic.25d.cellchat)

Dynamic.25d.cellchat <- computeCommunProb(Dynamic.25d.cellchat, type = "triMean")
Dynamic.25d.cellchat <- filterCommunication(Dynamic.25d.cellchat, min.cells = 10)
Dynamic.25d.cellchat <- computeCommunProbPathway(Dynamic.25d.cellchat)
Dynamic.25d.cellchat <- aggregateNet(Dynamic.25d.cellchat)

saveRDS(Dynamic.25d.cellchat, file = file.path(project_root, "data", "processed", "Dynamic.25d.cellchat.rds"))
Dynamic.25d.cellchat <- readRDS(file.path(project_root, "data", "processed", "Dynamic.25d.cellchat.rds"))


#merge cellchat obj of different group
object.list <- list(Static = Static.25d.cellchat, Dynamic = Dynamic.25d.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


#figure1
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#figure2
pdf(file = file.path(project_root, "results", "figures", "Diff_Interaction_Network.pdf"), width = 14, height = 7)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

dev.off()

#figure3
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


#figure4
object.list <- lapply(object.list, function(x) {
  x <- netAnalysis_computeCentrality(x, slot.name = "netP")
  return(x)
})

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i],
                                               weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

cellchat@netP$Static$pathways
cellchat@netP$Dynamic$pathways

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave(file.path(project_root, "results", "figures", "compare_pathway_strength.pdf"), p, width = 10, height = 12)


#figure5

netAnalysis_signalingChanges_scatter(cellchat,
                                     idents.use = "Endothelial Cells",
                                     signaling.exclude = "NOTCH",
                                     comparison = c(2, 1))


netVisual_bubble(cellchat, sources.use = "Endothelial Cells",
                 targets.use = c("HSC", "Erythroid","Megakaryocyte"),
                 comparison = c(1, 2), angle.x = 45)

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_cellchat.txt")
)
