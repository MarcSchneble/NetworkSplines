# adjust path
setwd("~/NetworkSplines")

# load required packages
library(igraph)
library(colorspace)
library(graphics)
library(latex2exp)


# definition of considered graph by connected edges
A = matrix(c(1,4,1,2,2,3,4,3,4,7,3,7,4,6,6,7,6,5), byrow = T, ncol = 2)
G = graph_from_edgelist(A, directed = FALSE)


# graph representation
pdf(file = "Plots/ExampleNetworkGraph.pdf", width = 5, height = 5)
par(mar = c(0, 0, 1.25, 0))
set.seed(6)
plot.igraph(G, vertex.color = "black", vertex.size = 10, vertex.frame.color = "black",
            vertex.label.font = 4, vertex.label.cex = 2, vertex.label.color = "black",
            vertex.label.degree = c(-pi/4, 0.8*pi, pi/2, -pi/4, -pi/6, -pi/6, pi),
            vertex.label = c(TeX("$\\mathbf{v_1}$"), TeX("$\\mathbf{v_2}$"), TeX("$\\mathbf{v_3}$"),
                             TeX("$\\mathbf{v_4}$"), TeX("$\\mathbf{v_5}$"), TeX("$\\mathbf{v_6}$"), TeX("$\\mathbf{v_7}$")),
            edge.color = rainbow_hcl(9), edge.width = 5,
            edge.label = c(TeX("$\\mathbf{e_1}$"), TeX("$\\mathbf{e_2}$"), TeX("$\\mathbf{e_3}$"), 
                           TeX("$\\mathbf{e_4}$"), TeX("$\\mathbf{e_5}$"), TeX("$\\mathbf{e_6}$"),
                           TeX("$\\mathbf{e_7}$"), TeX("$\\mathbf{e_8}$"), TeX("$\\mathbf{e_9}$")),
            edge.label.color = rainbow_hcl(9), edge.label.cex = 2,
            edge.label.x = c( 0.59,  0.71, -0.32, -0.15, -0.48, -0.97, -0.12, -0.84, -0.69),
            edge.label.y = c(-0.31, -0.91, -0.88, -0.41, -0.21, -0.36,  0.22,  0.30,  0.71),
            vertex.label.dist = c(2,2,2,2,2,2,2.5),
            layout = layout_nicely, margin = c(0, 0.1, 0, 0.2))
title(main = "Network Graph Representation", cex.main = 1.25)
dev.off()


# network representation
v = matrix(c(0.2,  0.6, 1.3,  0.0, 0.35,  0,  1,
             1.5, 1.40, 1.25, 1.0,  0.6,  0,  0), ncol = 2)

pdf(file = "Plots/ExampleNetworkGeometric.pdf", width = 5, height = 5)
par(mar = c(0, 0, 1.25, 0))
plot.igraph(G, vertex.size = 0, vertex.shape = "none",
            vertex.label.font = 4, vertex.label.cex = 2, vertex.label.color = "black",
            vertex.label.degree = c(pi, pi/2, -pi/6, -pi, -pi/2, pi/6, pi/6),
            vertex.label.dist = 1.5,
            vertex.label = c(TeX("$\\mathbf{v_1}$"), TeX("$\\mathbf{v_2}$"), TeX("$\\mathbf{v_3}$"),
                             TeX("$\\mathbf{v_4}$"), TeX("$\\mathbf{v_5}$"), TeX("$\\mathbf{v_6}$"), TeX("$\\mathbf{v_7}$")),
            edge.color = rainbow_hcl(9), edge.width = 5, 
            edge.label = c(TeX("$\\mathbf{e_1}$"), TeX("$\\mathbf{e_2}$"), TeX("$\\mathbf{e_3}$"), 
                           TeX("$\\mathbf{e_4}$"), TeX("$\\mathbf{e_5}$"), TeX("$\\mathbf{e_6}$"),
                           TeX("$\\mathbf{e_7}$"), TeX("$\\mathbf{e_8}$"), TeX("$\\mathbf{e_9}$")), edge.curved = TRUE,
            edge.label.cex = 2, edge.label.color = rainbow_hcl(9),
            edge.label.x = c(-0.82, -0.41, 0.45,  0.10, -0.09,  0.94, -0.88, -0.21, -0.39),
            edge.label.y = c( 0.68,  0.96, 0.87,  0.09, -0.12, -0.25, -0.31, -0.82, -0.43),
            layout = v, margin = c(0.1, 0.1, 0, 0.25))
#title(main = "Geometric Network Representation", cex.main = 1.75)
dev.off()