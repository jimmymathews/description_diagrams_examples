
library(igraph)

N <- 13
middle <- 5
outer <- (N-middle)/2
L <- list(c(1:outer), c((outer+1):(outer+middle)), c((outer+middle+1):N))
N <- sum(length(unlist(L)))
U <- lapply(L, function(u) {u+N})
M <- N

N_color <- "#3030cc"
M_color <- "#61a347"
concept_color <- "#fff9ab"
edges_color <- "#202020"
vertex_size <- 10

big_cluster <- function() {
    edges <- c()
    for(i in 1:N) {
        for(j in 1:M) {
            edges <- c(edges, c(i,N+j))
        }
    }
    graph <- make_bipartite_graph(c(rep(0, N), rep(1, M)), edges, directed = TRUE)
    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))
    locations <- rbind(locations1, locations2)
    vertex_colors <- c(rep(N_color, N), rep(M_color, M))
    return(list(graph, locations, vertex_colors))
}

big_cluster_binding <- function() {
    binding <- N+M+1
    edges <- c()
    for(i in 1:N) {
        l <- i
        edges <- c(edges, c(l, binding))
    }
    for(j in 1:M) {
        u <- j + N
        edges <- c(edges, c(binding, u))
    }
    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))
    locationsb <- c((1+N)/2, 0.5)
    locations <- rbind(locations1, locations2, locationsb)
    graph <- make_graph(edges, directed = TRUE)
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 4))
    return(list(graph, locations, vertex_colors))
}

get_location <- function(label=c("12", "23", "mid_L", "mid_U")) {
    d <- ( (middle+outer-1) ) /2
    i <- (middle-1)
    a <- (outer)
    e <- i / (2*(a+i))

    if(label == "12") {
        return(c( 1 +         d, 0.50))
    }
    if(label == "23") {
        return(c( 1 + (N-1) - d, 0.50))
    }
    if(label == "mid_L") {
        return(c( 1 + (N-1)/2, e))
    }
    if(label == "mid_U") {
        return(c( 1 + (N-1)/2, 1-e))   
    }
}

single_binding <- function(label=c("12", "23", "mid_L", "mid_U")) {
    edges <- c()

    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))

    locations <- rbind(locations1, locations2, ...)
    graph <- make_graph(edges, directed = TRUE)
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 1))
    return(list(graph, locations, vertex_colors))    
}

bind_all <- function() {
    pair_binding_L <- matrix(0, nrow=3, ncol=3)
    pair_binding_U <- matrix(0, nrow=3, ncol=3)
    pair_binding_L[1,2] <- N+M+1
    pair_binding_U[1,2] <- N+M+1
    pair_binding_L[2,3] <- N+M+2
    pair_binding_U[2,3] <- N+M+2
    mid_L <- N+M+3
    mid_U <- N+M+4

    edges <- c()

    for(I in 1:3) {
        for(J in 1:3) {
            if(pair_binding_U[I,J] != 0) {
                UIJ <- setdiff(union(U[[I]], U[[J]]), U[[2]])
                for(u in UIJ) {
                    edges <- c(edges, c(pair_binding_U[I,J], u))
                }
            }
            if(pair_binding_L[I,J] != 0) {
                LIJ <- setdiff(union(L[[I]], L[[J]]), L[[2]])
                for(l in LIJ) {
                    edges <- c(edges, c(l, pair_binding_L[I,J]))
                }
            }
        }
    }

    for(I in 1:3) {
        if(I == 2) {
            for(u in U[[I]]) {
                edges <- c(edges, c(mid_U, u))
            }
            for(l in L[[I]]) {
                edges <- c(edges, c(l, mid_L))
                c(l, mid_L)
            }
        }
    }

    edges <- c(edges, c(mid_L,               pair_binding_L[1,2]))
    edges <- c(edges, c(mid_L,               pair_binding_L[2,3]))
    edges <- c(edges, c(pair_binding_U[1,2], mid_U))
    edges <- c(edges, c(pair_binding_U[2,3], mid_U))

    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))

    # d <- ( (middle+outer-1) ) /2
    # locations12 <- c( 1 +         d, 0.50)
    # locations23 <- c( 1 + (N-1) - d, 0.50)

    # i <- (middle-1)
    # a <- (outer)
    # e <- i / (2*(a+i))
    # locationsL2 <- c( 1 + (N-1)/2, e)
    # locationsU2 <- c( 1 + (N-1)/2, 1-e)
    locations <- rbind(locations1, locations2, get_location("12"), get_location("23"), get_location("mid_L"), get_location("mid_U"))
    graph <- make_graph(edges, directed = TRUE)
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 4))
    return(list(graph, locations, vertex_colors))
}

plot_result <- function(l) {
    plot(l[[1]], layout=l[[2]], vertex.label=NA, vertex.color=l[[3]], edge.color=edges_color, vertex.size=vertex_size)    
}

graph_all <- big_cluster()
png("graph_all.png")
par(mar=c(0,0,0,0)+.1)
plot_result(graph_all)
dev.off()

graph_big_cluster_binding <- big_cluster_binding()
png("graph_big_cluster_binding.png")
par(mar=c(0,0,0,0)+.1)
plot_result(graph_big_cluster_binding)
dev.off()

graph_bind_all <- bind_all()
png("graph_bind_all.png")
par(mar=c(0,0,0,0)+.1)
plot_result(graph_bind_all)
dev.off()

# Generates PDFs when used with Rscript
# dev.new()
# plot_result(graph_all)
# dev.new()
# plot_result(graph_big_cluster_binding)
# dev.new()
# plot_result(graph_bind_all)





