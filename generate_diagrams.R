
suppressMessages(library(igraph))

N <- 9
middle <- 3

N_color <- "#3030cc"
M_color <- "#61a347"
concept_color <- "#fff9ab"
edges_color <- "#202020"
vertex_size <- 10

outer <- (N-middle)/2
L <- list(c(1:outer), c((outer+1):(outer+middle)), c((outer+middle+1):N))
N <- sum(length(unlist(L)))
U <- lapply(L, function(u) {u+N})
M <- N

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

get_all_edges <- function(set1, set2) {
    edges <- c()
    for(s1 in set1) {
        for(s2 in set2) {
            edges <- c(edges, c(s1,s2))
        }
    }
    return(edges)
}

single_binding <- function(label=c("12", "23", "mid_L", "mid_U")) {
    edges <- c()
    additional <- N+M+1

    for(I in 1:3) {
        for(J in 1:3) {
            if(label == "12" && I %in% c(1,2) && J %in% c(1,2)) {
                next
            }
            if(label == "23" && I %in% c(2,3) && J %in% c(2,3)) {
                next
            }
            if(label == "mid_L" && J == 2) {
                next
            }
            if(label == "mid_U" && I == 2) {
                next
            }
            edges <- c(edges, get_all_edges(U[[I]], L[[J]]))
        }
    }

    if(label == "12") {
        for(I in 1:2) {
            for(u in U[[I]]) {
                edges <- c(edges, c(additional, u))
            }
        }
        for(J in 1:2) {
            for(l in L[[J]]) {
                edges <- c(edges, c(l, additional))
            }
        }
    }
    if(label == "23") {
        for(I in 2:3) {
            for(u in U[[I]]) {
                edges <- c(edges, c(additional, u))
            }
        }
        for(J in 2:3) {
            for(l in L[[J]]) {
                edges <- c(edges, c(l, additional))
            }
        }
    }
    if(label == "mid_L") {
        for(l in L[[2]]) {
            edges <- c(edges, c(l, additional))
        }
        for(u in unlist(U)) {
            edges <- c(edges, c(additional, u))
        }
    }
    if(label == "mid_U") {
        for(l in unlist(L)) {
            edges <- c(edges, c(l, additional))
        }
        for(u in U[[2]]) {
            edges <- c(edges, c(additional, u))
        }
    }

    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))
    locations <- rbind(locations1, locations2, get_location(label=label))
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

    locations <- rbind(locations1, locations2, get_location("12"), get_location("23"), get_location("mid_L"), get_location("mid_U"))
    graph <- make_graph(edges, directed = TRUE)
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 4))
    return(list(graph, locations, vertex_colors))
}

plot_result <- function(l) {
    plot(l[[1]], layout=l[[2]], vertex.label=NA, vertex.color=l[[3]], edge.color=edges_color, vertex.size=vertex_size)    
}

to_png <- function(l, name) {
    png(name)
    par(mar=c(0,0,0,0)+.1)
    plot_result(l)
    dev.off()
}

graph_all <- big_cluster()
graph_big_cluster_binding <- big_cluster_binding()
graph_single_12 <- single_binding("12")
graph_single_23 <- single_binding("23")
graph_single_mid_L <- single_binding("mid_L")
graph_single_mid_U <- single_binding("mid_U")
graph_bind_all <- bind_all()

to_png(graph_all, "graph_all.png")
to_png(graph_big_cluster_binding, "graph_big_cluster_binding.png")
to_png(graph_single_12, "graph_single_12.png")
to_png(graph_single_23, "graph_single_23.png")
to_png(graph_single_mid_L, "graph_single_mid_L.png")
to_png(graph_single_mid_U, "graph_single_mid_U.png")
to_png(graph_bind_all, "graph_bind_all.png")

# Generates PDFs when used with Rscript
# dev.new()
# plot_result(graph_all)
# dev.new()
# plot_result(graph_big_cluster_binding)
# dev.new()
# plot_result(graph_bind_all)





