
suppressMessages(library(igraph))

N <- 9
middle <- 3

N_color <- "#3030cc"
M_color <- "#61a347"
concept_color <- "#fff9ab"
edges_color <- "#202020"
vertex_size <- 15

outer <- (N-middle)/2
L <- list(c(1:outer), c((outer+1):(outer+middle)), c((outer+middle+1):N))
N <- sum(length(unlist(L)))
U <- lapply(L, function(u) {u+N})
M <- N

focus_one_cluster <- function() {
    edges <- c()
    N12 <- outer + middle
    for(i in 1:N12) {
        for(j in 1:N12) {
            edges <- c(edges, c(i,N12+j))
        }
    }
    graph <- make_bipartite_graph(c(rep(0, N12), rep(1, N12)), edges, directed = TRUE)
    add_vertices(graph, 2*outer)
    locations1 <- t(sapply(1:N12, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:N12, function(j) { c(j,1) }))
    locations3 <- t(sapply(((N12+1):N), function(j){ c(j,0) }))
    locations4 <- t(sapply(((N12+1):M), function(j){ c(j,1) }))
    locations <- rbind(locations1, locations2, locations3, locations4)
    vertex_colors <- c(rep(N_color, N12), rep(M_color, N12), rep(N_color, outer), rep(M_color, outer))
    return(list(graph, locations, vertex_colors))
}

focus_one_cluster_binding <- function() {
    binding <- N+M+1
    edges <- c()
    N12 <- outer + middle
    for(i in 1:N12) {
        l <- i
        edges <- c(edges, c(l, binding))
    }
    for(j in 1:N12) {
        u <- j + N12
        edges <- c(edges, c(binding, u))
    }
    graph <- make_graph(edges, directed = TRUE)
    add_vertices(graph, 2*outer)
    locations1 <- t(sapply(1:N12, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:N12, function(j) { c(j,1) }))
    locations3 <- t(sapply(((N12+1):N), function(j){ c(j,0) }))
    locations4 <- t(sapply(((N12+1):M), function(j){ c(j,1) }))
    locations <- rbind(locations1, locations2, locations3, locations4, get_location("12"))
    vertex_colors <- c(rep(N_color, N12), rep(M_color, N12), rep(N_color, outer), rep(M_color, outer), concept_color)
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

single_binding <- function(label=c("12", "23", "mid_L", "mid_U", "none")) {
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
            if(abs(I-J) > 1) {
                next
            }
            edges <- c(edges, get_all_edges(L[[J]], U[[I]]))
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
    if(label == "none") {
        locations <- rbind(locations1, locations2)
        vertex_colors <- c(rep(N_color, N), rep(M_color, M))
    } else {
        locations <- rbind(locations1, locations2, get_location(label=label))
        vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 1))
    }
    graph <- make_graph(edges, directed = TRUE)
    return(list(graph, locations, vertex_colors))
}

double_binding <- function(labels=c("12", "23", "mid_L", "mid_U", "none")) {
    edges <- c()
    additionals <- c(N+M+1, N+M+2)

    for(I in 1:3) {
        for(J in 1:3) {
            if("12" %in% labels && I %in% c(1,2) && J %in% c(1,2)) {
                next
            }
            if("23" %in% labels && I %in% c(2,3) && J %in% c(2,3)) {
                next
            }
            if("mid_L" %in% labels && J == 2) {
                next
            }
            if("mid_U" %in% labels && I == 2) {
                next
            }
            if(abs(I-J) > 1) {
                next
            }
            edges <- c(edges, get_all_edges(L[[J]], U[[I]]))
        }
    }

    count = 1
    for(label in labels) {
        additional <- additionals[count]
        count <- count + 1
        if(label == "12") {
            for(I in 1:2) {
                if(labels[count] == "mid_U" && I == 2) {
                    next
                }
                for(u in U[[I]]) {
                    edges <- c(edges, c(additional, u))
                }
            }
            for(J in 1:2) {
                if(labels[count] == "mid_L" && J == 2) {
                    next
                }
                for(l in L[[J]]) {
                    edges <- c(edges, c(l, additional))
                }
            }
        }
        if(label == "23") {
            for(I in 2:3) {
                if(labels[count] == "mid_U" && I == 2) {
                    next
                }
                for(u in U[[I]]) {
                    edges <- c(edges, c(additional, u))
                }
            }
            for(J in 2:3) {
                if(labels[count] == "mid_L" && J == 2) {
                    next
                }
                for(l in L[[J]]) {
                    edges <- c(edges, c(l, additional))
                }
            }
        }
        if(label == "mid_L") {
            for(l in L[[2]]) {
                edges <- c(edges, c(l, additional))
            }
            Uset <- unlist(U)
            if(labels[count-2] == "12") {
                Uset <- U[[3]]
            }
            if(labels[count-2] == "23") {
                Uset <- U[[1]]
            }
            for(u in Uset) {
                edges <- c(edges, c(additional, u))
            }
        }
        if(label == "mid_U") {
            if(labels[count-2] == "12") {
                Lset <- L[[3]]
            }
            if(labels[count-2] == "23") {
                Lset <- L[[1]]
            }
            for(l in Lset) {
                edges <- c(edges, c(l, additional))
            }
            for(u in U[[2]]) {
                edges <- c(edges, c(additional, u))
            }
        }
    }
    if(labels[2] %in% c("mid_L")) {
        additionals <- c(additionals[2], additionals[1])
    }
    edges <- c(edges, additionals)

    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))
    locations <- rbind(locations1, locations2, get_location(label=labels[1]), get_location(label=labels[2]))
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 2))
    graph <- make_graph(edges, directed = TRUE)
    return(list(graph, locations, vertex_colors))
}

triple_binding <- function(missing=c("12", "23", "mid_L", "mid_U")) {
    edges <- c()
    all4 <- c("12", "23", "mid_L", "mid_U")
    remaining <- setdiff(all4, missing)
    additionals <- c(N+M+1, N+M+2, N+M+3)

    for(I in 1:3) {
        for(J in 1:3) {
            if("12" %in% remaining && I %in% c(1,2) && J %in% c(1,2)) {
                next
            }
            if("23" %in% remaining && I %in% c(2,3) && J %in% c(2,3)) {
                next
            }
            if("mid_L" %in% remaining && J == 2) {
                next
            }
            if("mid_U" %in% remaining && I == 2) {
                next
            }
            if(abs(I-J) > 1) {
                next
            }
            edges <- c(edges, get_all_edges(L[[J]], U[[I]]))
        }
    }

    ranks <- c(0, 0, -1, 1)
    names(ranks) <- all4
    for(i in 1:3) {
        for(j in 1:3) {
            ci <- remaining[i]
            cj <- remaining[j]
            if(ranks[ci] > ranks[cj] && abs(ranks[ci]-ranks[cj]) < 2) {
                edges <- c(edges, c(additionals[j], additionals[i]))
            }
        }
    }

    names(additionals) <- remaining
    for(label in remaining) {
        if(label == "mid_U") {
            for(u in U[[2]]) {
                edges <- c(edges, c(additionals[label], u))
            }
        }
        if(label == "mid_L") {
            for(l in L[[2]]) {
                edges <- c(edges, c(l, additionals[label]))
            }
        }
        if(label == "12") {
            for(l in L[[1]]) {
                edges <- c(edges, c(l, additionals[label]))
            }
            for(u in U[[1]]) {
                edges <- c(edges, c(additionals[label], u))
            }
        }
        if(label == "23") {
            for(l in L[[3]]) {
                edges <- c(edges, c(l, additionals[label]))
            }
            for(u in U[[3]]) {
                edges <- c(edges, c(additionals[label], u))
            }
        }
    }

    if(missing == "mid_U") {
        for(u in c(U[[2]])) {
            edges <- c(edges, c(additionals["12"], u))
        }
        for(u in c(U[[2]])) {
            edges <- c(edges, c(additionals["23"], u))
        }
    }
    if(missing == "mid_L") {
        for(l in c(L[[2]])) {
            edges <- c(edges, c(l, additionals["12"]))
        }
        for(l in c(L[[2]])) {
            edges <- c(edges, c(l, additionals["23"]))
        }
    }
    if(missing == "12") {
        for(u in c(U[[1]])) {
            edges <- c(edges, c(additionals["mid_L"], u))
        }
        for(l in c(L[[1]])) {
            edges <- c(edges, c(l, additionals["mid_U"]))
        }
    }
    if(missing == "23") {
        for(u in c(U[[3]])) {
            edges <- c(edges, c(additionals["mid_L"], u))
        }
        for(l in c(L[[3]])) {
            edges <- c(edges, c(l, additionals["mid_U"]))
        }
    }

    locations1 <- t(sapply(1:N, function(j) { c(j,0) }))
    locations2 <- t(sapply(1:M, function(j) { c(j,1) }))
    locations <- rbind(locations1, locations2, get_location(remaining[1]), get_location(remaining[2]), get_location(remaining[3]))
    vertex_colors <- c(rep(N_color, N), rep(M_color, M), rep(concept_color, 3))
    graph <- make_graph(edges, directed = TRUE)
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
    plot(l[[1]], layout=l[[2]], vertex.label=NA, vertex.color=l[[3]], edge.color=edges_color, vertex.size=vertex_size, edge.arrow.size=0.2)
}

to_png <- function(l, name) {
    png(name)
    par(mar=c(0,0,0,0)+.1)
    plot_result(l)
}

graph_focus_one_cluster <- focus_one_cluster()
graph_focus_one_cluster_binding <- focus_one_cluster_binding()
graph_bicluster <- single_binding("none")
graph_single_12 <- single_binding("12")
graph_single_23 <- single_binding("23")
graph_single_mid_L <- single_binding("mid_L")
graph_single_mid_U <- single_binding("mid_U")
graph_12_mid_U <- double_binding(labels=c("12", "mid_U"))
graph_23_mid_U <- double_binding(labels=c("23", "mid_U"))
graph_12_mid_L <- double_binding(labels=c("12", "mid_L"))
graph_23_mid_L <- double_binding(labels=c("23", "mid_L"))
graph_triple_12 <- triple_binding("12")
graph_triple_23 <- triple_binding("23")
graph_triple_mid_L <- triple_binding("mid_L")
graph_triple_mid_U <- triple_binding("mid_U")
graph_bind_all <- bind_all()

# to_png(graph_focus_one_cluster, "graph_focus_one_cluster.png")
# to_png(graph_focus_one_cluster_binding, "graph_focus_one_cluster_binding.png")
# to_png(graph_bicluster, "graph_bicluster.png")
# to_png(graph_single_12, "graph_single_12.png")
# to_png(graph_single_23, "graph_single_23.png")
# to_png(graph_single_mid_L, "graph_single_mid_L.png")
# to_png(graph_single_mid_U, "graph_single_mid_U.png")
# to_png(graph_bind_all, "graph_bind_all.png")

dev.new()
par(mfrow = c(3, 3), mai=c(0.1,0.1,0.1,0.1))
plot_result(graph_12_mid_U)
plot_result(graph_single_mid_U)
plot_result(graph_23_mid_U)
plot_result(graph_single_12)
plot_result(graph_bicluster)
plot_result(graph_single_23)
plot_result(graph_12_mid_L)
plot_result(graph_single_mid_L)
plot_result(graph_23_mid_L)

dev.new()
par(mfrow = c(3, 3), mai=c(0.1,0.1,0.1,0.1))
plot.new()
plot_result(graph_triple_mid_L)
plot.new()
plot_result(graph_triple_23)
plot_result(graph_bind_all)
plot_result(graph_triple_12)
plot.new()
plot_result(graph_triple_mid_U)
plot.new()

