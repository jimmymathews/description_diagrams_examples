
library(igraph)

N <- 10
M <- 10

edges <- c()
for(i in 1:N) {
    for(j in 1:M) {
        edges <- c(edges, c(i,N+j))
    }
}
g1 <- make_bipartite_graph(c(rep(0, N), rep(1, M)), edges, directed = TRUE)
locs1 <- t(sapply(1:N, function(j) { c(j,0) }))
locs2 <- t(sapply(1:M, function(j) { c(j,1) }))
locs_g1 <- rbind(locs1, locs2)
png("g1.png")
par(mar=c(0,0,0,0)+.1)
plot(g1, layout=locs_g1, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M)), edge.color="#202020", vertex.size=20)
dev.off()



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
locs1 <- t(sapply(1:N, function(j) { c(j,0) }))
locs2 <- t(sapply(1:M, function(j) { c(j,1) }))
locsb <- c((1+N)/2, 0.5)
locs_g2 <- rbind(locs1, locs2, locsb)
g2 <- make_graph(edges, directed = TRUE)
png("g2.png")
par(mar=c(0,0,0,0)+.1)
plot(g2, layout=locs_g2, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M), rep("#CC3030", 4)), edge.color="#202020", vertex.size=20)
dev.off()


edges <- c()
# U <- list(c(1,2), c(3,4,5), c(6,7))
# L <- lapply(U, function(u) {u+7})
L <- list(c(1,2,3), c(4,5,6,7), c(8,9,10))
U <- lapply(L, function(u) {u+N})
bindingL <- matrix(0, nrow=4, ncol=4)
bindingU <- matrix(0, nrow=4, ncol=4)
bindingL[2,4] <- N+M+1
bindingL[1,2] <- N+M+2
bindingL[2,3] <- N+M+3
bindingU[2,4] <- N+M+4
bindingU[1,2] <- N+M+2
bindingU[2,3] <- N+M+3
for(I in 1:3) {
    for(J in 1:3) {
        if(bindingU[I,J] != 0) {
            UIJ <- setdiff(union(U[[I]], U[[J]]), U[[2]])
            for(u in UIJ) {
                edges <- c(edges, c(bindingU[I,J], u))
            }
        }
        if(bindingL[I,J] != 0) {
            LIJ <- setdiff(union(L[[I]], L[[J]]), L[[2]])
            for(l in LIJ) {
                edges <- c(edges, c(l, bindingL[I,J]))
            }
        }
    }
    if(I == 2) {
        for(u in U[[I]]) {
            edges <- c(edges, c(bindingU[2,4], u))
        }
        for(l in L[[I]]) {
            edges <- c(edges, c(l, bindingL[2,4]))
            c(l, bindingL[2,4])
        }
    }
}
edges <- c(edges, c(bindingL[2,4], bindingL[1,2]))
edges <- c(edges, c(bindingL[2,4], bindingL[2,3]))
edges <- c(edges, c(bindingU[1,2], bindingU[2,4]))
edges <- c(edges, c(bindingU[2,3], bindingU[2,4]))

locs1 <- t(sapply(1:N, function(j) { c(j,0) }))
locs2 <- t(sapply(1:M, function(j) { c(j,1) }))
locsL2 <- c(  (1+N)/2, 0.25)
locs12 <- c(  (1+7)/2, 0.50)
locs23 <- c( (4+10)/2, 0.50)
locsU2 <- c(  (1+N)/2, 0.75)
locs_g3 <- rbind(locs1, locs2, locsL2, locs12, locs23, locsU2)
g3 <- make_graph(edges, directed = TRUE)
png("g3.png")
par(mar=c(0,0,0,0)+.1)
plot(g3, layout=locs_g3, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M), rep("#CC3030", 4)), edge.color="#202020", vertex.size=20)
dev.off()

dev.new()
plot(g1, layout=locs_g1, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M)), edge.color="#202020", vertex.size=20)
dev.new()
plot(g2, layout=locs_g2, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M), rep("#CC3030", 4)), edge.color="#202020", vertex.size=20)
dev.new()
plot(g3, layout=locs_g3, vertex.label=NA, vertex.color=c(rep("#3030CC", N), rep("#61a347", M), rep("#CC3030", 4)), edge.color="#202020", vertex.size=20)





