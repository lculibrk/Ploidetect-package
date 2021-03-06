compressdata <- function(t, x, segmentation_threshold = segmentation_threshold){
  # Arrange data by position and remove tibble-ness
  t <- t %>% arrange(pos) %>% as.data.frame()
  # If this is the first iteration, set npoints to one (npoints records how many original vertices were compressed into one vertex)
  if(is.null(t$npoints)){
    t$npoints <- 1
  }
  # Get differences between neighbouring points
  diffs <- abs(diff(t$residual/t$npoints))
  # Initialize graph from data
  graph <- graph(edges = c(row.names(t)[1], rep(row.names(t[-c(1, nrow(t)),]), each = 2), row.names(t)[nrow(t)]), directed = F)
  # Set edges to have the diffs attribute
  graph <- set_edge_attr(graph, name = "diff", value = diffs)
  # Give vertices appropriate attributes
  graph <- set_vertex_attr(graph, name = "residual", value = t$residual)
  graph <- set_vertex_attr(graph, name = "npoints", value = t$npoints)
  graph <- set_vertex_attr(graph, name = "from", value = t$pos)
  graph <- set_vertex_attr(graph, name = "to", value = t$end)
  # Holdover
  toy <- graph
  # loop over all vertices
  sort_by_diff <- data.frame("vertex" = V(toy)$name)
  sort_by_diff$diff <- 0
  for(row in 1:nrow(sort_by_diff)){
    sort_by_diff$diff[row] <- max(edge_attr(toy, "diff", incident(toy, sort_by_diff$vertex[row])))
  }
  sort_by_diff <- sort_by_diff %>% arrange(diff)
  for(vertex in sort_by_diff$vertex){
    if(length(incident(toy, vertex)) == 0){
      next
    }
    # If vertex is an outlier (diffs are over some threshold fraction of what we expect for a copy change) then break all edges
    if(all(edge_attr(toy, "diff", incident(toy, vertex)) > segmentation_threshold*x)){
      toy <- delete_edges(toy, incident(toy, vertex))
      next
    }
    # If vertex has two edges, break the one with larger "diff"
    if(length(incident(toy, vertex)) == 2){
      toy <- delete_edges(toy, incident(toy, vertex)[which.max(get.edge.attribute(toy, "diff", incident(toy, vertex)))])
    }
  }
  #print(get.vertex.attribute(toy, "npoints", V(toy)))
  #print(toy)
  # Get list of all vertex pairs to merge
  tomerge <- ends(toy, E(toy))
  # Get all vertices
  vertices <- V(toy)
  # Coerce vertices into a format where value is the vertex value and name is vertex name
  vertnames <- names(vertices)
  vertices <- as.numeric(vertices)
  names(vertices) <- vertnames
  # Change "tomerge" from names to values
  tomerge[,2] <- vertices[which(names(vertices) %in% tomerge[,2])]
  tomerge[,1] <- vertices[which(names(vertices) %in% tomerge[,1])]
  # Not needed I think
  #todelete <- vertices[which(vertices %in% tomerge[,2])]
  # Change pairs of vertices to repeat the same vertex twice (used in contract.vertices() to map which vertices to contract)
  vertices[which(vertices %in% tomerge[,2])] <- tomerge[,1]
  mode(vertices) <- "integer"
  lint <- vertices[1]
  for(i in 2:length(vertices)){
    if(vertices[i] - 1 > lint){
      vertices[vertices == vertices[i]] <- lint + 1
    }
    lint <- vertices[i]
  }
  
  # Merge connected vertices
  toy <- contract.vertices(toy, mapping = vertices, vertex.attr.comb = list("residual" = "sum", "npoints" = "sum", "from" = "min", "to" = "max", "name" = "first"))
  # Delete all old vertices
  #toy <- delete.vertices(toy, which(names(V(toy)) == "character(0)"))
  # Reconstruct the data.frame we began with
  dat <- data.frame("residual" = get.vertex.attribute(toy, "residual"), 
                    "npoints" = get.vertex.attribute(toy, "npoints"), 
                    "pos" = get.vertex.attribute(toy, "from"),
                    "end" = get.vertex.attribute(toy, "to"))
  #print(dat[which(dat$npoints == 1),])
  if(T){
    print(paste0("Iteration complete, segment count = ", nrow(dat)))
  }
  #print(dat)
  return(dat)
}
