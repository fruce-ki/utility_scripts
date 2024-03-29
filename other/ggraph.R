library(ggraph)
library(tidygraph)

set_graph_style(plot_margin = margin(1,1,1,1))


# https://ggraph.data-imaginist.com/articles/Layouts.html
# https://ggraph.data-imaginist.com/articles/Nodes.html
# https://ggraph.data-imaginist.com/articles/Edges.html

# Layout is a projection of the nodes onto X-Y coordinates
# All layouts from the >> graphlayouts << and >> igraph << packages are available, and ggraph itself also provides some of the more specialised layouts itself.

graph <- as_tbl_graph(highschool)

# Not specifying the layout - defaults to "auto"
ggraph(graph) + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()

ggraph(graph, layout = 'kk') + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()

ggraph(graph, layout = 'kk', maxiter = 100) + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()

layout <- create_layout(graph, layout = 'eigen'). # really just a df
head(layout)
attributes(layout)
ggraph(layout) + 
  geom_edge_link(aes(colour = factor(year))) + 
  geom_node_point()




# An arc diagram
ggraph(graph, layout = 'linear') + 
  geom_edge_arc(aes(colour = factor(year)))
# A coord diagram
ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(colour = factor(year))) + 
  coord_fixed()


graph <- tbl_graph(flare$vertices, flare$edges)
# An icicle plot
ggraph(graph, 'partition') + 
  geom_node_tile(aes(fill = depth), size = 0.25)
# A sunburst plot
ggraph(graph, 'partition', circular = TRUE) + 
  geom_node_arc_bar(aes(fill = depth), size = 0.25) + 
  coord_fixed()



graph <- as_tbl_graph(highschool) %>% 
  mutate(degree = centrality_degree())
lapply(c('stress', 'fr', 'lgl', 'graphopt'), function(layout) {
  ggraph(graph, layout = layout) + 
    geom_edge_link(aes(colour = factor(year)), show.legend = FALSE) +
    geom_node_point() + 
    labs(caption = paste0('Layout: ', layout))
})



graph <- graph %>% 
  mutate(friends = ifelse(
    centrality_degree(mode = 'in') < 5, 'few',
    ifelse(centrality_degree(mode = 'in') >= 15, 'many', 'medium')
  ))
ggraph(graph, 'hive', axis = friends, sort.by = degree) + 
  geom_edge_hive(aes(colour = factor(year))) + 
  geom_axis_hive(aes(colour = friends), size = 2, label = FALSE) + 
  coord_fixed()



ggraph(graph, 'focus', focus = node_is_center()) + 
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), data.frame(r = 1:5), colour = 'grey') + 
  geom_edge_link() + 
  geom_node_point() + 
  coord_fixed()



graph <- tbl_graph(flare$vertices, flare$edges)
set.seed(1)
ggraph(graph, 'circlepack', weight = size) + 
  geom_node_circle(aes(fill = depth), size = 0.25, n = 50) + 
  coord_fixed()

set.seed(1)
ggraph(graph, 'circlepack', weight = size) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = depth)) +
  coord_fixed()



ggraph(graph, 'treemap', weight = size) + 
  geom_node_tile(aes(fill = depth), size = 0.25)

ggraph(graph, 'treemap', weight = size) + 
  geom_edge_link() + 
  geom_node_point(aes(colour = depth))



ggraph(graph, 'cactustree') + 
  geom_node_circle(aes(fill = depth), size = 0.25) + 
  coord_fixed()

importFrom <- match(flare$imports$from, flare$vertices$name)
importTo <- match(flare$imports$to, flare$vertices$name)

ggraph(graph, 'cactustree') + 
  geom_node_circle(aes(fill = depth), 
                   size = 0.25, 
                   alpha = 0.2
  ) + 
  geom_conn_bundle(aes(colour = after_stat(index)),
                   data = get_con(importFrom, importTo),
                   edge_alpha = 0.25
  ) +
  theme(legend.position = "none") +
  coord_fixed()


ggraph(graph, 'tree') + 
  geom_edge_diagonal()

dendrogram <- hclust(dist(iris[, 1:4]))
ggraph(dendrogram, 'dendrogram', height = height) + 
  geom_edge_elbow()

ggraph(dendrogram, 'dendrogram', circular = TRUE) + 
  geom_edge_elbow() + 
  coord_fixed()

tree <- create_tree(100, 2, directed = FALSE) %>% 
  activate(edges) %>% 
  mutate(length = runif(n()))
ggraph(tree, 'unrooted', length = length) + 
  geom_edge_link()


graph <- create_notable('zachary')
ggraph(graph, 'matrix', sort.by = node_rank_leafsort()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed()

ggraph(graph, 'matrix', sort.by = node_rank_spectral()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed()



ggraph(graph, 'fabric', sort.by = node_rank_fabric()) + 
  geom_node_range(colour = 'grey') + 
  geom_edge_span(end_shape = 'square') + 
  coord_fixed()

ggraph(graph, 'fabric', sort.by = node_rank_fabric(), shadow.edges =TRUE) + 
  geom_node_range(colour = 'grey') + 
  geom_edge_span(aes(filter = shadow_edge), colour ='lightblue' , end_shape = 'square') + 
  geom_edge_span(aes(filter = !shadow_edge), end_shape = 'square') + 
  coord_fixed()





graph <- create_notable('meredith') %>% 
  mutate(group = sample(c('A', 'B'), n(), TRUE))

ggraph(graph, 'stress') + 
  geom_node_voronoi(aes(fill = group), max.radius = 1) + 
  geom_node_point() + 
  geom_edge_link() + 
  coord_fixed()



gr <- tbl_graph(flare$vertices, flare$edges)

ggraph(gr, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal() + 
  geom_node_point(aes(filter = leaf)) + 
  coord_fixed()



l <- ggraph(gr, layout = 'partition', circular = TRUE)

l + geom_node_arc_bar(aes(fill = depth)) + 
  coord_fixed()

l + geom_edge_diagonal() + 
  geom_node_point(aes(colour = depth)) + 
  coord_fixed()






library(purrr)
library(rlang)

set_graph_style(plot_margin = margin(1,1,1,1))
hierarchy <- as_tbl_graph(hclust(dist(iris[, 1:4]))) %>% 
  mutate(Class = map_bfs_back_chr(node_is_root(), .f = function(node, path, ...) {
    if (leaf[node]) {
      as.character(iris$Species[as.integer(label[node])])
    } else {
      species <- unique(unlist(path$result))
      if (length(species) == 1) {
        species
      } else {
        NA_character_
      }
    }
  }))

hairball <- as_tbl_graph(highschool) %>% 
  mutate(
    year_pop = map_local(mode = 'in', .f = function(neighborhood, ...) {
      neighborhood %E>% pull(year) %>% table() %>% sort(decreasing = TRUE)
    }),
    pop_devel = map_chr(year_pop, function(pop) {
      if (length(pop) == 0 || length(unique(pop)) == 1) return('unchanged')
      switch(names(pop)[which.max(pop)],
             '1957' = 'decreased',
             '1958' = 'increased')
    }),
    popularity = map_dbl(year_pop, ~ .[1]) %|% 0
  ) %>% 
  activate(edges) %>% 
  mutate(year = as.character(year))



ggraph(hairball, layout = 'stress') + 
  geom_edge_link(aes(colour = year))

ggraph(hairball, layout = 'stress') + 
  geom_edge_fan(aes(colour = year))

# let's make some of the student love themselves
loopy_hairball <- hairball %>% 
  bind_edges(tibble::tibble(from = 1:5, to = 1:5, year = rep('1957', 5)))

ggraph(loopy_hairball, layout = 'stress') + 
  geom_edge_link(aes(colour = year), alpha = 0.25) + 
  geom_edge_loop(aes(colour = year))

ggraph(hairball, layout = 'stress') + 
  geom_edge_density(aes(fill = year)) + 
  geom_edge_link(alpha = 0.25)



ggraph(hierarchy, layout = 'dendrogram', height = height) + 
  geom_edge_elbow()

ggraph(hierarchy, layout = 'dendrogram', height = height) + 
  geom_edge_bend()

ggraph(hierarchy, layout = 'dendrogram', height = height) + 
  geom_edge_diagonal()



ggraph(hierarchy, layout = 'dendrogram', height = height) + 
  geom_edge_elbow2(aes(colour = node.Class))



# Random names - I swear
simple <- create_notable('bull') %>% 
  mutate(name = c('Thomas', 'Bob', 'Hadley', 'Winston', 'Baptiste')) %>% 
  activate(edges) %>% 
  mutate(type = sample(c('friend', 'foe'), 5, TRUE))

ggraph(simple, layout = 'graphopt') + 
  geom_edge_link(arrow = arrow(length = unit(4, 'mm'))) + 
  geom_node_point(size = 5)

ggraph(simple, layout = 'graphopt') + 
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(size = 5)

ggraph(simple, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name),
                     end_cap = label_rect(node2.name)), 
                 arrow = arrow(length = unit(4, 'mm'))) + 
  geom_node_text(aes(label = name))


ggraph(simple, layout = 'graphopt') + 
  geom_edge_link(aes(label = type), 
                 angle_calc = 'along',
                 label_dodge = unit(2.5, 'mm'),
                 arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(size = 5)









library(data.table)



E <- rbind(data.table( expand.grid(c('a', 'b', 'c', 'd'), c('a', 'e')) ),
           data.table( expand.grid(c('k'), c('g', 'h', 'i', 'j')) ) )
setnames(E, c('from', 'to'))
l = union(levels(E$from), levels(E$to))
#E[, from := factor(from, levels = l)]
#E[, to := factor(to, levels = l)]
E[, from := as.character(from)]
E[, to := as.character(to)]
E <- E[from != to, ]
E[, score := sample.int(n = 50, size = .N, replace = TRUE)]

V <- data.table(node = 1:length(l), tag = l)


gdt <- tbl_graph(nodes = V, edges = E, node_key = 'tag')
groups <- igraph::components(gdt)  


ggraph(gdt, layout = 'tree') + 
  geom_edge_link(aes(colour = score,
                     start_cap = label_rect(node1.tag),
                     end_cap = label_rect(node2.tag))) + 
  geom_node_text(aes(label = tag,
                     colour = as.factor(groups$membership)))


