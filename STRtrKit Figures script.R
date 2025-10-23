#set linkgroup= [something]

library(tidyverse)
library(ggplot2)

##read and prep data
data <- read_csv(file="STRtrKit.csv")
data$categorical_score <- factor(
  data$ScoreCat,
  levels = c(
    "Contradictory",
    "Limited",
    "Moderate",
    "Strong",
    "Supportive",
    "Definitive"
  ),
  ordered = TRUE)
unique(data$ScoreCat)
data %>% filter(
  Group %in% c(
    "Clingen",
    "Ambry",
    "Labcorp", 
    "STRtrKit",
    "G2P"
  )
) %>% group_by(
  Gene
) %>% summarise(
  num_groups_present = n_distinct (
    Group
  )
) %>% filter(
  num_groups_present >= 4
) -> shared.genes
##Jitter
ggplot(
  data = data,
  aes(
    x=Gene,
    y=Group,
    color=ScoreCat,
    size=ScoreCat
  )
) +  
  geom_jitter(
    width = 0.2,
    height = 0.2
  )
##Bar
ggplot(
  data = data,
  aes(
    x=ScoreCat,
    fill = Group
  )
) +  
  geom_bar(
    position = "dodge"
  )

##simple plots
ggplot(
  data = shared.genes, 
  aes(
    x=Group, 
    fill=Group
    )
  ) + 
  geom_bar(
    position = "dodge"
    )

ggplot(
  data = shared.genes, 
  aes(
    x=Gene
    , fill=Gene
    )
  ) + 
  geom_bar(
    position = "dodge"
    )

##UpSet Plot
library(UpSetR)
library(data.table)
dt <- fread("STRtrKit.csv")
# Make sure column names are consistent
setnames(dt, old = names(dt), new = tolower(names(dt)))
# Drop duplicates of Gene–Group pairs
dt_unique <- unique(dt[, .(gene, group)])
# Build list of genes per group
sets <- split(dt_unique$gene, dt_unique$group)
sets <- lapply(sets, unique)
# Convert to incidence matrix
incidence <- UpSetR::fromList(sets)
# Plot
upset(
  incidence,
  nsets = min(10, ncol(incidence)),   # show up to 10 groups
  nintersects = 19,                   # show top 30 intersections
  order.by = c("degree", "freq"),
  decreasing = c(TRUE, TRUE),
  empty.intersections = "on",
  mainbar.y.label = "Intersection size",
  sets.x.label = "Genes per group"
)

##Sankey Plot
library(networkD3)
install.packages("networkD3")
library(dplyr)
install.packages("dplyr")
library(htmlwidgets)
install.packages("htmlwidgets")
library(RColorBrewer)
install.packages("RColorBrewer")
library(jsonlite)
install.packages("jsonlite")
#Link 1: Gene -> Group
links1 <- aggregate(
  x = list(value = rep(1, nrow(data))),
  by = list(source = data$Gene,target = data$Group),
  FUN = length
)

#Flow 2: Group -> Score
links2 <- aggregate(
  x = list(value = rep(1, nrow(data))),
  by = list(source = data$Group,target = data$ScoreCat),
  FUN = length
)

#Combine them
links_df <- bind_rows(links1, links2)

#nodes table
nodes_df <- data.frame(
  name = unique(c(links_df$source,links_df$target)),
  stringsAsFactors = FALSE
)

#Prep the data
links_df$source_id <- match(links_df$source,nodes_df$name) - 1
links_df$target_id <- match(links_df$target,nodes_df$name) - 1
links_df$LinkGroup <- nodes_df$name[
  links_df$target_id + 1
]
nodes_df$NodeGroup <- nodes_df$name
sources <- unique(links_df$source)
targets <- unique(links_df$target)
first_target_nodes <- setdiff(targets,sources)
first_target_names <- nodes_df$name[first_target_nodes + 1]

#Colors
unique_groups <- sort(unique(first_target_names))
palette <- RColorBrewer::brewer.pal(n = min(length(unique_groups), 12), "Set3")
my_colors <- sprintf(
  'd3.scaleOrdinal().domain(%s).range(%s)',
  jsonlite::toJSON(unique_groups),
  jsonlite::toJSON(palette)
)
#Plot
sankey <- sankeyNetwork(
  Links = links_df,
  Nodes = nodes_df,
  Source = "source_id",
  Target = "target_id",
  Value = "value",
  NodeID = "name",
  LinkGroup = "LinkGroup",
  NodeGroup = "NodeGroup",
  fontSize = 12,
  nodeWidth = 30
)
sankey
saveWidget(sankey,"sankey_plot.html", selfcontained = TRUE)
