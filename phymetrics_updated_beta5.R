#####Load Libraries-------------------------------------------------------------

suppressMessages(require(ape))
suppressMessages(require(tidyverse))
suppressMessages(require(phangorn))
suppressMessages(require(foreach))
suppressMessages(require(doParallel))
suppressMessages(library(getopt))
suppressMessages(library(tibble))

#####Set Up---------------------------------------------------------------------

#Options
spec <- matrix(c(
  "help",      "h", 0, "logical",    # flag to show help
  "tree",      "t", 1, "character",  # path to tree file (required)
  "metadata",  "m", 1, "character",  # path to metadata csv (required)
  "nthreads",  "n", 1, "integer",    # number of threads (optional)
  "ranks",     "k", 1, "character",  # comma-separated list of ranks (required)
  "outgroups", "o", 1, "character",  # comma-separated list of outgroup ids (required)
  "root",      "r", 1, "character",  # comma-separated list of tip names to root on (optional)
  "pre",       "p", 1, "character"   # optional prefix for output files
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

#Check Required Arguments
if (!is.null(opt$help) || is.null(opt$tree) || is.null(opt$metadata) ||
    is.null(opt$ranks) || is.null(opt$outgroups)) {
  cat(getopt(spec, usage = TRUE))
  quit(save = "no", status = 1)
}

help_message <- "
Usage: Rscript PhyMetrics.R [options]

Options:
  -t, --tree         Path to input phylogenetic tree file (Newick format) [required]
  -m, --metadata     Path to metadata CSV file with taxonomic assignments [required]
  -k, --ranks        Comma-separated list of taxonomic ranks to analyse (e.g. genus,species) [required]
  -o, --outgroups    Comma-separated list of tip IDs to define outgroups [required]
  -r, --root         Comma-separated tip names to use for rooting the tree (default: outgroups)
  -n, --threads      Number of parallel threads to use (default: number of available cores - 1)
  -p, --pre          Prefix to use for output files (optional)
  -h, --help         Show this help message and exit

Description:
  This script analyses a phylogenetic tree with associated metadata to compute
  cluster metrics such as the taxonomic retention index (tRI), taxonomic consistency index (tCI),
  taxonomic Simpson diversity index (tSDI), identifies clusters, and detects intruders.
"
if (any(c("-h", "--help") %in% commandArgs(trailingOnly = TRUE))) {
  cat(help_message)
  quit(save = "no", status = 0)
}

#Parse Arguments
tree_file <- opt$tree
metadata_file <- opt$metadata
num_threads <- ifelse(is.null(opt$nthreads), parallel::detectCores() - 1, opt$nthreads)

ranks <- strsplit(opt$ranks, ",")[[1]] %>% trimws()
outgroups <- strsplit(opt$outgroups, ",")[[1]] %>% trimws()

#If root is not provided, default to outgroup
root_tips <- if (!is.null(opt$root)) {
  strsplit(opt$root, ",")[[1]] %>% trimws()
} else {
  outgroups
}

#Make Parallel
cl <- makeCluster(num_threads)
doParallel::registerDoParallel(cl)
on.exit(stopCluster(cl))

#Load Arguments
tree <- ape::read.tree(tree_file)
tree <- root(tree, outgroup = root_tips, resolve.root = TRUE)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
colnames(metadata)[1] <- "mt_id"  # Force first column to be 'mt_id'
metadata <- metadata %>% filter(mt_id %in% tree$tip.label)
metadata[metadata == ""] <- NA
metadata[] <- lapply(metadata, trimws)



output_prefix <- ifelse(is.null(opt$pre), "", paste0(opt$pre, "_"))
#####Functions------------------------------------------------------------------


compute_indices <- function(binary_row, tree) {
  tip_labels <- colnames(binary_row)
  in_group <- tip_labels[which(binary_row == 1)]
  
  if (length(in_group) < 2 || length(in_group) >= length(tip_labels)) {
    return(list(tRI = NA, tCI = NA, parsimony_score = NA))
  }
  
  dat <- phyDat(binary_row, type = "USER", levels = c("0", "1"))
  s <- parsimony(tree, dat)
  
  g <- sum(binary_row == 1)
  m <- 1
  
  if ((g - m) == 0 || is.na(s) || s == 0) {
    return(list(tRI = NA, tCI = NA, parsimony_score = s))
  }
  
  tRI <- (g - s) / (g - m)
  tCI <- m / s
  
  return(list(tRI = tRI, tCI = tCI, parsimony_score = s))
}


get_transition_clusters <- function(tree, transition_nodes, node_states) {
  get_cluster <- function(node) {
    tips <- c()
    queue <- list(node)
    
    while (length(queue) > 0) {
      current <- queue[[1]]
      queue <- queue[-1]
      
      if (current <= length(tree$tip.label)) {
        if (node_states[[current]] == "1") {
          tips <- c(tips, tree$tip.label[current])
        }
      } else {
        children <- tree$edge[tree$edge[, 1] == current, 2]
        for (child in children) {
          if (("1" %in% node_states[[child]]) & !(child %in% transition_nodes)) {
            queue <- c(queue, child)
          }
        }
      }
    }
    unique(tips)
  }
  lapply(transition_nodes, get_cluster)
}

find_transitions <- function(tree, states) {
  tree <- reorder(tree, "postorder")
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  total_nodes <- n_tips + n_nodes
  
  node_states <- vector("list", total_nodes)
  node_states[1:n_tips] <- lapply(states[tree$tip.label], function(x) x)
  
  transitions <- list()
  
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    if (child > n_tips && is.null(node_states[[child]])) {
      children <- tree$edge[tree$edge[, 1] == child, 2]
      states1 <- node_states[[children[1]]]
      states2 <- node_states[[children[2]]]
      
      intersect_states <- intersect(states1, states2)
      if (length(intersect_states) > 0) {
        node_states[[child]] <- intersect_states
      } else {
        node_states[[child]] <- union(states1, states2)
        transitions[[length(transitions) + 1]] <- child
      }
    }
  }
  list(transitions = transitions, node_states = node_states)
}

merge_contiguous_sisters <- function(tree, clusters) {
  merged <- clusters
  changed <- TRUE
  
  parent_of <- function(node) {
    row <- which(tree$edge[,2] == node)
    if(length(row) == 0) return(NA)
    tree$edge[row, 1]
  }
  
  while(changed) {
    changed <- FALSE
    n <- length(merged)
    if(n <= 1) break
    
    to_remove <- integer(0)
    
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        
        # Representative node for cluster i
        node_i <- if(length(merged[[i]]) == 1) {
          tip_index <- which(tree$tip.label == merged[[i]])
          parent_of(tip_index)   # parent of tip
        } else {
          getMRCA(tree, merged[[i]])
        }
        
        # Representative node for cluster j
        node_j <- if(length(merged[[j]]) == 1) {
          tip_index <- which(tree$tip.label == merged[[j]])
          parent_of(tip_index)
        } else {
          getMRCA(tree, merged[[j]])
        }
        
        if(is.na(node_i) || is.na(node_j)) next
        
        # Parents and grandparents
        parent_i <- parent_of(node_i)
        parent_j <- parent_of(node_j)
        parent_2i <- if(!is.na(parent_i)) parent_of(parent_i) else NA
        parent_2j <- if(!is.na(parent_j)) parent_of(parent_j) else NA
        
        # --- Merge conditions ---
        
        # (1) Sisters: same parent
        is_sister <- !is.na(parent_i) && !is.na(parent_j) && parent_i == parent_j
        
        # (2) Directly connected: one MRCA is parent of the other
        is_directly_connected <- parent_i == node_j || parent_j == node_i
        
        # (3) Paraphyletic: one cluster’s grandparent equals the other's parent
        is_paraphyletic_cluster <- (!is.na(parent_2i) && !is.na(parent_j) && parent_2i == parent_j) ||
          (!is.na(parent_2j) && !is.na(parent_i) && parent_2j == parent_i)
        
        # (4) Paraphyletic: one cluster’s grandparent equals the other's node (tip parent)
        is_paraphyletic_tip <- (!is.na(parent_2i) && parent_2i == node_j) ||
          (!is.na(parent_2j) && parent_2j == node_i)
        
        if(is_sister || is_directly_connected || is_paraphyletic_cluster || is_paraphyletic_tip) {
          merged[[i]] <- unique(c(merged[[i]], merged[[j]]))
          to_remove <- c(to_remove, j)
          changed <- TRUE
        }
      }
    }
    
    if(length(to_remove) > 0) merged <- merged[-unique(to_remove)]
  }
  
  return(merged)
}


flag_nested_clusters <- function(tree, clusters) {
  if (length(clusters) == 0) {
    return(data.frame(
      cluster_id = integer(0),
      is_nested = logical(0),
      nested_in = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  results <- data.frame(
    cluster_id = seq_along(clusters),
    is_nested = FALSE,
    nested_in = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(clusters)) {
    focal <- clusters[[i]]
    if (length(focal) < 2) next
    mrca_node <- getMRCA(tree, focal)
    descendant_tips <- phangorn::Descendants(tree, mrca_node, type = "tips")[[1]]
    descendant_labels <- tree$tip.label[descendant_tips]
    for (j in seq_along(clusters)) {
      if (i == j) next
      other <- clusters[[j]]
      if (all(other %in% descendant_labels)) {
        results$is_nested[j] <- TRUE
        results$nested_in[j] <- i
      }
    }
  }
  return(results)
}

calc_tSDI <- function(clusters, total_terminals) {
  s <- length(clusters)
  if (s == 0 || total_terminals <= 1) return(NA_real_)
  
  ni <- sapply(clusters, length)
  numerator <- sum(ni * (ni - 1))
  denominator <- total_terminals * (total_terminals - 1)
  
  tSDI <- numerator / denominator
  return(tSDI)
}

compute_node_depths <- function(tree) {
  tree <- reorder(tree, "cladewise")

  n_nodes <- max(tree$edge)
  depth <- numeric(n_nodes)

  for(i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]
    child  <- tree$edge[i,2]
    depth[child] <- depth[parent] + 1
  }

  depth
}

build_parent_map <- function(tree) {
  parent <- integer(max(tree$edge))
  parent[tree$edge[,2]] <- tree$edge[,1]
  parent
}

fast_mrca <- function(a, b, parent, depth) {
  while(depth[a] > depth[b]) a <- parent[a]
  while(depth[b] > depth[a]) b <- parent[b]

  while(a != b) {
    a <- parent[a]
    b <- parent[b]
  }
  a
}

# distance between two nodes
dist_nodes <- function(a, b, depth, parent) {
  m <- fast_mrca(a, b, parent, depth)
  depth[a] + depth[b] - 2 * depth[m]
}

compute_tree_diameter <- function(tree, depth, parent) {
  tips <- seq_along(tree$tip.label)

  start <- tips[1]

  d1 <- sapply(tips, function(t) {
    dist_nodes(start, t, depth, parent)
  })

  tip1 <- tips[which.max(d1)]

  d2 <- sapply(tips, function(t) {
    dist_nodes(tip1, t, depth, parent)
  })

  tip2 <- tips[which.max(d2)]

  D_max <- dist_nodes(tip1, tip2, depth, parent)

  list(D_max = D_max, tip1 = tip1, tip2 = tip2)
}

# Function to calculate CSI for one group
compute_tCSI_for_group <- function(tree, clusters, depth, parent, D_max, gamma = 0.25) {

  if(length(clusters) == 0) {
    return(list(C = NA, S = NA, tCSI = NA, s = 0, N = 0))
  }

  s <- length(clusters)
  ni <- sapply(clusters, length)
  N <- sum(ni)

  if(s == 1) {
    return(list(C = 1, S = 0, tCSI = 1, s = s, N = N))
  }

  # ---- Concentration ----
  p <- ni / N
  H <- -sum(p * log(p))
  C <- 1 - (H / log(s))

  # ---- cluster nodes (MRCAs) ----
  cluster_nodes <- sapply(clusters, function(cl) {
    if(length(cl) == 1) {
      idx <- which(tree$tip.label %in% cl)
      tree$edge[match(idx, tree$edge[,2]), 1]
    } else {
      getMRCA(tree, cl)
    }
  })

  # ---- NEW separation: mean pairwise distance between clusters ----
  k <- length(cluster_nodes)
  dists <- c()

  for(i in 1:(k - 1)) {
    for(j in (i + 1):k) {
      d <- dist_nodes(cluster_nodes[i], cluster_nodes[j], depth, parent)
      dists <- c(dists, d)
    }
  }

  S <- ifelse(length(dists) > 0 && D_max > 0, mean(dists) / D_max, 0)

  # ---- Final score ----
  tCSI <- (1 - gamma) * C + gamma * (1 - S)

  list(C = C, S = S, tCSI = tCSI, s = s, N = N)
}
#####MAIN SCRIPT----------------------------------------------------------------
start_time <- Sys.time()

#####Create Presence/Absence Tables---------------------------------------------

# For each rank, generate a binary presence/absence matrix
presence_absence_matrix <- list()
for (rank in ranks) {
  meta_clean <- metadata %>%
    filter(is.na(.data[[rank]]) & mt_id %in% outgroups | !is.na(.data[[rank]])) %>%
    select(mt_id, group = all_of(rank))
  
  if (nrow(meta_clean) == 0) next
  
  rank_pa_matrix <- meta_clean %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = mt_id, values_from = present, values_fill = 0) %>%
    filter(!is.na(group)) %>%
    column_to_rownames("group")
  
  presence_absence_matrix[[rank]] <- rank_pa_matrix
}

print("Step 1: Presence/Absence matrix created")

#####Compute tRI & tCI----------------------------------------------------------

tRI_tCI_results <- foreach(rank = names(presence_absence_matrix), .combine = bind_rows, .packages = c("phangorn", "ape", "tibble", "dplyr")) %:%
  foreach(i = seq_len(nrow(presence_absence_matrix[[rank]])), .combine = bind_rows) %dopar% {
    rank_pa_matrix <- presence_absence_matrix[[rank]]
    row <- rank_pa_matrix[i, , drop = FALSE]
    group_name <- rownames(rank_pa_matrix)[i]
    valid_tips <- colnames(rank_pa_matrix)
    pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))
    
    indices_result <- compute_indices(row, pruned_tree)
    n_tips <- sum(row)
    
    if (!is.na(indices_result$tRI) && !is.na(indices_result$tCI)) {
      tibble(
        rank = rank,
        group = group_name,
        tRI = indices_result$tRI,
        tCI = indices_result$tCI,
        parsimony_score = indices_result$parsimony_score,
        n_tips = n_tips
      )
    } else {
      tibble()
    }
  }

write.csv(tRI_tCI_results, paste0(output_prefix, "tRI_tCI_per_group.csv"), row.names = FALSE)
print("Step 2: tRI and tCI calculated")

#####Compute Ensemble for Indices-----------------------------------------------

#filter outgroups
outgroup_groups_by_rank <- list()
for (rank in ranks) {
  outgroup_groups <- metadata %>%
    filter(mt_id %in% outgroups, !is.na(.data[[rank]])) %>%
    pull(!!sym(rank)) %>%
    unique()
  outgroup_groups_by_rank[[rank]] <- outgroup_groups
}

filtered_results <- tRI_tCI_results
for (rank in names(outgroup_groups_by_rank)) {
  outgroup_groups <- outgroup_groups_by_rank[[rank]]
  filtered_results <- filtered_results %>%
    filter(!(rank == rank & group %in% outgroup_groups))
}

#calculate mean
ensemble_tRI <- filtered_results %>%
  group_by(rank) %>%
  summarise(mean_tRI = mean(tRI, na.rm = TRUE), 
            mean_tCI = mean(tCI, na.rm = TRUE), 
            .groups = "drop")


write.csv(ensemble_tRI, paste0(output_prefix, "tRI_tCI_ensemble.csv"), row.names = FALSE)
rm(tRI_tCI_results, ensemble_tRI); gc()
print("Step 3: Ensemble statistics calculated")

#####Cluster Detection----------------------------------------------------------

cluster_results <- foreach(rank = names(presence_absence_matrix), .packages = c("ape", "tibble", "dplyr")) %dopar% {
  rank_pa_matrix <- presence_absence_matrix[[rank]]
  if (nrow(rank_pa_matrix) == 0) return(NULL)
  
  valid_tips <- colnames(rank_pa_matrix)
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))
  
  group_clusters <- list()
  
  for (i in seq_len(nrow(rank_pa_matrix))) {
    group_name <- rownames(rank_pa_matrix)[i]
    binary_row <- rank_pa_matrix[i, , drop = FALSE]
    
    states <- setNames(as.list(as.character(binary_row[1, ])), colnames(binary_row))
    
    res <- find_transitions(pruned_tree, states)
    transition_nodes <- res$transitions
    node_states <- res$node_states
    
    raw_clusters <- get_transition_clusters(pruned_tree, transition_nodes, node_states)
    
    tip_labels <- colnames(binary_row)
    allowed_tips <- tip_labels[which(binary_row == 1)]
    filtered_list <- lapply(raw_clusters, function(x) intersect(x, allowed_tips))
    filtered_list <- Filter(function(x) length(x) > 0, filtered_list)
    unique_clusters <- filtered_list[!duplicated(lapply(filtered_list, function(x) paste(sort(x), collapse = "_")))]
    
    if(length(unique_clusters) > 1) {
      unique_clusters <- merge_contiguous_sisters(tree, unique_clusters)
    }
    
    nested_flags <- flag_nested_clusters(pruned_tree, unique_clusters)
    
    # Merge nested clusters with their parents
    if (nrow(nested_flags) > 0 && any(nested_flags$is_nested)) {
      merged_clusters <- unique_clusters
      clusters_to_remove <- integer(0)
      
      for (j in which(nested_flags$is_nested)) {
        nested_cluster <- merged_clusters[[j]]
        parent_index <- nested_flags$nested_in[j]
        
        if (!is.na(parent_index) && parent_index <= length(merged_clusters)) {
          parent_cluster <- merged_clusters[[parent_index]]
          combined_cluster <- unique(c(parent_cluster, nested_cluster))
          merged_clusters[[parent_index]] <- combined_cluster
          clusters_to_remove <- c(clusters_to_remove, j)
        }
      }
      
      if (length(clusters_to_remove) > 0) {
        merged_clusters <- merged_clusters[-clusters_to_remove]
      }
      
      unique_clusters <- merged_clusters
    }
    
    group_clusters[[group_name]] <- list(
      clusters = unique_clusters,
      nested_flags = nested_flags
    )
  }
  return(group_clusters)
}
names(cluster_results) <- names(presence_absence_matrix)

cluster_df <- list()

for (rank in names(cluster_results)) {
  group_clusters <- cluster_results[[rank]]
  if (is.null(group_clusters)) next
  
  for (group_name in names(group_clusters)) {
    clusters <- group_clusters[[group_name]]$clusters
    if (length(clusters) == 0) next
    
    cluster_df_tmp <- tibble(
      rank = rank,
      group = group_name,
      cluster_id = seq_along(clusters),
      tips = sapply(clusters, function(x) paste(sort(x), collapse = ",")),
      n_tips = sapply(clusters, length)
    )
    
    cluster_df[[length(cluster_df) + 1]] <- cluster_df_tmp
  }
}

cluster_df <- bind_rows(cluster_df)

write.csv(cluster_df, paste0(output_prefix, "clusters.csv"), row.names = FALSE)
rm(cluster_df); gc()
print("Step 4: Clusters identified")

#####Calculate tSDI-------------------------------------------------------------

tSDI_results <- foreach(rank = names(cluster_results), .combine = bind_rows) %dopar% {
  rank_groups <- cluster_results[[rank]]  # List of groups
  
  # pa table for this rank (presence/absence)
  rank_pa_matrix <- presence_absence_matrix[[rank]]
  
  res_list <- list()
  
  for (group_name in names(rank_groups)) {
    clusters <- rank_groups[[group_name]]$clusters
    total_terminals <- sum(rank_pa_matrix[group_name, ])
    
    tSDI_val <- calc_tSDI(clusters, total_terminals)
    
    res_list[[length(res_list) + 1]] <- tibble(
      rank = rank,
      group = group_name,
      tSDI = tSDI_val,
      n_clusters = length(clusters),
      total_terminals = total_terminals
    )
  }
  
  bind_rows(res_list)
}

filtered_tSDI <- tSDI_results
for (rank_name in names(outgroup_groups_by_rank)) {
  outgroup_groups <- outgroup_groups_by_rank[[rank_name]]
  filtered_tSDI <- filtered_tSDI %>%
    filter(!(rank == rank_name & group %in% outgroup_groups))
}

ensemble_tSDI <- filtered_tSDI %>%
  group_by(rank) %>%
  summarise(mean_tSDI = mean(tSDI, na.rm = TRUE), .groups = "drop")

write.csv(tSDI_results, paste0(output_prefix, "tSDI_per_group.csv"), row.names = FALSE)
write.csv(ensemble_tSDI, paste0(output_prefix, "tSDI_ensemble.csv"), row.names = FALSE)
rm(tSDI_results, ensemble_tSDI); gc()
print("Step 5: tSDI scores computed")


#####Calculate new tCSI------------------------------------------------------
tCSI_results <- foreach(rank = names(cluster_results),
                        .combine = bind_rows,
                        .packages = c("ape","tibble")) %dopar% {

  rank_groups <- cluster_results[[rank]]
  if(is.null(rank_groups)) return(NULL)

  rank_pa <- presence_absence_matrix[[rank]]
  valid_tips <- colnames(rank_pa)
  if(length(valid_tips) == 0) return(NULL)

  # prune tree
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))

  # unit tree
  unit_tree <- pruned_tree
  unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))

  # ---- precompute structures ----
  depth  <- compute_node_depths(unit_tree)
  parent <- build_parent_map(unit_tree)

  # ONLY keep diameter for normalization
  diam <- compute_tree_diameter(unit_tree, depth, parent)
  D_max <- diam$D_max

  # ---- compute per group ----
  res_list <- vector("list", length(rank_groups))
  i <- 1

  for(group_name in names(rank_groups)) {

    clusters <- rank_groups[[group_name]]$clusters
    if(length(clusters) == 0) next

    res <- compute_tCSI_for_group(
      tree = unit_tree,
      clusters = clusters,
      depth = depth,
      parent = parent,
      D_max = D_max,
      gamma = 0.25
    )

    res_list[[i]] <- tibble(
      rank = rank,
      group = group_name,
      tCSI = res$tCSI,
      C = res$C,
      S = res$S,
      n_clusters = res$s,
      total_terminals = res$N
    )

    i <- i + 1
  }

  gc()
  bind_rows(res_list)
}

#####Filter outgroups if needed-------------------------------------------------

for(rank_name in names(outgroup_groups_by_rank)) {
  outgroup_groups <- outgroup_groups_by_rank[[rank_name]]
  tCSI_results <- tCSI_results %>%
    filter(!(rank == rank_name & group %in% outgroup_groups))
}

#####Calculate ensemble statistics----------------------------------------------

ensemble_tCSI <- tCSI_results %>%
  group_by(rank) %>%
  summarise(mean_tCSI = mean(tCSI, na.rm = TRUE),
            mean_C = mean(C, na.rm = TRUE),
            mean_S = mean(S, na.rm = TRUE),
            .groups = "drop")


write.csv(tCSI_results, paste0(output_prefix, "tCSI_per_group.csv"), row.names = FALSE)
write.csv(ensemble_tCSI, paste0(output_prefix, "tCSI_ensemble.csv"), row.names = FALSE)
rm(tCSI_results, ensemble_tCSI); gc()
print("Step 6: tCSI scores computed")

################################################################################
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(paste("Total time elapsed:", round(as.numeric(elapsed_time, units = "secs"), 2), "seconds"))
