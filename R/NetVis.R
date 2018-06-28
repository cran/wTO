#' @title NetVis
#' @param Node.1 Names of the Nodes.1 that are connected to the Nodes.2. It's the output from wTO.Complete or Consensus.
#' @param Node.2 Names of the Nodes.2 that are connected to the Nodes.1. It's the output from wTO.Complete or Consensus.
#' @param wTO weight of the links, the wTO output from wTO.Complete or wTO.Consensus.
#' @param pval p-values for the wTO value. By default it is NULL.
#' @param padj Adjusted p-values for the wTO value. By default it is NULL.
#' @param cutoff It's a list containing the kind of cutoff to be used (pval, Threshold or pval.adj)and it's value. Example: cutoff= list(kind = "Threshold", value = 0.5)
#' @param layout a layout from the igraph package.
#' @param smooth.edges If the edges should be smoothed or not.
#' @param path If the graph should be saved specify the name of the file.
#' @param Cluster TRUE or FALSE if the nodes should be clustered (double click to uncluster).
#' @param MakeGroups algorithm to find clusters. One of the followings: walktrap, optimal, spinglass, edge.betweenness, fast_greedy, infomap, louvain, label_prop, leading_eigen. Default to FALSE.
#' @param legend TRUE or FALSE if the legend should appear.
#' @param manipulation TRUE or FALSE if the graph should be editable.
#' @param shape a list shape=list(shape = "triangle", names = NULL), with the shape and the IDs that should have a different shape, shape can be: diamond, star, triangle, triangleDown or square.
#' @description Given a set of Nodes and the weight of the edges, a cutoff for the edges, it draws the networks. Returns a list with the nodes and edges  attributes. And plots the network.
#' @importFrom  visNetwork visNetwork visInteraction visEdges visOptions visClusteringByGroup visLegend visPhysics visIgraphLayout visOptions visSave visExport
#' @importFrom plyr arrange join
#' @importFrom igraph graph_from_data_frame degree E
#' @importFrom data.table as.data.table
#' @importFrom  magrittr "%>%"
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>

#' @export
#'
#' @examples
#'
#'  X =  wTO.Complete( k =1, n = 5, Data = Microarray_Expression1,
#'  Overlap = ExampleGRF$x, method = "p",  plot = FALSE)
#' # Plot with the default aguments.
#'  NetVis(Node.1 = X$wTO$Node.1, Node.2 = X$wTO$Node.2,
#'  wTO = X$wTO$wTO_sign, cutoff = list(kind =
#' "Threshold", value = 0.50))
#'
#'\dontrun{
#' # Plotting just the edges with p-value < 0.05, with straight edges, nodes clustered,
#' # no legend and mapipulation of the graph enabled.
#'   NetVis(Node.1 = X$wTO$Node.1, Node.2 = X$wTO$Node.2,
#'  wTO = X$wTO$wTO_sign, pval = X$wTO$pval_sign,
#'  padj = X$wTO$pval_sign,
#'   cutoff= list(kind = "pval", value = 0.05),
#'   smooth.edges = FALSE,
#'  Cluster = TRUE, legend = FALSE, manipulation = TRUE)
#' # Plotting just the edges with wTO > 0.50, no legend and the nodes:
#' # "ZNF738", "ZNF677" with triagle shape,
#' # no legend and mapipulation of the graph enabled.
#'  NetVis(Node.1 = X$wTO$Node.1, Node.2 = X$wTO$Node.2,
#'  wTO = X$wTO$wTO_sign, pval = X$wTO$pval_sign,
#'  padj = X$wTO$pval_sign, cutoff= list(kind = "Threshold", value = 0.5),legend = FALSE,
#'  shape = list(shape = "triangle", names = c("ZNF738", "ZNF677")))
#'
#'  }

NetVis =
  function (Node.1, Node.2, wTO, pval= NULL, MakeGroups = FALSE,
            padj= NULL,  cutoff = list(kind = "Threshold", value = 0.5),
            layout = NULL, smooth.edges = T, path = NULL, Cluster = F,
            legend = T, shape=list(shape = "triangle", names= NULL), manipulation = F)
  {
    input_vis = data.frame (Node.1 = Node.1, Node.2 = Node.2, wTO = as.numeric(wTO))
    if(!is.null(pval)){
      input_vis$pval = pval
    }
    if(!is.null(padj) ){
      input_vis$padj = padj
    }
    `%ni%` <- Negate(`%in%`)
    `%>%` <- magrittr::`%>%`
    if (cutoff$kind %ni% c("Threshold", "pval", "pval.adj")) {
      stop("cutoff kind must be \"Threshold\", \"pval\" or \"pval.adj\".")
    }
    if (is.numeric(cutoff$value) == F) {
      stop("cutoff value must be numeric.")
    }
    MakeGroups_pos = c('walktrap',
                       'optimal', 'spinglass', 'edge.betweenness',
                       'fast_greedy', 'infomap', 'louvain', 'label_prop',
                       'leading_eigen', FALSE)

    if(MakeGroups %ni% MakeGroups_pos){
      stop("MakeGroups should be FALSE or one of the following options: 'walktrap',
           'optimal', 'spinglass', 'edge.betweenness',
           'fast_greedy', 'infomap', 'louvain', 'label_prop',
           'leading_eigen'.")
    }
    if (Cluster %ni% c(T, F)) {
      stop("Cluster must be T / F.")
    }
    if (smooth.edges %ni% c(T, F)) {
      stop("smooth.edges must be T / F.")
    }
    input_vis = subset(input_vis, abs(input_vis$wTO) > 0.01)
    if (cutoff$kind == "Threshold") {
      input_vis = subset(input_vis, abs(input_vis$wTO) >= cutoff$value)
    }
    else if (cutoff$kind == "pval") {
      input_vis = subset(input_vis, input_vis$pval <= cutoff$value)
    }
    else if (cutoff$kind == "pval.adj") {
      input_vis = subset(input_vis, input_vis$pval.adj <= cutoff$value)
    }
    if (nrow(input_vis) <= 2) {
      stop("There is less than 2 nodes on your network. Choose a lower cutoff.")
    }
    if (smooth.edges == T) {
      smooth.edges = "enabled"
    }
    input_vis = input_vis[!is.na(input_vis$wTO), ]
    input_vis = plyr::arrange(input_vis, input_vis$Node.1, input_vis$Node.2)
    nodes <- data.frame(id = sort(unique(c(as.character(input_vis$Node.1),
                                           as.character(input_vis$Node.2)))))
    g = igraph::graph_from_data_frame(input_vis, directed = F)
    DEGREE = as.data.frame(igraph::degree(g))
    igraph::E(g)$weight = abs(input_vis$wTO)
    names(DEGREE) = "degree"
    DEGREE$id = row.names(DEGREE)
    nodes = suppressMessages(plyr::join(nodes, DEGREE))

    nodes$shape = ifelse(nodes$id %in% shape$names, shape$shape, "dot")

    nodes$value = (nodes$degree - min(nodes$degree))/(max(nodes$degree) -
                                                        min(nodes$degree))
    nodes$value = nodes$value * 2 + 1
    nodes$size = nodes$value

    if(MakeGroups == FALSE){
      group = 1
    }
    if (MakeGroups == 'infomap'){

      group = igraph::cluster_infomap(g)$membership
    }
    else if (MakeGroups == 'walktrap'){
      group = igraph::cluster_walktrap(g)$membership
    }
    else if (MakeGroups == 'leading_eigen'){
      group = igraph::cluster_leading_eigen(g)$membership
    }
    else if (MakeGroups == 'louvain'){
      group = igraph::cluster_louvain(g)$membership
    }
    else if (MakeGroups == 'label_prop'){
      group = igraph::cluster_label_prop(g)$membership
    }
    else if (MakeGroups == 'fast_greedy'){
      group = igraph::cluster_fast_greedy(g)$membership
    }
    else if (MakeGroups == 'optimal'){
      group = igraph::cluster_optimal(g)$membership
    }
    else if (MakeGroups == 'spinglass'){
      group = igraph::cluster_spinglass(g)$membership
    }
    else if (MakeGroups == 'edge.betweenness'){
      group = igraph::edge.betweenness.community(g)$membership
    }
    nodes = plyr::join(nodes, data.frame(id = igraph::V(g)$name, group = group))

    nodes$label = nodes$id
    # print(nodes)
    # nodes$shape = "circle"
    #

    # nodes$value = ifelse(nodes$id %in% shape$names, 3, 1)

    nodes$title = paste0("<p> Node ID: ", nodes$id, "<br>Degree: ",
                         nodes$degree, "</p>")
    edges <- data.frame(from = input_vis$Node.1, to = input_vis$Node.2)
    wto = abs(input_vis$wTO)
    edges$width = 0.5 + 5 * abs((wto - min(wto))/(max(wto) -
                                                    min(wto)))
    edges$color = ifelse(input_vis$wTO > 0, "violetred", "springgreen")
    edges$title = paste0("<p> wTO: ", round(input_vis$wTO, 2),
                         "</p>")
    ledges <- data.frame(color = c("violetred", "springgreen"),
                         label = c("+ wTO", "- wTO"), arrows = c("", ""))

    network <- visNetwork::visNetwork(nodes, edges) %>%
      visNetwork::visInteraction(navigationButtons = TRUE) %>%
      visNetwork::visEdges(smooth = smooth.edges) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                     degree = 1, hover = T), nodesIdSelection = list(enabled = TRUE,
                                                                                                     style = "width: 200px; height: 26px;\n   background: #f8f8f8;\n   color: darkblue;\n   border:none;\n   outline:none;"),
                             manipulation = F)  %>%
      visNetwork::visPhysics(enabled = F) %>%
      visNetwork::visExport(type = "pdf",
                            name = "networkpdf",
                            float = "left",
                            label = "Save pdf",
                            background = "transparent",
                            style= "")%>%
    visNetwork::visExport(type = "png",
                          name = "networkpng",
                          float = "right",
                          label = "Save png",
                          background = "transparent",
                          style= "")
    if (Cluster == T) {
      network <- network %>% visNetwork::visClusteringByGroup(groups = unique((nodes$group)))
    }
    if (legend == T) {
      network <- network %>% visNetwork::visLegend(width = 0.3,
                                                   position = "right", main = "Group", addEdges = ledges,
                                                   ncol = 2)

    }
    if (!is.null(layout)) {
      network <- network %>% visNetwork::visIgraphLayout(layout = layout)
    }
    if(manipulation == T){
      network <- network %>% visNetwork::visOptions(manipulation = TRUE)
    }
    if (is.null(path)) {
      network
    }  else if (!is.null(path)) {
      visNetwork::visSave(network, file = path)
      message(path)
    }
    nodesout =data.frame(id = nodes$id, group = nodes$group, degree = nodes$degree)
    return(list(Nodes = nodesout,  network = network))

    }
