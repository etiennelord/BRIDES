
info_network<-function(g1,g2=NULL, attributes='') {	
	if (is.null(g2)) g2=g1
	if (attributes=='') {
		#we take all the node node in graph1			
		vertex_of_g2<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;			
	} else {
		vertex_of_g2=V(g2)[V(g2)$tax==as.factor(attributes)]$name;		
	}	
	vertex_of_g1=length(V(g1)$name)
	cat("Network characteristics:\n");
	if (attributes=='') {
		cat("Total of new nodes in network Y:", length(vertex_of_g2),"\n"); 
	} else {
		str1=paste('Selected nodes "',attributes,'" in network Y   :');
		cat(str1, length(vertex_of_g2),"\n"); 
	}
	path_to_investigate=(vertex_of_g1*(vertex_of_g1-1))/2;
	if (is.directed(g1)||is.directed(g2)) {
	path_to_investigate=(vertex_of_g1*vertex_of_g1)-vertex_of_g1;
	}
	
	cat("Number of edges in network Y:", length(E(g2)), "\n");
	cat("Number of nodes in network Y:", length(V(g2)$name), "\n"); 
	cat("Number of nodes in network X:", length(V(g1)$name), "\n"); 
	cat("Total of pathways to investigate:", path_to_investigate,"\n");
	cat("Clustering coefficient network Y:", transitivity(g2),"\n");
	cat("Clustering coefficient network X:", transitivity(g1),"\n");
	cat("Average degree \u00B1 std in network Y:", mean(degree(g2)),"\u00B1",sd(degree(g2)),"\n");
	cat("Average degree \u00B1 std in network X:", mean(degree(g1)),"\u00B1",sd(degree(g1)),"\n");
	cat("Average path length in network Y:",average.path.length(g2), "\n");	
	cat("Average path length in network X:",average.path.length(g1), "\n");
	cat("Number of clusters in network Y:",clusters(g2)$no, "\n");
	cat("Number of clusters in network X:",clusters(g1)$no, "\n");
	cat("Average cluster size \u00B1 std in network Y:",mean(clusters(g2)$csize),"\u00B1",sd(clusters(g2)$csize),"\n");	
	cat("Average cluster size \u00B1 std in network X:",mean(clusters(g1)$csize),"\u00B1",sd(clusters(g1)$csize),"\n");	
	#cat("Nodes distribution in network Y (first row taxa, second row count):\n");
	#print(summary(as.factor(V(g2)$tax)));
}

