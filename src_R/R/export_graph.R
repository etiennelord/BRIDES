export_graph<-function(g,file,attributes=F){
	if (file.exists(file)) file.remove(file)
	file_attr=paste(file,".attr.txt",sep="")
	file.create(file)
	t=get.edgelist(g)
	for (i in 1:nrow(t)) {
		 	cat(t[i,1],t[i,2],"\n",sep="\t",file=file, append=T);
	}
	if (attributes) {
		if (file.exists(file_attr)) file.remove(file_attr)
		file.create(file_attr)
		for (i in 1:length(V(g))) {
			cat(V(g)[i]$name,V(g)[i]$tax,"\n",sep="\t",file=file_attr, append=T);
		}
	}
}