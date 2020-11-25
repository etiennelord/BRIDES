#======================================================================================
#| Program : BRIDES 2020                                                              |
#| Authors : Etienne Lord (UdeM),Vladimir Makarenkov (Universite du Quebec a Montreal)|
#|      and Cindy Bouchard, Francois-Joseph Lapointe (Universite de Montreal)         |
#| This program implements a fast algorithm for characterizing evolving               |
#| similarity networks using Breakthroughs, Roadblocks, Impasses, Detours,            |
#| Equals and Shortcuts.                                                              |
#| 2019. Addition of adaptive search for best paths                                   | 
#| 2020. Addition of weighted scenarios                                               |
#| Revision : 29 jan. 2019                                                            |
#======================================================================================
#
# Usage: 
# source("BRIDES_v2020.R")
# set.seed(1)
# g=random_network(10)
# BRIDES(g$g1,g$g2)
# BRIDES(g$g1,g$g2,runmode="incremental",min_additional = 1, max_additional = 2)
# BRIDES(g$g1,g$g2,runmode="exhaustive", max_additional = 2)
# BRIDES(g$g1,g$g2,runmode="exhaustive", max_additional = 2, wt=c(-3,-1,1,1,1,1))
# 
# *Note: wt refers to the B R I D E S scoring scheme.

library(igraph)
library(SDDE)
	

  ###################################################
  # Function match_matrix
  # Helper function for genetic search
  # Giving a chromosome, look if it is contained 
  # in a given matrix (matr)
  match_matrix=function(chromosome, matr) {
    # Don't allow empty chromosome
    if (sum(chromosome)==0) return (TRUE)
    for (i in 1:nrow(matr)) {
      if (all(matr[i,]==chromosome)) return(TRUE)
    }
    return (FALSE)
  }
  
   ###################################################
  # Function allsolution
  # if some results are equals, return the combn results 
  # BRIDES 2020  
  # size is still the number of added nodes
  allsolution=function(size=2, evalFunc=NULL, maxgenes=0) {
	if (is.null(evalFunc)) {
      warning("No eval function defined. Note, we try to minimize this function.");
      return(NULL);
    }
	
    if (maxgenes>size||maxgenes==0) {
      maxgenes=size;
    }
    
    results=list()
    # population = chromosomes note: unique
    # size = size of chromosome
    # iters = number of iteration
    # mutationChance = ...
    # evaluations = call of each chromosome
    # best_chromosome = best chromosome for each round
    # final = best chromosome solution
    # final_score = score for the final function
    # final_chromosome = final chromosome index
    # brides = brides results
    # additional_nodes = additional_nodes for each chromosome
    
    results$population=matrix(0,nrow=0,ncol=size)
    results$brides=matrix(0,nrow=0, ncol=8)
    results$size=size
    results$best=c()
    results$evaluations=c()
    results$additional_nodes=c()
    
	chromosome_comb=combn(size,maxgenes)
	# Test each chromosomes
	cat("Exhaustive search for ",ncol(chromosome_comb)," iterations with max.:",maxgenes," additional nodes.\n",sep="")
	for (idx in 1:ncol(chromosome_comb)){
		best_parents=which(results$evaluations==min(results$evaluations))
		chromosome=array(0,size)
		ch_template=chromosome_comb[,idx]
		chromosome[ch_template]=1
		res=evalFunc(chromosome)
		score=res$score
		cat("[Iteration ",idx,"/",ncol(chromosome_comb),"] (score:",score,")\n", sep="")
		results$brides=rbind(results$brides,res$brides)
		results$population=rbind(results$population,chromosome)
		results$evaluations=c(results$evaluations,score)
		results$additional_nodes=c(results$additional_nodes,res$nodes)
		results$best_chromosome=c(results$best_chromosome,best_parents)
	}
    rownames(results$population)=1:nrow(results$population)
    best_parents=which(results$evaluations==min(results$evaluations))
    results$final=results$population[best_parents,]
    results$final_score=results$evaluations[best_parents]
    results$final_chromosome=best_parents;
    cat("Exhaustive search done.\n")
    return(results)
  
  }
  
  ###################################################
  # Function optimize_ga
  # Genetic Algorithm to optimize the solution of 
  # BRIDES 2019
  # This is a 0,1 chromosome model which try to 
  # MINIMIZE the evalFunc.
  # size:  is the size(len) of the binary array
  # iters: number of iterations to run
  # mutationChance: chance of having a completely random chromosome
  optimize_ga=function(size=2,iters=50,evalFunc=NULL, mutationChance=0.10, maxgenes=0, mingenes=1) {
    if (is.null(evalFunc)) {
      warning("No eval function defined. Note, we try to minimize this function.");
      return(NULL);
    }
	if (mingenes<=0) mingenes=1
    if (maxgenes>size||maxgenes==0) {
      maxgenes=size;
    }
    if (mingenes>=maxgenes-1) {
      mingenes=maxgenes-1
    }
    results=list()
    # population = chromosomes note: unique
    # size = size of chromosome
    # iters = number of iteration
    # mutationChance = ...
    # evaluations = call of each chromosome
    # best_chromosome = best chromosome for each round
    # final = best chromosome solution
    # final_score = score for the final function
    # final_chromosome = final chromosome index
    # brides = brides results
    # additional_nodes = additional_nodes for each chromosome
    
    total=iters+size # Total chromosome
    results$population=matrix(0,nrow=0,ncol=size)
    results$brides=matrix(0,nrow=0, ncol=8)
    results$size=size
    results$iters=iters
    results$mutationChance=mutationChance
    results$best=c()
    results$evaluations=c()
    results$additional_nodes=c()
    
    ## Pretraining i.e. test all unique chromosomes
    cat("Pre-training. Evaluating ",size," solutions.\n", sep="")
    for (i in 1:size) {
      chromosome=array(0,size)
      chromosome[i]=1
      res=evalFunc(chromosome)
      score=res$score
      results$brides=rbind(results$brides,res$brides)
      results$population=rbind(results$population,chromosome)
      results$evaluations=c(results$evaluations,score)
      results$additional_nodes=c(results$additional_nodes,res$nodes)
    }
    ## Training
    i=1
    retry=0
    cat("Optimization for ",iters," iterations with min.:",mingenes," max.:",maxgenes," additional nodes.\n",sep="")
    while (i<=iters) {
      best_parents=which(results$evaluations==min(results$evaluations))
      other_parents=c(1:nrow(results$population))
      other_parents=other_parents[-best_parents]
	  if (length(other_parents)==0) other_parents=best_parents; #Just in case there is no winner
      #print(best_parents)
      #print(other_parents)
      p1=sample(best_parents,1)      
	  p2=sample(other_parents,1)
      ## mix (additive)
      #print(p1)
      #print(p2)
      chromosome=results$population[p1,]|results$population[p2,]
      # generate a complete random chromosome
      if (sample(1:100,1)<(mutationChance*100)) {
        # Constrained mutations
		possible=rep(0,size)
		possible[1:mingenes]=1
		possible[mingenes+1:maxgenes]=sample(0:1,maxgenes-mingenes,replace=T)
		chromosome=sample(possible,size,replace=F)
		#chromosome=sample(0:1,size,replace = T) # Completely random chromosome
      }
      if (retry>4*iters) {
        parents=sample(other_parents,2,replace=F)
        chromosome=results$population[parents[1],]|results$population[parents[2],]
      }
      ## Test is contained already
      ##sum(chromosome)<=maxgenes&&sum(chromosome)>mingenes&&
      if (sum(chromosome)<=maxgenes&&sum(chromosome)>mingenes&&!match_matrix(chromosome,results$population)) {
        
        res=evalFunc(chromosome)
        score=res$score
        results$brides=rbind(results$brides,res$brides)
        results$additional_nodes=c(results$additional_nodes,res$nodes)
        cat("[Iteration ",i,"/",iters,"] (score:",score,")\n", sep="")
        results$population=rbind(results$population,chromosome)
        results$evaluations=c(results$evaluations,score)
        results$best_chromosome=c(results$best_chromosome,best_parents)
        i=i+1
        retry=0
      } else {
        retry=retry+1
      }
      if (retry>5*iters) {
        cat("No more suitable solutions found.\n")
        break;
      }
    }
    rownames(results$population)=1:nrow(results$population)
    best_parents=which(results$evaluations==min(results$evaluations))
    results$final=results$population[best_parents,]
    results$final_score=results$evaluations[best_parents]
    results$final_chromosome=best_parents;
    cat("Optimization done.\n")
    return(results)
    
  }
	
	###################################################
	#function split_sample
	#Split a sample x into equals parts of maxsize
	split_sample<-function(x, maxsize=1000) {
		#note: since we want both node, we multiply by 2
		maxsize<-maxsize*2;
		return (split(x, ceiling(seq_along(x)/maxsize)));
	}
	
	create_sample_path <- function(total_n, offset, size,directed=FALSE){
			if (!directed) 
				return(.C("createSamplePath", as.integer(total_n), as.integer(offset), as.integer(size), path = integer(size*2), as.integer(directed))$path);
			if (directed) {
				 path=c();
				 t=0
				 for(j in 1:total_n) {
					for(i in 1:total_n){
						if (i!=j) {			
							if (as.integer((t/size))+1==offset) {			
								path=c(path,i);
								path=c(path,j);
							}
							t=t+1
						}
						if ((as.integer(t/size)+1)>offset) break;
					}
				}
				return (path)
			}
		}

	###################################################
	#function multicore	
	# This detect the number of available cores in a 
	# system. On Windows, return 1.
		multicore<- function(nc=0) {
			cores <- if (.Platform$OS.type == "windows")
				1
			else
				min(8L, ceiling(detectCores()/2))
			getOption("mc.cores", cores)
			if (nc!=0) return (nc);
			return (cores)
		}
	###################################################
	#function to export networks to file
	#Added support for edge.weights
	export_network<-function(g,file,attributes=F){
		if (!is.igraph(g)) {
			export_network(g$g1,paste(file,"_g1.txt",sep=""), attributes);
			export_network(g$g2,paste(file,"_g2.txt",sep=""),file, attributes);
		} else {
			if (file.exists(file)) file.remove(file)
			file_attr=paste(file,".attr.txt",sep="")
			file.create(file)
			if (is.weighted(g)) {
				ta=data.frame(get.edgelist(g), E(g)$weight)
			} else {
				ta=data.frame(get.edgelist(g))
			}
			cat(nrow(ta));
			write.table(ta,sep="\t",file=file, row.names = F,col.names = F,quote = F);
			gor=g
			v_attribute=names(vertex_attr(g))
			if (attributes&& ("tax"%in%v_attribute)) {
					g=delete.vertices(g,degree(g)==0) # remove degree 0
					if (file.exists(file_attr)) file.remove(file_attr)
					file.create(file_attr)
					df <-data.frame(V(g)$name,V(g)$tax) # Will fail if vertex with degree 0
					write.table(df,sep="\t",file=file_attr, row.names = F,col.names = F,quote = F);		
				}
			# Add degree 0 nodes
			g=gor
			for (n in names(V(g)[degree(g)==0])) cat(n,"\t\n", file=file, sep="",append=T); #Slow
			if (attributes&& ("tax"%in%v_attribute)) {
				for (nid in V(g)[degree(g)==0]) {
					cat(V(g)[nid]$name,"\t",V(g)[nid]$tax,"\t\n",sep="", file=file_attr,sep="",append=T);
				}
			}
		}
	}		
	###################################################
	#function load network from edgelist
	load_network<-function(filename_or_df, filename_tax_or_df='', edge_weight='equal',directed=FALSE) {	
		if (!is.data.frame(filename_or_df)) {
			dS=read.table(filename_or_df,sep='\t');
		} else {
			dS=filename_or_df;
		}
		gS=graph.data.frame(dS,directed=directed);
		if (filename_tax_or_df!='') {
			if (!is.data.frame(filename_tax_or_df)) {
				dS2=read.table(filename_tax_or_df,sep='\t',row.names=1);
			} else {
				dS2=filename_tax_or_df;
			}
			V(gS)$tax=as.character(dS2[V(gS)$name,1]) 
		}
		#equal, proportional, or inverse
		if (edge_weight=='equal') E(gS)$weight=1.0;
		if (edge_weight==''||edge_weight=='proportional') {
			E(gS)$weight=dS[,3];		
		}
		if (edge_weight=='inverse') {
			E(gS)$weight=dS[,3];	
			E(gS)$weight=1/E(gS)$weight;
		}
		#Remove edge and node without names 
		gS=delete.vertices(gS,V(gS)[V(gS)$name==""]);
		
		return (gS);
	}

###################################################
#This is the main compute function for one pathway
pathBRIDES<-function(g1,g2_without_k,g1names,g2_unique_names_primed,g2_unique_number_primed,node1_number,node2_number,no_new_node,maxdistance,maxtime,maxnode,t0) {
	options(warn=-1); #disable warnings since some vertex could become unreachable
	print("HELLO\n")
	###################################################
    ## Flags
	verbose<-FALSE;
	trace<-FALSE;
	directed<-is.directed(g1)
	###################################################
    ## constant 
  
	c_equal=5;
	c_shortcut=6;
	c_detour=4;
	c_roadblock=2;
	c_impasse=3;
	c_breakthrough=1;
	#c_dead_end_or_detour=99;
	
	###################################################
    ## Trace type
	
	trace_len=c();
	trace_len_origin=c();
	trace_path=c(); #tax
	trace_type=c(); #taxnames 
	trace_path_type=0; #
	
	i=node1_number;
	j=node2_number;
			
			require(igraph); #Ensure that the library is loaded in each thread on some systems...
			
			rac=0;
			inf=0;
			detour=0;
			egal=0;
			dead=0;
			error=0;
			total_done=0;
			brk=0;
			deadend_or_detour=0;
				iso_g3short_ij=shortest.paths(g2_without_k, V(g2_without_k)[g1names[i]], V(g2_without_k)[g1names[j]], algorithm = "dijkstra", mode="out");
				g1short_ij=shortest.paths(g1, V(g1)[g1names[i]],V(g1)[g1names[j]],algorithm = "dijkstra", mode="out")
				gp=good_path2(g2_without_k,g1names[i],g1names[j],g2_unique_number_primed);
				# Original len
				trace_len_origin=g1short_ij;
				if(!is.finite(iso_g3short_ij)&&is.finite(g1short_ij)) {
					#error=error+1; #We are missing some edges in g2
					dead=dead+1;
					trace_path_type=c_roadblock;				
					trace_len=Inf;	
				} 			
				if(!is.finite(g1short_ij)) 
				{
					if(is.finite(iso_g3short_ij)&&gp){
						brk=brk+1;
						trace_path_type=c_breakthrough;
						#if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);
						#To trace this path we must call the shortest path method.
						trace_len=iso_g3short_ij;						
						if (trace) {
							 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j],mode="out")$res;	
							 #cat("\nTrace for Shortcut: ", g1names[i]," to ",g1names[j], "\n", sep="\t");							
							 if (length(paths)>0) {
								# Note: for now, we only get the first paths
								for (p in paths) {
									if (any((p[2:(length(p)-1)]) %in% g2_unique_number_primed)) {
										trace_path=c();
										trace_type=c();
										for (k in p) {    
											trace_path=c(trace_path,V(g2_without_k)[k]$name);
											trace_type=c(trace_type,V(g2_without_k)[k]$tax);
										 }
									}
								}
							 }						
						}
						
						
					}
					else {
						inf=inf+1;
						trace_path_type=c_impasse;
						trace_len=iso_g3short_ij;	
						#if (verbose) cat(g1names[i],g1names[j],"Disconnected",0,"\n",sep="\t",file=filename, append=TRUE);

					}
				} 
				else			
				{
					if(g1short_ij>iso_g3short_ij&&gp){
						rac=rac+1;
						trace_path_type=c_shortcut;
						trace_len=iso_g3short_ij;
						#if (verbose) cat(g1names[i],g1names[j],"Shortcut",iso_g3short,"\n",sep="\t",file=filename, append=TRUE);
						 if (trace) {
						 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j], mode="out")$res;	
							# cat("\nTrace for Shortcut: ", g1names[i]," to ",g1names[j], "\n", sep="\t");
						 	if (length(paths)>0) {
								# Note: for now, we only get the first paths
								for (p in paths) {
									if (any((p[2:(length(p)-1)]) %in% g2_unique_number_primed)) {
										trace_path=c();
										trace_type=c();
										for (k in p) {    
											trace_path=c(trace_path,V(g2_without_k)[k]$name);
											trace_type=c(trace_type,V(g2_without_k)[k]$tax);
										 }
									}
								}	
							 }	
					 	  }
					
					} 
					else if (iso_g3short_ij>g1short_ij&&gp) {
						detour=detour+1;
						trace_path_type=c_detour;
						#if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);				
						 if (trace) {
							#cat("Detour\n");	
							 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j], mode="out")$res;									
								trace_len=iso_g3short_ij;	
								if (length(paths)>0) {
										# Note: for now, we only get the first paths
										for (p in paths) {
											if (any((p[2:(length(p)-1)]) %in% g2_unique_number_primed)) {
												trace_path=c();
												trace_type=c();
												for (k in p) {    
													trace_path=c(trace_path,V(g2_without_k)[k]$name);
													trace_type=c(trace_type,V(g2_without_k)[k]$tax);
												 }
											}
										}	
								}	
							}	
					} else {					
						  
						# Test if the a and b can reach any k
						nopath_to_k=TRUE;								
						iso_g3short_i=shortest.paths(g2_without_k, g1names[i], algorithm="dijkstra")
						if (!directed) {
							iso_g3short_j=shortest.paths(g2_without_k, g1names[j], algorithm="dijkstra")
						} 
						if (length(g2_unique_names_primed)>0)
							for (kname in g2_unique_names_primed) {							
								
								if (!directed&&is.finite(iso_g3short_i[g1names[i],kname])&&is.finite(iso_g3short_j[g1names[j],kname])) {
									nopath_to_k=FALSE;	
								}
								#If directed, we only check i to k 
								if (directed&&is.finite(iso_g3short_i[g1names[i],kname])) {
									nopath_to_k=FALSE;	
								} 
							}
						if (no_new_node) {
							#Special case (Roadblock)
							if (g1short_ij==iso_g3short_ij) {
								dead=dead+1;
								trace_path_type=c_roadblock;
								trace_len=iso_g3short_ij;	
							}						
						} else if (nopath_to_k) {
							dead=dead+1;
							trace_path_type=c_roadblock;
							trace_len=Inf;	
							#if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=file, append=TRUE);		
						
						# Test if  there is a path from i to j with a k
						} else {						  
							paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j], mode="out")$res;					
							found=FALSE;						
							for (path in paths) {	
								for (k in path) {
									name=V(g2_without_k)[k]$name;
									if (!(name %in% g1names)&&name!=g1names[i]&&name!=g1names[j]) {
										found=TRUE;																	
									}
								}
							}
							if (found) {
								
								if (iso_g3short_ij>g1short_ij) {
									detour=detour+1;
									trace_path_type=c_detour;														
									#if (verbose) cat(g1names[i],g1names[j],"Detour",iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE); 
									 if (trace) {
										 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j], mode="out")$res;	
										 if (length(paths)>0) {									
											# We need to calculate the pathlen
											for (e in E(g2_without_k, path=paths[[1]])) {
												 w=ifelse(is.finite(E(g2_without_k)[e]$weight),E(g2_without_k)[e]$weight,1);
												trace_len=trace_len+w;
											}
											# Note: for now, we only get the first paths that might not be the good one here...											
											for (p in paths) {
												if (any((p[2:(length(p)-1)]) %in% g2_unique_number_primed))
													trace_path=c();
													trace_type=c();
													for (k in p) {    
														trace_path=c(trace_path,V(g2_without_k)[k]$name);
														trace_type=c(trace_type,V(g2_without_k)[k]$tax);
													 }									
												}
											}	
									 }
									
								} else {									
									egal=egal+1;
									trace_path_type=c_equal;
									#if (verbose) cat(g1names[i],g1names[j],"Equal", iso_g3short_ij,"\n",sep="\t", file=filename, append=TRUE);
									 if (trace) {
										 paths<-get.all.shortest.paths(g2_without_k,g1names[i],g1names[j], mode="out")$res;	
										trace_len=iso_g3short_ij;
										if (length(paths)>0) {											
											# Note: for now, we only get the first paths that might not be the good one here...
											# Longer buyt found the good path
											for (p in paths) {
												if (any((p[2:(length(p)-1)]) %in% g2_unique_number_primed)) {
													trace_path=c();
													trace_type=c();
													for (k in p) {    
														trace_path=c(trace_path,V(g2_without_k)[k]$name);
														trace_type=c(trace_type,V(g2_without_k)[k]$tax);
													}
												}	
											}										 
										}	
									 }
								}
							} else {
								#Possible dead-end
								deadend=TRUE;																
								g2_without_k_and_j=delete.vertices(g2_without_k,g1names[j]);
								g2_without_k_and_i=delete.vertices(g2_without_k,g1names[i]);
								iso_g3short_i=shortest.paths(g2_without_k_and_j, g1names[i],algorithm = "dijkstra")
								if (!directed) {
									iso_g3short_j=shortest.paths(g2_without_k_and_i, g1names[j],algorithm = "dijkstra")
								}
								list_of_knames<-c(); #for later use valid k for i and j
								order_of_list<-c(); #sort for faster								
								total_access_k=0;
								for (kname in g2_unique_names_primed) {
										# Note: we only add k if it's distance to i AND j is SMALLER than max_distance
									if (!directed) {
										if (is.finite(iso_g3short_i[g1names[i],kname])&&is.finite(iso_g3short_j[g1names[j],kname])) {																				
											total_access_k=total_access_k+1;
											if (maxdistance==0||(iso_g3short_i[g1names[i],kname]<maxdistance&&iso_g3short_j[g1names[j],kname]<maxdistance)) {
												aprox_len= max(iso_g3short_i[g1names[i],kname],iso_g3short_j[g1names[j],kname]);							
												list_of_knames<-c(list_of_knames, kname);								
												order_of_list<-c(order_of_list, aprox_len);
											}
										}
									} else {
										if (is.finite(iso_g3short_i[g1names[i],kname])) {																				
											total_access_k=total_access_k+1;
											if (maxdistance==0||(iso_g3short_i[g1names[i],kname]<maxdistance)) {
												aprox_len= iso_g3short_i[g1names[i],kname];							
												list_of_knames<-c(list_of_knames, kname);								
												order_of_list<-c(order_of_list, aprox_len);
											}
										}
									}
								}
								# Order the list of knames by distance 
								if (length(order_of_list)>0) {
									list_of_knames<-list_of_knames[order(order_of_list)];
									#order_of_list<-order(order_of_list); -- faster
								}
								# prime list of node if maxnode is specified
								if (maxnode>0&&length(list_of_knames)>0) {
									list_of_knames<-split_sample(list_of_knames, floor(maxnode/2))[[1]]; 
								}
																
								# look if we can acces
								ddlen=0; #length of detour
								tp0 <- proc.time(); #maxtime to search
								tp1<-proc.time()-tp0;
								if (length(list_of_knames)>0)
									for (kname in list_of_knames) {
								
									#find the shortest path between a and k 
									#and b and k 
									
									if (deadend) {
										#Check for maxtime.
										if (maxtime!=3600) {
											tp1<-(proc.time()-tp0)[[3]];
											if (tp1>maxtime) break;
										}	 
										# This is the really demanding (ressoure) question
										# Should we flag it and do it latter in a parallel ?
										# distinct thread?
									
										paths1<-get.all.shortest.paths(g2_without_k_and_j,g1names[i],kname, mode="out")$res;
										if (!directed) paths2<-get.all.shortest.paths(g2_without_k_and_i,g1names[j],kname, mode="out")$res;
										if (directed) paths2<-get.all.shortest.paths(g2_without_k_and_i,kname,g1names[j], mode="out")$res;
										#Ensure no intersection of paths
										if (length(paths1)>0&&length(paths2)>0) {
											#We look if the 2 shortest path do not intersect (have a common vertex here)
											#to avoid path like:
											#
											#                 /(j)
											#       (i)--a---b-----(k) (a and b in this case will (break) the path to k
											#
											#
											intersect=TRUE;
											for (p in paths1) { 										
												if (intersect) {												
													to_remove=c(V(g2_without_k_and_j)[kname]);
													trace_len=0;
													 for (e in E(g2_without_k_and_j, path=p)) {
															 w=ifelse(!is.null(E(g2_without_k_and_j)[e]$weight),E(g2_without_k_and_j)[e]$weight,1);
															 trace_len=trace_len+w;
													 }
													p = p[! p %in% to_remove]
													
													for (q in paths2) {
														to_remove=c(V(g2_without_k_and_i)[kname]);														
														q2=q;
														q = q[! q %in% to_remove]
														#Note that we only test from k to p 
														if (intersect) {							
														if(any(V(g2_without_k_and_j)[p] %in% V(g2_without_k_and_i)[q])) {
															vertex=V(g2_without_k_and_j)[p]$name %in% V(g2_without_k_and_i)[q]$name;
															p2= p[vertex]
															p2=V(g2_without_k_and_j)[p2]$name;														
															g2_without_k_and_p2=delete.vertices(g2_without_k,p2)
															if (!directed) paths3<-get.all.shortest.paths(g2_without_k_and_p2,g1names[j],kname, mode="out")$res;
															if (directed) paths3<-get.all.shortest.paths(g2_without_k_and_p2,kname,g1names[j], mode="out")$res;
															if (length(paths3)>0) {
																
																for (tp in p) {
																	trace_path=c(trace_path,V(g2_without_k_and_j)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_j)[tp]$tax);
																}																
																for (e in E(g2_without_k_and_p2, path=paths3[[1]])) {
																			 w=ifelse(!is.null(E(g2_without_k_and_p2)[e]$weight),E(g2_without_k_and_p2)[e]$weight,1);
																			 trace_len=trace_len+w;
																}
																for (tp in rev(paths3[[1]])) {		
																	trace_path=c(trace_path,V(g2_without_k_and_p2)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_p2)[tp]$tax);
																}
																intersect=FALSE;
															}													
														} else {														
															# Trace using p +k and q
																for (tp in p) {
																	trace_path=c(trace_path,V(g2_without_k_and_j)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_j)[tp]$tax);
																}							
																for (e in E(g2_without_k_and_i, path=q2)) {
																	 w=ifelse(!is.null(E(g2_without_k_and_i)[e]$weight),E(g2_without_k_and_i)[e]$weight,1);															
																	trace_len=trace_len+w;																			 
																}
																for (tp in rev(q2)) {															
																	trace_path=c(trace_path,V(g2_without_k_and_i)[tp]$name);
																	trace_type=c(trace_type,V(g2_without_k_and_i)[tp]$tax);
																}

															#trace_len=iso_g3short_i[g1names[i],kname]+iso_g3short_j[g1names[j],kname];
															intersect=FALSE;
														}
														}
													}	
												}	
											}
											if (!intersect) {																							
												deadend=FALSE;			
											}							
										}									
									} #if dead end
								} #end for k
								
								if (!deadend) {
									trace_path_type=c_detour;
									detour=detour+1;
									#if (verbose) cat(g1names[i],g1names[j],"Detour",ddlen,"\n",sep="\t", file=filename, append=TRUE);
								} else {
									#its a real dead-end if we have evaluated all possibilities...
									if ((maxdistance==0||total_access_k==length(list_of_knames))&&tp1<maxtime) {
										#trace_len=0;
										trace_path_type=c_roadblock;
										dead=dead+1;
										#if (verbose) cat(g1names[i],g1names[j],"Dead",0,"\n",sep="\t", file=filename, append=TRUE);
									} else {
										dead=dead+1;
										deadend_or_detour=deadend_or_detour+1;
										trace_path_type=c_roadblock;
									}	
								}
							}	
						}
					}
				}	
			total_done=total_done+1;
 temps1<-proc.time()-t0
 
sdis=inf;
segal=egal;
sdetour=detour;
srac=rac;
sdead=dead;
sdd=deadend_or_detour;
sbrk=brk;

utime=temps1[1]
stime=temps1[2]
rtime=temps1[3]

stotal<-sdis+srac+segal+sdetour+sdead+sdd+sbrk;
r=data.frame(sdis,srac,segal,sdetour,sdead,sdd, stotal, utime, stime, rtime)
# if (verbose) {
	# cat('Disconnected nodes      :', sdis,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Shortcuts               :', srac,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Equals                  :', segal,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Detours                 :', sdetour,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Dead ends               :', sdead,"\n", sep="\t", file=filename, append=TRUE);
	# str=paste('Dead ends or detour (maxdistance>',maxdistance,'):');
	# cat(str, sdd,"\n", sep="\t", file=filename, append=TRUE);
	# cat('Total                    :', stotal, "\n", sep="\t", file=filename, append=TRUE);
	# cat('Real Time                :', rtime, "\n", sep="\t", file=filename, append=TRUE);
# }
#ddname=paste('Dead ends or detour');
#colnames(r)=c('disconnected nodes','shortcuts','equals','detours','dead ends',ddname,'total','user time','system time','real time')
if (trace) {
	cat(g1names[i],"->",g1names[j],"\n");
	cat("original:",trace_len_origin,"\n");
	cat("augmented:",trace_len,"\n");
	print(trace_path);
	print(trace_type);
}
# Type of paths
if (trace_path_type==c_equal) ptype="Equal";
if (trace_path_type==c_shortcut) ptype="Shortcut";
if (trace_path_type==c_detour) ptype="Detour";
if (trace_path_type==c_roadblock) ptype="Roadblock";
if (trace_path_type==c_impasse) ptype="Impasse";
if (trace_path_type==c_breakthrough) ptype="Breakthrough";
#if (trace_path_type==c_dead_end_or_detour) ptype="Dead end or Detour";
if (trace_path_type==0) ptype="Undefined";

# Return a list
result<-list("from"=g1names[i], "to"=g1names[j], "path_type"=ptype,"path_type0"=trace_path_type,"original_path_length"=as.numeric(trace_len_origin), "augmented_path_length"=as.numeric(trace_len), "path"=trace_path, "path_visited_taxa"=trace_type);
return(result)
}

# New version Feb 2016
# This version record the different pattern of path
# e.g. Euk->Pla->Pla->Vir->Euk
# New version Jan 2019, added the option for incremental graph analysis
BRIDES<-function(X,Y=NULL,src="default", dest="default",A='default',random=0, maxdistance=100, maxtime=100,maxnode=100,maxcores=1, outfile="",size=1000, first=0, last=0, sample_paths=c(), weighted=FALSE, runmode="default", verbose=TRUE, min_additional=1, max_additional=0, max_iters=100, mutationChance=0.1, wt =c(-2,1,2,-1,-1,-2)) 
{
    #####################################################
	## Options
	if (!runmode %in% c("incremental","default","exhaustive")) {
	   warning("The runmode parameter must be either: 'default','incremental' or 'exhaustive'.\n");
	   return (NULL);
	}
	
	#####################################################
	## Wrapper
	g1=X
	g2=Y
	taxnames=A;
	npath=random;
	node1=src
	node2=dest
	randomize=F
	
	if (!is.null(g1$g1)) {
	 g2=g1$g2
	 g1=g1$g1
	}
	if (is.null(g2)&&taxnames=="default") {
		print(is.null(g2))
		print(taxnames=="default")
		warning("No augmented network set. Either supply two networks or supply one network with some attributes.\n");
	   return (NULL);
	}
	if (is.null(g2)) {
		g2=g1;
	}
	if ((is.directed(g1)&&!is.directed(g2))||(is.directed(g2)&&!is.directed(g1))) {
	  warning("Both networks must be either directed or undirected.\n");
	  return (NULL);
	}	
	
	######################################################
	## LOCAL VARIABLES
	directed=is.directed(g1);
	g1names<-V(g1)$name;    #list of vertex taxnames in g1
	g2names<-V(g2)$name;    #list of vertex taxnames in g2
	node1_number=0;		    #A single vertex from
	node2_number=0;		    #A single vertex to
	no_new_node=FALSE;       #Flag, if true, we only report the changed network topology changes.						
	
	#Test if all names in g1 are also in g2
	if (all(g1names %in% g2names)!=TRUE) {
		warning("! Warning ! Not all nodes in network X are in Y\n");	
		len_remove=length(g1names[!g1names %in% g2names]);
		g1=delete.vertices(g1,g1names[!g1names %in% g2names]);	
		g1names<-V(g1)$name;
	}
	#Replace edges with NA weight with 1
	E(g1)[is.na(weight)]$weight=1
	E(g2)[is.na(weight)]$weight=1
	# Be sure to remove weights
	if (weighted==FALSE) {
		E(g1)$weight=1
		E(g2)$weight=1
	}
	
	#####################################################
	# Incremental 2019 (call multiple time BRIDES       #
	# To try to find the 'BEST' additionnal nodes       #

	if (runmode!="default"&&is.null(g2))  {
		warning("In incremental or exhaustive mode, a source graph X (smaller) and destination graph Y(bigger) must be specified.\n");
		return (NULL);
	}
	
	if (runmode!="default") {
	  
		# 1. vertex in first graph
		# 2. vertex in second graph
		g2_unique_names<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;
		if (length(g2_unique_names)<2) {
			warning("In incremental mode, a source graph X (smaller) and destination graph Y(bigger) must be specified.\n");
			return (NULL);
		}
		#Weight for the optimization model score
		#      B R I D E  S  [Note, additionnal zero since vector contains other parameters]
		#wt =c(-2,1,2,-1,-1,-2,0,0) # scoring scheme
		
		if (length(wt)<8) {
			for (i in ((length(wt)+1):8)) {
				wt[i]=0;
			}
		}
		cat("Using the weighing scheme [ B R I D E S ]: ",wt[1:6],"\n")
		maxn=length(g2_unique_names);
		
		# Evaluation function 
		evalFunc=function(x){
		   if (sum(x)<1) return (0);
		   g2n=g2_unique_names[x>0]
		   g2n_names=c(g1names,g2n)
		   ing2=induced_subgraph(g2,V(g2)[g2n_names])
		   
			brides=BRIDES(g1,ing2,maxdistance=maxdistance, maxtime=maxtime,maxnode=maxnode,outfile="",size=size,verbose=FALSE, maxcores=maxcores)
			score=sum(wt*brides)
		  results=list()
		  results$score=score;
		  results$brides=brides;
		  results$nodes=paste(g2n, collapse=",");
  		return(results)
		}
		
		# Genetic algorithm 
		if (runmode=="incremental") {
			results=optimize_ga(size = maxn, iters = max_iters,mingenes=min_additional,maxgenes=max_additional,mutationChance=mutationChance,evalFunc = evalFunc)
		}
		if (runmode=="exhaustive") {
			results=allsolution(size = maxn, maxgenes=max_additional,evalFunc = evalFunc)
		}
		# Create the human-readable outputs
		results$nodes=g2_unique_names;
		df=data.frame(results$brides)
		rownames(df)=results$additional_nodes
		df$score=results$evaluations
		df$best=rep("",nrow(df))
		df$best[results$final_chromosome]="*"
		print(df)
		cat("* Best solutions found using criteria.\n")
		results$final_table=df;
		## Create the best networks
		results$networks=list()
		ni=1;
		for (idx in results$final_chromosome) {
		  chromosome=results$population[idx,]
		  g2n=g2_unique_names[chromosome>0]
		  g2n_names=c(g1names,g2n)
		  ing2=induced_subgraph(g2,V(g2)[g2n_names])
		  results$networks[[ni]]=ing2
		  ni=ni+1
		  plot(ing2,main=paste("Addition of ",results$additional_nodes[idx]))
		}
		
		# return the optimization results
		return(results);
	}

	
	
	################################################
	## Selection of the k node in the augmented graph
	if (taxnames=='default') {
		#we take all the node node in graph1			
		g2_unique_names<-V(g2)[!(V(g2)$name %in% V(g1)$name)]$name;		
	} else {
		g2_unique_names=V(g2)[V(g2)$tax==as.factor(taxnames)]$name;		
	}
	
	# Handle taxnames if not found
	if (is.null(V(g1)$tax)) V(g1)$tax=0;
	if (is.null(V(g2)$tax)) {
		V(g2)$tax=0;
		V(g2)[!(V(g2)$name %in% V(g1)$name)]$tax=1;
	}
	
################################################
	## If we have a node1, we only take this node1
	
	if (is.numeric(node1)) {
		node1_number=node1;
	} else if (node1!='default'){
		if (length(V(g1)[V(g1)$name==as.factor(node1)]$name)>0) {
			node1_number=match(node1, V(g1)$name)
		} else {
			warning(paste("Node with name :",node1," not found in X!\n"));
			return(c());
		}
	}
	###################################################
	## Look if node2 is specified
	if (is.numeric(node2)) {
		node2_number=node2;
	} else if (node2!='default'){
		if (length(V(g1)[V(g1)$name==as.factor(node2)]$name)>0) {
			node2_number=match(node2, V(g1)$name)
		} else {
			warning(paste("Node with name :",node2," not found in X!\n"));
			return(c());
		}	
	}
	if (node2_number!=0&&(node2_number==node1_number)) {
		#cat("Warning! Same number of nodes in network g1 and network g2\n");
		warning("Warning! Same number of nodes in network g1 and network g2\n");
		#return(c());
	}
	
	

	
	#################################################
	## Start of calculations
	##
	t0 <- proc.time()
	g2_degree_one=c()
	g2_unique_names_primed=c()  #Name of unique vertex in g2 without the degree one
	g2_unique_number_primed=c()  #Number of unique vertex
	# First prime not connected k
	if (verbose) cat("Priming unconnected nodes...\n");
	for (name in g2_unique_names) {
	
		if(degree(g2,name)==1) {
			g2_degree_one=c(g2_degree_one, name)
		} else {
			iso_g3short=shortest.paths(g2, v=name, algorithm="dijkstra");
			# Are we connected to any non k node 
			if(any(is.finite(iso_g3short[name,V(g1)$name]))) {
				g2_unique_names_primed=c(g2_unique_names_primed,name)	
			} else {
				g2_degree_one=c(g2_degree_one, name)
			}	
		}
	}	
	
	if (length(g2_degree_one)>0) {
		g2_without_k=delete.vertices(g2,which(V(g2)$name %in% g2_degree_one))
	} else {
		g2_without_k=g2;
	}
	g2_unique_number_primed=as.numeric(V(g2_without_k)[g2_unique_names_primed]);	
	#Test si tous les noms dans g1 sont aussi dans g3
	if (all(g1names %in% g2names)!=TRUE) {
		#cat("! Warning ! Not all name in g2 are in g1.\n");
		#if (verbose) cat("! Warning ! Not all name in g2 are in g1", file=filename, append=TRUE);	
		return(c());
	}

	if (length(g2_unique_names_primed)==0) {
			#cat("! Warning ! No new nodes accessibles in g2 from g1.\n");
			warnings("! Warning ! No new nodes accessibles in g2 from g1.\n");
			#, file=file, append=TRUE);
			# We call the new s
			no_new_node=TRUE;
			#return(c());
		}
	###################################################
	## Variables for dispatching to the different functions
	g1names=g1names[!g1names %in% g2_unique_names_primed]
	total_n=length(g1names);
	
	total_paths=(total_n*(total_n-1))/2;
	if (npath!=0&&npath<size) size=npath;
	if (directed)  total_paths=total_n*(total_n-1);
	if (npath!=0&&npath<total_paths) total_paths=npath;
	total_group=as.integer(total_paths/size)+1;
	if (total_paths %% size==0) total_group=total_group-1; #Correct for limit case
	if (last==0||last>total_group) last=total_group;
	if (first==0||first<1) first=1;
	
	
	####################################################
	## MAIN COMPUTING FUNCTION CALLS
	if (node2_number==0||node1_number==0||npath!=0) {
	##################################################
	## Implemented in a different function 
	## Added verbose for incremental (2019)
	if (verbose) {
		cat("==========================================================\n");
		cat("Total",total_paths,"pathways divided into", total_group, "groups.\n"); 
		cat("==========================================================\n");
		cat("Run parameters:\n");
		if (directed) {
		cat("Networks                    : directed\n");
		} else {
		cat("Networks                    : undirected\n");
		}
		cat("Nodes in networkX           :", total_n,"\n");
		cat("Nodes in networkY           :", length(V(g2)),"\n");
		cat("Total added nodes (K)       :", length(g2_unique_names),"\n");
		cat("Attributes for added nodes  :",taxnames,"\n");
		cat("Total paths                 :",total_paths,"\n");
		#cat("Randomize paths        :",randomize,"\n");
		cat("Weighted edges              :",weighted,"\n");
		cat("Group size                  :",size,"\n");
		cat("Start group                 :",first,"\n");
		cat("End group                   :",last,"\n");
		cat("Maxdistance                 :",maxdistance,"\n");
		if (maxtime!=100) cat("Maxtime                     :",maxtime,"\n");
		cat("Maxnode                     :",maxnode,"\n");
		cat("Maxcores                    :",maxcores,"\n");  
		cat("==========================================================\n");
	}
	if (outfile!="") write("# src\tdest\tdist_x\tdist_y\tBRIDES\tpath\ttaxa\n",outfile)
		if (npath==0) {
			npath=(total_n*(total_n-1))/2;
			if (directed) npath=npath*2;
		}
		
		cl <- makeCluster(multicore(maxcores))
			registerDoParallel(cl=cl);
			
		total_s=array(0,8);
		if (verbose) cat("#  ","B","R","I","D","E","S","(utime","stime)","\n", sep="\t")
		
		for (p in first:last) {
			pathways=c();
			if (!directed) {
				sample_paths=create_sample_path(total_n, p, size);
			} else {
				sample_paths=create_sample_path(total_n, p, size,TRUE);
			}
			sample_paths=sample_paths[sample_paths != 0] #Trim the sample_path of zero (0) indice
			npath=length(sample_paths)/2;
			total=0
			t0=proc.time();
			s<-foreach(hi=1:npath, .combine=function(i,j, .export=c("pathways")) { return(i+j)}) %dopar% {
			#For debug
			#s=c()
			#for (hi in 1:npath) {
				if(!exists("pathBRIDES", mode="function")) source("BRIDES_v2019.R")
				i=sample_paths[(hi-1)*2+1];
				j=sample_paths[(hi-1)*2+2];				
				b=pathBRIDES(g1,g2_without_k,g1names,g2_unique_names_primed,g2_unique_number_primed,i,j,no_new_node,maxdistance,maxtime,maxnode,proc.time());
				a=array(0,6)
				a[b$path_type0]=1;
				pathways=c(pathways,b);
				pa=paste(b$path,collapse=",");
				pav=paste(b$path_visited_taxa,collapse=",");
				cat(g1names[i],g1names[j],b$original_path_length,b$augmented_path_length,b$path_type,pa,pav,"\n",sep="\t",file=outfile, append=T);
				a
				#s=c(s,a)
			}
			stopCluster(cl)
			ttime=(proc.time()-t0)
			s=c(s,as.numeric(ttime[1]),as.numeric(ttime[3]))			
			total_s = total_s +s
			names(s)<-c("B","R","I","D","E","S","(utime","stime)")
			if (verbose) cat("",s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],"\n",sep="\t");
		} # End first to last group 
		if (verbose)  cat("====================== RESULTS ===========================\n");
		names(total_s)<-c("B","R","I","D","E","S","(utime","stime)")
		return(total_s)
	# CASE 2. We have one sample 
	} else {
		pathBRIDES(g1,g2_without_k,g1names,g2_unique_names_primed,g2_unique_number_primed,node1_number,node2_number,no_new_node,maxdistance,maxtime,maxnode,proc.time());
	}
	
}

###################################################
	#Any good path
	good_path2<-function(g,node1,node2, additional_node) {
		path<-get.all.shortest.paths(g,node1,node2, mode="out")$res;	
		for (p in path) {
			if (length(p)<3) return (FALSE);
			if (length(unique(p))!=length(p)) return (FALSE);			
			if (any((p[2:(length(p)-1)]) %in% additional_node)) return (TRUE);  
		}
		return (FALSE);
	}	
	
	
