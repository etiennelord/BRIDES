source("BRIDES_2022.R")
options(error = function() traceback())
# This generate some random networks and then, execute the C++ version of BRIDES
# to evaluate different paths using some heuristics.
# See: https://doi.org/10.1371/journal.pone.0161474
# Repeat 1000
nreplicate=1000;
time_h1=c();
time_h2=c();
time_h3=c();
time_h4=c();
time_h5=c();

for (nt in c("erdos","barabasi","watts")) {
for (nr in c(5,25,50,100)) {
		for (i in 1:nreplicate) {
			r=random_network(100,nr,type=nt);
			set.seed(i);
			j=1;
			g1names=paste(nt,"_",i,"_",nr,"_g1.txt",sep="");
			g2names=paste(nt,"_",i,"_",nr,"_g2.txt",sep="");
			export_network(r$g1,g1names,T);
			export_network(r$g2,g2names,T);
			#Run some heuristic found in the CPP version (1-BRIDES, 2-BRIDES_Y, 3-BRIDES_YC, 4-BRIDES_EC,5-DFS)
			for (h in 1:5) {
				tp0 <- proc.time();
				grnames=paste0(nt,"_",i,"_",nr,"_h",h,"_result.txt");
				#Note, we do 100 random path from each network
				#Note: this use the cpp version 
				gcommand=paste("./brides -X=",g1names," -Y=",g2names," -heuristic=",h," -seed=",i," -random=100 -output=",grnames,sep="");
				t1 <- try(system(gcommand, intern = TRUE))				
				temps1<-proc.time()-tp0
				cat(h,temps1[3],"\n")
				# save as nt,nr,rep,system time,user time
				if (h==1) time_h1=c(time_h1,c(nt,nr,i,temps1[2],temps1[3]));
				if (h==2) time_h2=c(time_h2,c(nt,nr,i,temps1[2],temps1[3]));
				if (h==3) time_h3=c(time_h3,c(nt,nr,i,temps1[2],temps1[3]));
				if (h==4) time_h4=c(time_h4,c(nt,nr,i,temps1[2],temps1[3]));
				if (h==5) time_h5=c(time_h5,c(nt,nr,i,temps1[2],temps1[3]));
				
				write(time_h1, ncolumns=5,file="ttime1.txt",sep=",");
				write(time_h2, ncolumns=5,file="ttime2.txt",sep=",");
				write(time_h3, ncolumns=5,file="ttime3.txt",sep=",");
				write(time_h4, ncolumns=5,file="ttime4.txt",sep=",");
				write(time_h5, ncolumns=5,file="ttime5.txt",sep=",");
			}
		}	
	}
}  

