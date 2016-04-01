/* 
 * File:   main.cpp
 * Author: Etienne Lord
 * Created on 4 February 2016, 07:28
 */

#include <iostream>       // std::cout
#include <queue>          // std::priority_queue
#include <map>            // std::map
#include <set>            // 
#include <cstring>        // std::string
#include <fstream>        // IO
#include <sstream>
#include <climits>
#include <algorithm>      // std::find
 #include <utility>
#include <stdlib.h>       // srand, rand 
#include <ctime>
#include <cmath>          // floor
#include <time.h>         // time 
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
   #define omp_get_max_threads() 0
   #define omp_get_num_procs() 0
#endif
using namespace std;

#define Inf 1e14
// MAX POSSIBLE NODE for SPEED 10*1024*1024
#define  maxnodes 10485760
// MAX POSSIBLE SPP between i and j investigated
//#define max_individual_path 100
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
		
#define c_equal         4
#define	c_shortcut      5
#define	c_detour        3
#define	c_roadblock     1
#define	c_impasse       2
#define	c_breakthrough  0


const char *description = "===============================================================================\n"
                          "| Program : BRIDES 2016                                                       |\n"
                          "| Authors : Etienne Lord,Vladimir Makarenkov (Universite du Quebec a Montreal)|\n"
                          "|           and Francois-Joseph Lapointe     (Universite de Montreal)         |\n"
                          "| This program computes the similarities in path between two related networks.|\n"
                          "===============================================================================\n";

const char *startMessage = "BRIDES V.1.0 - (2015) by Etienne Lord and Francois-Joseph Lapointe\n"; 

/////////////////////////////////////////////////////////////////////////////////////
//--Structure definition

template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
			out << v[i];
        if (i != last) 
            out << ", ";
    }
    out << "]";
    return out;
}


struct edge
{
    int from, to;
    float dist;
        
    edge() {dist=1.0f;}
    
    bool operator<(const edge& rhs) const
    {
        return dist < rhs.dist;
    }
};
   
struct Parameters {
	char graph1[1024];     // filename g1, X
	char graph2[1024];     // filename g2,Y
	char attributes[1024]; // filename attributes
	char graphinfo[1024];  // Node appartenance
	char outputfile[1024]; //outputfile	
	char debugfile[1024]; //outputfile	
	vector<std::string> K;//attributes to use for K nodes
	vector<std::string> nonK;//attributes to use for non K nodes 
	
        bool found_g1;//-found g1 file?
        bool found_g2;//-found g2 file?
        bool found_attributes; //--Annotations
        bool info;
		bool use_dist; //by default false, (dist set to 1)
		bool inv_dist; // use the inverse of the distance
		
        bool directed;
		int seed;  //--seed
		float random; //--number of random path
		int size;
        int first; //-- start group
        int last;   //--End verte
        float maxdistance; // default search length (default=20);		
		int maxtime;
        float maxnode;
		bool verbose; 
	
	int max_individual_path; 	
	int weight_method; //All weight as equal, inverse, default
		
	
	bool use_multiple_color; // do we required multiple color?
      
};

/////////////////////////////////////////////////////////////////////////////////////
//--Function definition
void constructor();
void destructor();

/////////////////////////////////////////////////////////////////////////////////////
//--Graph
void dijkstra (int S, vector< vector<int> > adj,vector< vector<float> > adj_dist, vector<float> &dist, vector< map<int,int> > &prev);
void dijkstra (int S, vector< vector<int> > adj,vector< vector<float> > adj_dist, vector<float> &dist, vector< map<int,int> > &prev, int ignored_nodes);
void dijkstra (int S, int T,vector< vector<int> > adj,vector< vector<float> > adj_dist, vector<float> &dist, vector< map<int,int> > &prev, map<int,int> ignored_nodes);

/////////////////////////////////////////////////////////////////////////////////////
//--Other functions
int readParameters(Parameters *param, char **argv, int nargc);
int ExtraireDonnees(const char * chaine, char *champs, char * contenu);
void help();
void output_header(char **argv);
void build_adjacency_list();
void load_info(char* filename);
int get_node(std::string name);
void load_attributes(char *filename);
std::map<int,int> load_edge(char *filename, std::vector<edge> &edge_list);
int add_node(std::string name);
void cout_map(vector< map<int, int> > p);
int randomInt (int min, int max);

/////////////////////////////////////////////////////////////////////////////////////
//--Variables
struct Parameters param; //--Input parameters

std::map<std::string, int> node_name_lookup; //--map to easy access to node name
std::vector<std::string> node_name;			 //--node_name, numbering start at (0..N)
std::vector<edge> edge_list_g1;              //--Edge list for graph1 (X): Note empty node is -1
std::vector<edge> edge_list_g2;				 //--Edge list for graph2 (Y)
std::map<int,int> node_id_g1;                //--Map of the nodes in g1 (faster)
std::map<int,int> node_id_g2;				 //--Map of the nodes in g2 (faster)
std::map<int, std::string> attributes;
std::map<std::string, int> uniques_attributes; //--Atrtibutes and their count


vector< vector<int> > undirected_adjlist_g1;
vector< vector<int> > undirected_adjlist_g2;
vector< std::map<int,float> > undirected_adjlist_g1_dist;
vector< std::map<int,float> > undirected_adjlist_g2_dist; 
vector< vector<float> > statistics; //statistics for each nodes
vector<int> random_paths; // unique random paths 

bool* node_graph1; //Node in node_name in graph1  
bool* node_graph2; //Node in node_name in graph2 // 
int current_node=0;
int total_n=0; //--Internal Counter 
int total_n_g1=0;
int total_n_g2=0;
bool* Nsb;   //--Vector from 0..n of node_appartenance (if true, found in g2 only) 
bool* Nsu;   //--Vector from 0..n of node_usage (default true, we use all nodes)
//--Internal informations
int total_group=0;
int size_of_groups=0;
float total_paths=0;
int total_k_node=0;
int total_nonk_node=0; //--shoudl be equal to g1 normally
	
ofstream FileOutput; //verbose path to output

void add_statistics(int src_id,int dest_id,int path_type, vector <int> path) {
	// Array statistics:
	// 0 1 2 3 4 5  6            7             8 
	// B R I D E S  inside_path  Min_len_to_k  Max_len_to_k
	statistics[src_id][path_type]++;
	statistics[dest_id][path_type]++;
	 for (int i=0; i<path.size();i++) {
	   if (path[i]!=src_id&&path[i]!=dest_id) statistics[path[i]][6]++;
	 }
	
}

void debug_info() {
	for (int i=0; i<total_n;i++) {
		cout<<node_name[i]<<" "<<Nsb[i]<<" "<<attributes[i]<<" "<<node_id_g1.count(i)<<" "<<node_id_g2.count(i)<<endl;
	}	
}

void dijkstra (int S, int T,vector< vector<int> > adj, vector< std::map<int,float> > adj_dist,vector<float> &dist, vector< map<int,int> > &prev, map<int,int> ignored_nodes) {			
			int len=adj.size();
			prev.clear();			
			for (int i=0; i<len; i++) {
				map<int,int> tmp;
				prev.push_back(tmp);
			}
			std::fill(dist.begin(), dist.end(), 1e14);
		    dist[S]=0.0f;
			multimap<float,int> pq;
			multimap<float,int>::iterator it;			
			pq.insert(multimap<float,int>::value_type(0.0,S));
			while (!pq.empty()) {
					it = pq.begin();
					int u=(*it).second;					
					float cost_to_u=(*it).first;
					pq.erase(it);
					if (u==T) break;
					if (cost_to_u <= dist[u]) {
						for ( int j=0;j<adj[u].size();j++ ) { 
							 int v=adj[u][j];
							 float Duv = adj_dist[u][v]; // Change here if edge is weighted
							if (ignored_nodes.count(v)==0) {
								float new_cost = cost_to_u + Duv;
								if (new_cost<dist[v]) {
									dist[v]=new_cost;
									prev[v].clear();
									prev[v].insert(std::pair<int,int>(u,u));
									pq.insert(multimap<float,int>::value_type(dist[v],v));
								} else if (new_cost==dist[v]) {
									prev[v].insert(std::pair<int,int>(u,u));
								}
							}
						}					
					}
					
			}					
	 }

void dijkstra (int S, vector< vector<int> > adj, vector< std::map<int,float> > adj_dist,vector<float> &dist, vector< map<int,int> > &prev, int ignored_nodes) {
			int len=adj.size();
			prev.clear();			
			for (int i=0; i<len; i++) {
				map<int,int> tmp;
				prev.push_back(tmp);
			}
			std::fill(dist.begin(), dist.end(), 1e14);
		    dist[S]=0.0f;
			multimap<float,int> pq;
			multimap<float,int>::iterator it;			
			pq.insert(multimap<float,int>::value_type(0.0,S));
			while (!pq.empty()) {
					it = pq.begin();								
					int u=(*it).second;					
					float cost_to_u=(*it).first;			
					pq.erase(it);		
					if (cost_to_u <= dist[u]) {
						for ( int j=0;j<adj[u].size();j++ ) {
							int v=adj[u][j];
							 float Duv = adj_dist[u][v]; // Change here if edge is weighted
							if (ignored_nodes!=v) {
								float new_cost = cost_to_u + Duv;
								if (new_cost<dist[v]) {
									dist[v]=new_cost;									
									prev[v].clear();									
									prev[v].insert(std::pair<int,int>(u,u));
									pq.insert(multimap<float,int>::value_type(dist[v],v));
								} else if (new_cost==dist[v]) {
									prev[v].insert(std::pair<int,int>(u,u));
								}
							}
						}					
					}
					
			}		
	 }	 
	 
 void dijkstra (int S, vector< vector<int> > adj,vector< std::map<int,float> > adj_dist, vector<float> &dist, vector< map<int,int> > &prev) {
			int len=adj.size();
			prev.clear();			
			for (int i=0; i<len; i++) {
				map<int,int> tmp;
				prev.push_back(tmp);
			}
			std::fill(dist.begin(), dist.end(), 1e14);
		    dist[S]=0.0f;
			multimap<float,int> pq;
			multimap<float,int>::iterator it;
			//for (int i=0; i<adj.size();i++) if (i!=S) pq.insert(multimap<float,int>::value_type(1e14,i));
			pq.insert(multimap<float,int>::value_type(0.0,S));
			while (!pq.empty()) {
					it = pq.begin();						
					int u=(*it).second;					
					float cost_to_u=(*it).first;
					pq.erase(it);			
					if (cost_to_u <= dist[u]) {
						for ( int j=0;j<adj[u].size();j++ ) {
							int v=adj[u][j];
							 float Duv = adj_dist[u][v]; // Change here if edge is weighted
							float new_cost = cost_to_u + Duv;
							if (new_cost<dist[v]) {
								dist[v]=new_cost;
								prev[v].clear();
								prev[v].insert(std::pair<int,int>(u,u));
								pq.insert(multimap<float,int>::value_type(dist[v],v));
							} else if (new_cost==dist[v]) {
								prev[v].insert(std::pair<int,int>(u,u));
							}
						}					
					}
					
			}			
	 }


 // use DFS
     void dfs_path(int s, int t, vector< map<int,int> > prev, vector< vector<int> > &paths, vector<int> &path, vector<bool> &onPath ) {		 
        if (paths.size()>param.max_individual_path) return;
		// found path from s to t
        if (s == t) {    			
			vector<int> tmp=path;
			paths.push_back(tmp);			
			return;
        } 
		  for(map<int,int>::iterator it = prev[t].begin(); it != prev[t].end(); ++it) {
			 int w=(*it).first; 	
			vector<int> tmp=path;			 
			if (!onPath[w]) {
				path.push_back(w);
				onPath[w] = true;
				  dfs_path(s, w, prev,paths, path, onPath);
				  path.pop_back();
				onPath[w] = false; 
			}
		  }		           	
    }


//--This will output only the first path
vector< vector<int> > get_path(int S, int T, vector< map<int,int> > prev) {
	vector< vector<int> > paths;
	vector<bool> onPath(prev.size(), false);
	onPath[T]=true;
	vector<int> path;	
	path.push_back(T);
	dfs_path(S,T, prev,paths, path,onPath);		
	
	for (int i=0; i<paths.size();i++ ) {
		std::reverse(paths[i].begin(), paths[i].end());
	}
	if (paths.size()==1&&paths[0].size()==0) paths.clear();
	return paths;
}
	 
	 void cout_map(vector< map<int, int> > p) {
		 cout<<"[ ";
		 for (int i=0; i<p.size();i++) {
			 vector<int> v;
			 for(map<int,int>::iterator it = p[i].begin(); it != p[i].end(); ++it) {
			  v.push_back(it->first);
			}
			cout<<v;
			if (i!=p.size()-1) cout<<", ";
		 }
		 cout<<" ]"<<endl;
	 }
	 
	 void cout_map( map<int, int>  p) {
		 vector<int> v;
			 for(map<int,int>::iterator it = p.begin(); it != p.end(); ++it) {
			  v.push_back(it->first);
			}
			cout<<v;			
	 }

float dist_path(vector<int> path,vector< std::map<int,float> > adj_dist) {
		float len=0.0f;
		
		if (path.size()<1) return Inf;
		// any duplicate and path contains one k nodes
		int src=path[0];
		for (int i=1; i<path.size();i++) {
			int dest=path[i];
			len+=(int)adj_dist[src][dest];
			src=dest;
		}
	return len;
}	 
	 
bool good_path(vector<int> path, int S, int T) {
		bool found=false;
		bool found_S=false;
		bool found_T=false;
		map<int,int> ht;
		
		if (path.size()<3) return false;
		// any duplicate and path contains one k nodes
		for (int i=0; i<path.size();i++) {
			if (Nsb[path[i]]) found=true;
			if (ht[path[i]]!=0) return false;
			ht[path[i]]=1;
			if (path[i]==T) found_T=true;
			if (path[i]==S) found_S=true;
		}
		
	return found&&found_S&&found_T;
}	

vector<int> good_path2(vector< vector<int> > paths, int S, int T) {
	vector<int> tmp;
	for (int i=0; i<paths.size();i++) {
		if (good_path(paths[i],S,T)) {
			for (int j=0; j<paths[i].size();j++) tmp.push_back(paths[i][j]);
			return tmp;
		}
	}
	return tmp;
}

void cout_paths(vector<int> path) {
		for (int j=0; j<path.size();j++) {
			cout<<node_name[path[j]];
			if (j<path.size()-1) cout<<",";
		}
		cout<<endl;
}

void cout_paths(vector< vector<int> > paths) {
	for (int i=0; i<paths.size();i++) {
		for (int j=0; j<paths[i].size();j++) cout<<node_name[paths[i][j]]<<" ";
		cout<<endl;		
	}
}

void fout_paths(vector<int> path) {
		for (int j=0; j<path.size();j++) {
			FileOutput<<node_name[path[j]];
			if (j<path.size()-1) FileOutput<<",";
		}
		FileOutput<<"\t";
		//TO DO, attribute here
		for (int j=0; j<path.size();j++) {
			FileOutput<<attributes[path[j]];
			if (j<path.size()-1) FileOutput<<",";
		}
		FileOutput<<endl;
}

vector<int> join(vector<int> v1, vector<int> v2) {
	vector<int> tmp;
	reverse(v2.begin(),v2.end());
		
	for (int i=0; i<v1.size();i++) tmp.push_back(v1[i]);
	for (int i=1; i<v2.size();i++) tmp.push_back(v2[i]);
	return tmp;
}

map<int,int> intersect(vector<int> v1, vector<int> v2,int k,int i2) {
	map<int,int> tmp;
	vector<int> to_erase;
	for (int i=1;i<v1.size();i++) {
		tmp.insert ( std::pair<int,int>(v1[i],0));
		tmp[v1[i]]++;
	}
	for (int i=1;i<v2.size();i++) {
		if (tmp[v2[i]]>0) {
			tmp[v2[i]]++;
			to_erase.push_back(v2[i]);
		} else {
			tmp.insert ( std::pair<int,int>(v2[i],0));
		}
	}
	tmp.erase(k);
	for (int i=0;i<to_erase.size();i++) tmp.erase(to_erase[i]);
	tmp.insert ( std::pair<int,int>(i2,0)); ///--becase we want to find a path from k to j 
	return tmp;
}

// Function createSamplePath
// This function return an array of paired node number 
// of size (size). The array returned is the one 
// specified by the offset
// 
// Arguments:
// stotal        = number of nodes in g1
// soffset (1:...) = starting position in pathways
//                   note: the total number of offset is ((total_n*(total_n-1)/2)/size)+1 --undirected
//                  total number of offset is ((total_n*(total_n-1))/size)+1
// size           = size of sampling (returned array)
// 

vector<int> createSamplePath(int soffset, float size, bool directed, bool random){
int i, j, k, t;

k=0;
t=0;
float total_path=(total_nonk_node*(total_nonk_node-1))/2;
if (directed) total_path*=2;
//--Handle the random path

vector<int> path;
if (random) {
	//Calculate total node in g2 found in g1 
	float total_g2_g1=0;
	for (int i=0; i<total_n_g1;i++) {
		if (node_id_g2.count(i)>0&&Nsu[i]) total_g2_g1++;
	}
	
	total_g2_g1=(total_g2_g1*(total_g2_g1-1))/2;
	if (directed) total_g2_g1*=2;
	
	//Now calculate the size of the random path 
	std::map<std::pair<int,int>, int> created;
	if (size>total_g2_g1) size=total_g2_g1;
	while (created.size()<size) {
		int s=randomInt(0,total_n_g1);
		int t=randomInt(0,total_n_g1);
		if (s>t&&!directed) {
			int w=s;
			s=t;t=w;
		}
		if (s!=t&&node_id_g2.count(s)>0&&node_id_g2.count(t)>0&&!Nsb[s]&&!Nsb[t]&&Nsu[s]&&Nsu[t]) created.insert(std::pair<std::pair<int,int>, int>(std::pair<int,int>(s,t),0));
	}
	for (std::map<std::pair<int,int>,int>::iterator it=created.begin(); it!=created.end(); ++it) {
		 path.push_back(((*it).first).first);
		  path.push_back(((*it).first).second);
		
	}
} else if (!directed) {
	 for(j = 0; j < total_n_g1; j++) {
		for(i = j+1; i < total_n_g1; i++){
			if (i>j&&node_id_g2.count(i)>0&&node_id_g2.count(j)>0&&!Nsb[i]&&!Nsb[j]&&Nsu[i]&&Nsu[j]) {			
				if ((int)((t/size))+1==soffset) {			
					path.push_back(i);
					path.push_back(j);
				}
				t++;
			}
			if (((int)(t/size)+1)>soffset) break;
		}
	}	
 } else {
	  for(j = 0; j < total_n_g1; j++) {
		for(i = 0; i < total_n_g1; i++){
			
			if (i!=j&&node_id_g2.count(i)>0&&node_id_g2.count(j)>0&&!Nsb[i]&&!Nsb[j]) {			
				if ((int)((t/size))+1==soffset) {			
					path.push_back(i);
					path.push_back(j);					
				}
				t++;
			}
			if (((int)(t/size)+1)>soffset) break;
		}
	}	 
 } 
 return path;
}
	 
vector<int> brides(int group) {	
	 bool verbose=param.verbose;
	 int total=0;	 //--Total path investigated
	 int B=0;
	 int R=0;
	 int I=0;
	 int D=0;
	 int E=0;
	 int S=0;
	 int total_time=0;
	int total_to_find=(total_nonk_node*(total_nonk_node-1))/2;
	 int last_p=0;
	 int iter_n=total_nonk_node;
	 vector<int> random_path;	 
	 //--Handle random paths or normal path
	 if (param.random==-1) {		 		 
		 random_path=createSamplePath(group,param.size,param.directed,false);
	 } else {
		 //calculate start node
		 int start_n=((group-1)*param.size)*2;
		 int last_n=start_n+(2*param.size);
		 if (last_n>random_paths.size()) last_n=random_paths.size();
		 for (int i=start_n; i<last_n;i++) {
			 random_path.push_back(random_paths[i]);
		 }
	 }	 
	iter_n=(random_path.size()/2);		 	 
	 #if defined(_OPENMP)
	 #pragma omp parallel for reduction(+:B,R,I,D,E,S,total,total_time) shared (random_path,statistics)
	 #endif
	 for (int ii=0;ii<iter_n;ii++) {
		 int i=random_path[2*ii];
		 vector < map<int,int> > prev_g1(total_n); //--Reserve total_n ;
		 vector<float> dist_g1(total_n);
		 dijkstra(i,undirected_adjlist_g1,undirected_adjlist_g1_dist, dist_g1, prev_g1);
		 vector < map<int,int> > prev_g2(total_n); //--Reserve total_n;
		 vector<float> dist_g2(total_n);
		 dijkstra(i,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2, prev_g2);
		 bool more=true; // we have more path to process?
		 // Min-Max path to k - Not good since we don't process all i (i>j for undirected)		 
		// if (param.verbose) {		
			// for (int l=0; l<total_n;l++) {
				// if (Nsb[l]) cout<<node_name[i]<<" "<<node_name[l]<<" "<<dist_g2[l]<<endl;
				// if (dist_g2[l]<Inf&&l!=i&&Nsb[l]&&dist_g2[l]<statistics[i][7]) statistics[i][7]=dist_g2[l];
				// if (dist_g2[l]<Inf&&l!=i&&Nsb[l]&&dist_g2[l]>statistics[i][8]) {
					// statistics[i][8]=dist_g2[l];
				// }
			 // }
		// }
		 std::clock_t start_time = std::clock();
		 while (more) {
			int  j=random_path[2*ii+1];  
			 int path_type=-1;
			 float real_dist_g2=0.0f;
 			
			 vector< vector<int> > path=get_path(i,j,prev_g2);
			 vector<int> gp=good_path2(path,i,j); 
			 real_dist_g2=dist_g2[j];
			 if (dist_g1[j]<Inf&&dist_g2[j]>=Inf) {
				path_type=c_roadblock;
			 } else			
			 if (dist_g1[j]>=Inf) {
			
				if (dist_g2[j]<Inf&&gp.size()>0) {
					path_type=c_breakthrough;
				} else {
					path_type=c_impasse;
				}				
			 } else {
				if (dist_g1[j]>dist_g2[j]&&gp.size()>0) {
					path_type=c_shortcut;
				} else if (dist_g1[j]<dist_g2[j]&&gp.size()>0) { 
					path_type=c_detour;
				} else if (dist_g1[j]==dist_g2[j]&&gp.size()>0) { 
					path_type=c_equal;
				} else {
					// Main big loop
					// test for k accessible nodes
					vector < map<int,int> > prev_g2_w_j(0); //--Reserve 10 ;
					vector < map<int,int> > prev_g2_w_i(0); //--Reserve 10 ;
					
					vector<float> dist_g2_w_j(total_n);
					vector<float> dist_g2_w_i(total_n);
					dijkstra(j,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_i, prev_g2_w_i,i);
					dijkstra(i,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_j, prev_g2_w_j,j);
					 std::vector<int> k_nodes;     //--K nodes (that we need to include in path)
					
					// TO DO here, test for distance and order the k by this distance
						multimap<float,int> k_nodes_dist;
						for (int l=0; l<total_n;l++) {
							if (dist_g2_w_j[l]<Inf&&dist_g2_w_i[l]<Inf&&Nsb[l]) k_nodes_dist.insert(multimap<float,int>::value_type(std::max(dist_g2_w_j[l],dist_g2_w_i[l]),l));
						}
				
												
						int knodes_size=0;
						float last_distance=0.0f;					
						for (std::multimap<float,int>::iterator it=k_nodes_dist.begin(); it!=k_nodes_dist.end(); ++it)
						{															
							if ((*it).first<param.maxdistance&&knodes_size<param.maxnode) {								
								k_nodes.push_back((*it).second);
								last_distance=(*it).first;
								knodes_size++;
							}
							if (last_distance>param.maxdistance||knodes_size>param.maxnode) break;
						}
						
					if (k_nodes.size()==0) {
						path_type=c_roadblock;	
					} else {
						 int ms = (std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000);
							// #We look if the 2 shortest path do not intersect (have a common vertex here)
							// #to havoid path like:
								// #
								// #                 /(j)
								// #       (i)--a---b-----(k) (a and b in this case will (break) the path to k
								// #
								// #
								bool nodeadend=false;																
								gp.clear();
								// 1. We try with only the path found for i to k and j to k
								//cout<<knodes_size;
								for (int l=0; l<k_nodes.size();l++) {
									int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
									if (elapsed>param.maxtime) {
										path_type=c_roadblock;
										break;
									}
									int k=k_nodes[l];
									
									//cout<<k<<" [k_nodes.size():"<<k_nodes.size()<<"]\n";
									vector< vector<int> > path_ik=get_path(i,k,prev_g2_w_j);
									//cout<<path_ik.size();
									vector< vector<int> > path_jk=get_path(j,k,prev_g2_w_i);																											
									// cout<<node_name[i]<<"-"<<node_name[k]<<"|\n";cout_paths(path_ik);
									 // cout<<node_name[j]<<"-"<<node_name[k]<<"|\n";cout_paths(path_jk);
									for (int pi=0; pi<path_ik.size();pi++) {
										for (int pj=0;pj<path_jk.size();pj++) {											 
											 vector<int> tmp=join(path_ik[pi],path_jk[pj]);
											 if (good_path(tmp,i,j)) {
												 if (gp.size()==0||gp.size()>tmp.size()) {
													 //gp.clear();
													 //for (int m=0; m<tmp.size();m++) gp.push_back(tmp[m]);
													 gp=tmp; 
													 //cout<<gp2<<endl;
												 }
												 nodeadend=true;												
											 } 
											 //if(nodeadend) break;
										}	
										//if(nodeadend) break;
									}									
									//--2. Try again with opt. but costly when we remove nodes duplicated in path i to k and j to k
									if (!nodeadend) {
										gp.clear();
										for (int pi=0; pi<path_ik.size();pi++) {
												int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
												if (elapsed>param.maxtime) {
													path_type=c_roadblock;
													break;
												}
												for (int pj=0;pj<path_jk.size();pj++) {											 
													 map<int,int> bad_nodes=intersect(path_ik[pi], path_jk[pj],k,i);
													 if (bad_nodes.size()>2) {														
														 vector < map<int,int> > prev_g2_w_i2(total_n); //--Reserve 10 ;	
														 vector<float> dist_g2_w_j2(total_n);														 
														 dijkstra(j,k,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_j2, prev_g2_w_i2,bad_nodes);														 
														 vector< vector<int> > path_jk2=get_path(j,k,prev_g2_w_i2);													
														  //cout<<node_name[j]<<"-"<<node_name[k]<<"|\n";cout_paths(path_jk2);
														 for (int pk=0; pk<path_jk2.size();pk++) {
															 vector<int> tmp2=join(path_ik[pi],path_jk2[pk]);													 
															 if (good_path(tmp2,i,j)) {
																if (gp.size()==0||gp.size()>tmp2.size()) {
																	gp=tmp2;																	
																}
																 nodeadend=true;
																 break;
															 } 
														 } //--End for 	pk													 
													 } //--Endif badnodes
													 if(nodeadend) break;													
												} //--End for pj	
												if(nodeadend) break;
										} //--End for pi
									} //--End no deadend	 
								} //--End for each node
								if (nodeadend) {
									//--We need the good dist here in real_dist_g2
									real_dist_g2=dist_path(gp,undirected_adjlist_g2_dist);
									//real_dist_g2=gp.size()-1;
									if (real_dist_g2==dist_g1[j]) {
										path_type=c_equal;
									} else {
										path_type=c_detour;
									}									
								} else {
									path_type=c_roadblock;
								}
								
											
					}
				} 
			 } //--End else
			#pragma omp atomic
			total++;			 		
			int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
			#pragma omp atomic
			total_time+=elapsed;
			if (verbose) {
				#pragma omp critical 
				{
					 string d1=SSTR(dist_g1[j]);
					 string d2=SSTR(real_dist_g2);
					 if (dist_g1[j]>=Inf) d1="Inf";
					 if (real_dist_g2>=Inf) d2="Inf";
					 FileOutput<< node_name[i] <<"\t"<< node_name[j]<<"\t"<<d1<<"\t"<<d2<<"\t";
					 //--cout << node_name[i] <<"\t"<< node_name[j]<<"\t"<<dist_g1[j]<<"\t"<<real_dist_g2<<"\t";				 
					 // switch(path_type) {
					// case c_breakthrough: cout<<"B"; break;
					// case c_roadblock:  cout<<"R"; break;
					// case c_impasse:  cout<<"I"; break;
					// case c_detour:  cout<<"D"; break;
					// case c_equal:  cout<<"E"; break;
					// case c_shortcut:  cout<<"S"; break;
					// }
					 switch(path_type) {
					case c_breakthrough: FileOutput<<"B"; break;
					case c_roadblock:  FileOutput<<"R"; break;
					case c_impasse:  FileOutput<<"I"; break;
					case c_detour:  FileOutput<<"D"; break;
					case c_equal:  FileOutput<<"E"; break;
					case c_shortcut:  FileOutput<<"S"; break;
					}
					add_statistics(i,j,path_type,gp);
					FileOutput<<"\t"<<elapsed<<"\t";
					// cout<<"\t"; 
					  fout_paths(gp);				
				}	
			}
			
			switch(path_type) {
				case c_breakthrough: B++; break;
				case c_roadblock: R++; break;
				case c_impasse: I++; break;
				case c_detour: D++; break;
				case c_equal: E++; break;
				case c_shortcut: S++; break;
			}
			//--Iterator for j
			j++;
		  if (j>=iter_n||random_path.size()>0) {
				more=false;
			}
		 } //-- End for j			
		
	 }	  
	 
    
	vector<int> results;
	results.push_back(B);
	results.push_back(R);
	results.push_back(I);
	results.push_back(D);
	results.push_back(E);
	results.push_back(S);
	results.push_back(total);
	results.push_back(total_time);
	return results;
}
	 
	 
	 
/*
 * MAIN
 */
int main(int nargc, char** argv) {
	constructor();
	//--test 
	
      if(readParameters(&param,argv,nargc)==-1){
             help();
			 destructor();
             exit(-1);
	 }
      cout<<description<<endl;
	  //--Random
	  if (param.seed!=-1) {
		  srand(param.seed); 
	  } else {
		  srand(time(NULL)); 
	  }
	  
      if (!param.found_g1) {
          cout<<"Error. Unable to locate the g1 file:"<<param.graph1<<endl;
		  destructor();
          exit(-1);
      }
	  
      if (!param.found_g2&&!param.found_attributes) {
          cout<<"Error. Unable to locate the g2 or attributes file g2:"<<param.graph2<<" attributes:"<<param.attributes<<endl;
		  destructor();
          exit(-1);
      }
	  if (param.verbose) {
		  FileOutput.open(param.outputfile);	  
		   if(!FileOutput.is_open())
		  {
			  printf("\n%s: verbose result file open failed...",param.outputfile);
			  destructor();
			  exit(1);
		  }
		
	  }
	  //--List option
	
    node_id_g1=load_edge(param.graph1, edge_list_g1);
	total_n_g1=node_id_g1.size();
	if (param.found_g2) {
		node_id_g2=load_edge(param.graph2, edge_list_g2);	
		total_n_g2=node_id_g2.size();
	}
	//--Test if we have the nodes of g1 in g2
	if (param.found_attributes) {
		  load_attributes(param.attributes);
	  }
	if (attributes.empty()&&!param.found_g2) {
		 cout<<"Error. Unable to load attributes file:"<<param.attributes<<endl;
		destructor();
		 exit(1);
	} else if (!attributes.empty()&&!param.found_g2) {
		 //copy g1 to g2
		node_id_g2=node_id_g1;
		total_n_g2=total_n_g1;
		edge_list_g2=edge_list_g1;
	}
	
	//--Create tha Adj. matrix
	build_adjacency_list();
   
	//--Seting some global variables
	total_paths=(total_nonk_node*(total_nonk_node-1.0))/2.0; 
	if (param.directed)  total_paths=1.0*total_nonk_node*(total_nonk_node-1);
	//-
	total_group=(total_paths/param.size)+1;
	if (fmod(total_paths,param.size)==0) total_group=total_group-1; //Correct for limit case
	
	
	if (param.random!=-1) {
		total_group=((param.random*total_paths)/param.size)+1;
		if (param.random>1.0) total_group=(param.random/param.size)+1;
		int to_create=(total_paths*param.random);
		if (param.random>1) to_create=param.random;
		if (to_create % param.size==0) total_group=total_group-1; 
		random_paths=createSamplePath(0,to_create,param.directed,true);		 
		//total_
	 }
	if (param.last==0||param.last>total_group) param.last=total_group;
	if (param.first==0||param.first<1) param.first=1;	
	output_header(argv);	
	 if (param.info) {
		//debug_info(); //--for debug
		destructor();
		exit(0);
	}
	time_t starttime = time(0);
	int B=0;
	int R=0;
	int I=0;
	int D=0;
	int E=0;
	int S=0;
	int total=0;
	int total_time=0;
	cout<<"\n============================ PARTIAL RESULTS ==================================\n";	
	cout<<"Group\tB\tR\tI\tD\tE\tS\tTotal\tCPU time (ms)"<<endl;
	if (param.verbose) {
		  FileOutput<<"\n============================ PARTIAL RESULTS ==================================\n";	
		  FileOutput<<"Node1\tNode2\tDist_X\tDist_Y\tType\tCPU time (ms)\tPath\tTaxa"<<endl;
	}
	
	#pragma omp parallel for ordered reduction(+:B,R,I,D,E,S,total,total_time) shared(random_paths)
	for (int grp=param.first; grp<=param.last;grp++) {
		
		vector<int> r=brides(grp);	
		B+=r[0];
		R+=r[1];
		I+=r[2];
		D+=r[3];
		E+=r[4];
		S+=r[5];
		total+=r[6];
		total_time+=r[7];
		#pragma omp critical 
		{		
			cout<<grp<<"\t"<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\t"<<r[4]<<"\t"<<r[5]<<"\t"<<r[6]<<"\t"<<r[7]<<"\t"<<endl;		
		}
	}
	#pragma omp critical 
	{		
		time_t endtime = time(0);
		double ttime=difftime(endtime, starttime);
		cout<<"\n================================= INFO =======================================\n";
		cout<<"(B) reakthrough : pathway impossible in network X but possible in network Y.\n";
		cout<<"(R) oadblock    : pathway possible in network X but impossible in network Y.\n";
		cout<<"(I) mpasse      : pathway impossible in both X and Y networks.\n";
		cout<<"(D) etour       : pathway shorter in network X than in network Y.\n";
		cout<<"(E) qual        : pathway of same length in networks X and Y.\n";
		cout<<"(S) hortcut     : pathway longer in network X than in network Y.\n";
		cout<<"\n================================ RESULTS ======================================\n";
		cout<<"\tB\tR\tI\tD\tE\tS\tTotal\tTime (s)"<<endl;	
		cout<<"\t"<<B<<"\t"<<R<<"\t"<<I<<"\t"<<D<<"\t"<<E<<"\t"<<S<<"\t"<<total<<"\t"<<ttime<<endl;
		cout<<"===============================================================================\n";	
		if (param.verbose) {
		    FileOutput<<"\n================================= INFO =======================================\n";
			FileOutput<<"(B) reakthrough : pathway impossible in network X but possible in network Y.\n";
			FileOutput<<"(R) oadblock    : pathway possible in network X but impossible in network Y.\n";
			FileOutput<<"(I) mpasse      : pathway impossible in both X and Y networks.\n";
			FileOutput<<"(D) etour       : pathway shorter in network X than in network Y.\n";
			FileOutput<<"(E) qual        : pathway of same length in networks X and Y.\n";
			FileOutput<<"(S) hortcut     : pathway longer in network X than in network Y.\n";
			FileOutput<<"\n============================== STATISTICS ====================================\n";
			FileOutput<<"Name (NK)\tB\tR\tI\tD\tE\tS\tInside\tAttribute"<<endl;
			for (int i=0; i<total_n;i++) {
				if (Nsu[i]&&!Nsb[i]) {
					 // string min1=SSTR(statistics[i][7]);
					 // string max1=SSTR(statistics[i][8]);
					 // // if (statistics[i][7]>=Inf) min1="Inf";
					 // if (statistics[i][8]>=Inf) max1="Inf";
					 // if (statistics[i][8]==-1) max1="Inf";
					 FileOutput<<node_name[i]<<"\t"<<statistics[i][0]<<"\t"<<statistics[i][1]<<"\t"<<statistics[i][2]<<"\t"<<statistics[i][3]<<"\t"<<statistics[i][4]<<"\t"<<statistics[i][5]<<"\t"<<statistics[i][6]<<"\t"<<attributes[i]<<endl;
				}
			}
			FileOutput<<endl<<"Name (K)\tInside\tAttribute"<<endl;
			
			for (int i=0; i<total_n;i++) {
				if (Nsb[i]) {
					FileOutput<<node_name[i]<<"\t"<<statistics[i][6]<<"\t"<<attributes[i]<<endl;
				}
			}
			FileOutput<<"*Inside: number of time a node is inside another path;K: K nodes;NK: non-K nodes\n";
			FileOutput<<"\n================================ RESULTS ======================================\n";
			FileOutput<<"\tB\tR\tI\tD\tE\tS\tTotal\tTime (s)"<<endl;
			FileOutput<<"\t"<<B<<"\t"<<R<<"\t"<<I<<"\t"<<D<<"\t"<<E<<"\t"<<S<<"\t"<<total<<"\t"<<ttime<<endl;
		}
	}
	if (param.verbose) FileOutput.close();
	
	destructor();
  
    return 0;
}

/*
 * Initialize some variables
 */
void constructor() {
     Nsb=new bool[maxnodes+1];
	 Nsu=new bool[maxnodes+1];
}


/*
 * Destructor
 */ 
void destructor() {    
    delete[] Nsb;
	delete[] Nsu;   
}


void build_adjacency_list() {
	total_k_node=0;
	for (int i=0; i<total_n;i++) {
		if (!param.K.empty()&&param.found_attributes) {
			//Is the node attributes part of K
			string attr=attributes[i];
			Nsb[i]=(std::find(param.K.begin(), param.K.end(), attr)!=param.K.end());
		} else {
			Nsb[i]=!(node_id_g1.count(i)>0);
		}
		//--Total
		if (Nsb[i]) {
			total_k_node++;
		} 
	}
	for (int i=0; i<total_n;i++) {
		Nsu[i]=true;
		if (!param.nonK.empty()&&param.found_attributes) {
			string attr=attributes[i];
			Nsu[i]=(std::find(param.nonK.begin(), param.nonK.end(), attr)!=param.nonK.end());		
		}
		if (Nsu[i]&&!Nsb[i]) {
			total_nonk_node++;
		}
	}
	//if (!param.K.empty()&&param.found_attributes) total_n_g1=total_nonk_node;
	for (int j=0; j<total_n; j++) {
			vector<int> tmp;
			vector<float> tmp2;
			std::map<int,float> tmp3;
			undirected_adjlist_g1.push_back(tmp);
			undirected_adjlist_g2.push_back(tmp);
			undirected_adjlist_g1_dist.push_back(tmp3);
			undirected_adjlist_g2_dist.push_back(tmp3);
		}

	for (int i=0; i<edge_list_g1.size();i++) {
		
		edge e=edge_list_g1[i];		
		if (e.to!=-1&&!Nsb[e.to]&&!Nsb[e.from]&&Nsu[e.to]&&Nsu[e.from]) {
			undirected_adjlist_g1[e.from].push_back(e.to);
			//cout<<node_name[e.from]<<" "<<node_name[e.to]<<" "<<e.dist<<endl;
			float dist=e.dist;
			if (param.inv_dist) dist=1.0f/dist;
			if (!param.use_dist) dist=1.0f;
			
			undirected_adjlist_g1_dist[e.from].insert(std::pair<int,float>(e.to,dist));
			if (!param.directed) {
				undirected_adjlist_g1[e.to].push_back(e.from);
				undirected_adjlist_g1_dist[e.to].insert(std::pair<int,float>(e.from,dist));
			}
		} 		
	}
	for (int i=0; i<edge_list_g2.size();i++) {
		edge e=edge_list_g2[i];
		if (e.to!=-1) {
			undirected_adjlist_g2[e.from].push_back(e.to);
			float dist=e.dist;
			if (param.inv_dist) dist=1.0f/dist;
			if (!param.use_dist) dist=1.0f;
			undirected_adjlist_g2_dist[e.from].insert(std::pair<int,float>(e.to,dist));
			//cout<<node_name[e.from]<<" "<<node_name[e.to]<<" "<<e.dist<<endl;
			if (!param.directed) {
				undirected_adjlist_g2[e.to].push_back(e.from);
				undirected_adjlist_g2_dist[e.to].insert(std::pair<int,float>(e.from,dist));
			}
		}  		
	}
	if (attributes.empty()) {
		for (int i=0;i<total_n;i++) {
			std::string s="0";
			if (Nsb[i]) s="1";
			attributes.insert(std::pair<int,std::string>(i,s));
		}
	}
	for (int i=0; i<total_n;i++) {
		vector<float> tmp(10,0.0f);
		statistics.push_back(tmp);		
		statistics[i][7]=Inf;
		statistics[i][8]=-1.0f;
	}
	
}



void load_attributes(char* filename) {
    
     std::ifstream file;
     std::string   line;
     
     file.open(filename);
     int count=0;
	 attributes.clear();
   try {
        while(std::getline(file, line))
        {
			if (line.find("#")!=0&&line.find("%")!=0&&line.length()>0) {
				std::stringstream   linestream(line);
				// std::getline(linestream, data, '\t');
				std::vector<std::string> tokens;
				std::string token;
				while(std::getline(linestream, token, '\t')) tokens.push_back(token); 
				int node_id=get_node(tokens[0]);
				if (node_id>-1) {
					attributes.insert(std::pair<int,string>(node_id,tokens[1]));
					//if (uniques_attributes.count(tokens[1])>0) {
					int count=uniques_attributes.count(tokens[1]);
					if (count==0) {
						uniques_attributes.insert(std::pair<string,int>(tokens[1],1));
					}	else {
						uniques_attributes.insert(std::pair<string,int>(tokens[1],uniques_attributes[tokens[1]]++));
					}
					
					//} else {
					//	uniques_attributes.insert(std::pair<string,int>(tokens[1],1));
					//}
					//attributes[node_id]=tokens[1];
				} 				
				//cout<<tokens[0]<<" "<<node_id<<" "<<tokens[1]<<endl;
			}
		}
	} catch(const std::exception& e) {cout<<"Error";}   
   file.close();
   // cout <<"Total: "<< attributes.size() << '\n';
    
}

int add_node(std::string name) {  
	//--Trim string
    if (name==""||name=="-") return -1;
	std::map<std::string,int>::iterator search =node_name_lookup.find(name);
    if (search==node_name_lookup.end()) {                
        node_name.push_back(name);    
           node_name_lookup[name]=current_node;                   
           current_node++;
		   total_n++;
           return (current_node-1);
       } else {
           return (search->second);
       }
}

int get_node(std::string name) {
      std::map<std::string,int>::iterator search =node_name_lookup.find(name);
    if (search==node_name_lookup.end()) {        
           //--Not found
           return (-1);
       } else {
           return (search->second);
       }
}

std::map<int,int> load_edge(char* filename, std::vector<edge> &edge_list) {
     std::ifstream file;
     std::string   line;
     std::string   data;
	 std::map<int,int> local_node;
	 edge_list.clear();
	 file.open(filename);
    
        std::string                 from;
        std::string                 to;
        int dist;

   try {
        while(std::getline(file, line))
        {
            std::stringstream   linestream(line);
            // std::getline(linestream, data, '\t');
            std::vector<std::string> tokens;
            std::string token;
            if (line.find("#")!=0&&line.find("%")!=0&&line.length()>0) {
                while(std::getline(linestream, token, '\t')) tokens.push_back(token); 
                //CASE 1. WE have space
				from="";
				to="";
				if (tokens.size()<2) {						
					std::stringstream linestream2(line);
					tokens.clear();
					linestream2 >> from >> to >> dist;
				} else {
				//CASE 2. WE have tab (default)
					from=tokens[0];
					to=tokens[1];					
				}
				dist=1.0f;
                if (tokens.size()>2) dist=atoi(tokens[2].c_str());
				if (dist<=0.0f) dist=1.0f; //--Don't permit negative distance
                 edge t;             
                 t.from=add_node(from);
                 t.to=add_node(to);
				 t.dist=dist; 
				 
				 local_node.insert ( std::pair<int,int>(t.from,0));
				 if (t.to!=-1) local_node.insert ( std::pair<int,int>(t.to,0));
				 
                 if (t.to!=t.from) edge_list.push_back(t);
            }
        }
   } catch(const std::exception& e) {}   
   cout<<"Done loading "<<filename<<"..."<<endl;
   file.close();
  return(local_node);
} 

int randomInt (int min, int max){
    int n = rand()%(max-min)+min; 
    return n;
}

/**
 * output some information to the file
 * Format will be similar to Nexus
 */
void output_header(char** argv) {
 
    //cout<<startMessage; 
	
	cout<<"\n=============================== PARAMETERS ====================================\n";		
	
	//-- OpenMP
	  cout<<"\n-=[Multicores support]=-"<<endl;
	  cout<<"Maximum threads  : "<<omp_get_max_threads()<<endl;
	  cout<<"Number of cores  : "<<omp_get_num_procs()<<endl;
	 //--Loaded network 
	  cout<<"\n-=[Input files]=-"<<endl;
	  if (param.directed) cout<<"Networks         : directed"<<endl;      
	  if (!param.directed) cout<<"Networks         : undirected"<<endl;      
	  cout<<"NetworkX         : "<<param.graph1<<endl;
      if (param.found_g2) cout<<"NetworkY         : "<<param.graph2<<endl;
      if (param.found_attributes) cout<<"Attributes       : "<<param.attributes<<endl; 
	  cout<<"Nodes in networkX: "<<total_nonk_node<<endl;   
      cout<<"Nodes in networkY: "<<total_n_g2<<endl; 
      if (param.found_attributes)  {
		  cout<<"Attributes|count : "<<uniques_attributes.size()<<endl;   
		  for(map<string,int>::iterator it = uniques_attributes.begin(); it != uniques_attributes.end(); ++it) {
			  cout<<"\t"<<it->first<<"|"<<it->second<<endl;
		  }
	  }	   
      if (!param.nonK.empty()) cout<<"non-K nodes attr.: "<<param.nonK<<endl;
	  cout<<"Total K nodes    : "<<total_k_node<<endl; 
	  if (!param.K.empty()) cout<<"K nodes attr.    : "<<param.K<<endl;
	  cout<<"Total paths      : "<<total_paths<<endl; 
      cout<<"\n-=[Run parameters]=-"<<endl;
	  if (param.random!=-1) {		
		cout<<"Running mode     : random"<<endl;  
		if (param.random<1.0) {
			cout<<"Investigated     : "<<(param.random*100)<<"% ("<<(int)(total_paths*param.random)<<")"<<endl; 
		} else {
			cout<<"Investigated     : "<<param.random<<endl;  
		}
	  } else {
	  cout<<"Running mode     : normal"<<endl;
	  }
	  if (param.use_dist&&param.inv_dist) {
		cout<<"Edge distance    : inverse weights"<<endl;      	  
	  } else if (param.use_dist) {
		cout<<"Edge distance    : edge weights"<<endl;    
	  } else {
		cout<<"Edge distance    : no weights used"<<endl;     
	  }
	  cout<<"Group size       : "<<param.size<<endl;      
	  cout<<"Number of groups : "<<total_group<<endl;      
	  cout<<"First group      : "<<param.first<<endl;      
      cout<<"Last group       : "<<param.last<<endl;
      cout<<"Maxdistance      : "<<param.maxdistance<<endl;       
      cout<<"Maxnode          : "<<param.maxnode<<endl;       	  
	  cout<<"Maxtime (s)      : "<<(param.maxtime/1000)<<endl;    
	  if (param.max_individual_path!=100) cout<<"Maxpath          : "<<param.max_individual_path<<endl; 
	  
	  cout<<"\n-=[Miscellaneous]=-"<<endl;
 	  if (param.seed<0) {
		  cout<<"Seed             : clock time"<<endl;      
	  } else {
		  cout<<"Seed             : "<<param.seed<<endl;      
	  }  
	  if (param.verbose) 
		  cout<<"Output file      : "<<param.outputfile<<endl;      
	  cout<<"\n===============================================================================\n";	
      
	  if (param.verbose) {
		  FileOutput<<"\n=============================== PARAMETERS ====================================\n";		
	
			//-- OpenMP
			  FileOutput<<"\n-=[Multicores support]=-"<<endl;
			  FileOutput<<"Maximum threads  : "<<omp_get_max_threads()<<endl;
			  FileOutput<<"Number of cores  : "<<omp_get_num_procs()<<endl;
			 //--Loaded network 
			  FileOutput<<"\n-=[Input files]=-"<<endl;
			  if (param.directed) FileOutput<<"Networks         : directed"<<endl;      
			  if (!param.directed) FileOutput<<"Networks         : undirected"<<endl;      
			  FileOutput<<"NetworkX         : "<<param.graph1<<endl;
			  if (param.found_g2) FileOutput<<"NetworkY         : "<<param.graph2<<endl;
			  FileOutput<<"Nodes in networkX: "<<total_nonk_node<<endl;   
			  FileOutput<<"Nodes in networkY: "<<total_n_g2<<endl; 			  
			  if (param.found_attributes)  {
				  FileOutput<<"Attributes|count : "<<uniques_attributes.size()<<endl;   
				  for(map<string,int>::iterator it = uniques_attributes.begin(); it != uniques_attributes.end(); ++it) {
					  FileOutput<<"\t"<<it->first<<"|"<<it->second<<endl;
				  }
			  }	   
			 
			  
			  
			  if (!param.nonK.empty()) FileOutput<<"non-K nodes attr.: "<<param.nonK<<endl;
			  FileOutput<<"Total K nodes    : "<<total_k_node<<endl; 
			  if (!param.K.empty()) FileOutput<<"K nodes attr.    : "<<param.K<<endl;
			  FileOutput<<"Total paths      : "<<total_paths<<endl; 
      
			  //--Running parameters
			  FileOutput<<"\n-=[Run parameters]=-"<<endl;
			  if (param.random!=-1) {		
				FileOutput<<"Running mode     : random"<<endl;  
				if (param.random<1.0) {
					FileOutput<<"Investigated     : "<<(param.random*100)<<"% ("<<(int)(total_paths*param.random)<<")"<<endl; 
				} else {
					FileOutput<<"Investigated     : "<<param.random<<endl;  
				}
			  } else {
			    FileOutput<<"Running mode     : normal"<<endl;
			  }
			  if (param.use_dist&&param.inv_dist) {
				FileOutput<<"Edge distance    : inverse weights"<<endl;      	  
			  } else if (param.use_dist) {
				FileOutput<<"Edge distance    : edge weights"<<endl;    
			  } else {
				FileOutput<<"Edge distance    : no weights used"<<endl;     
			  }
			  FileOutput<<"Group size       : "<<param.size<<endl;      
			  FileOutput<<"Number of groups : "<<total_group<<endl;      
			  FileOutput<<"First group      : "<<param.first<<endl;      
			  FileOutput<<"Last group       : "<<param.last<<endl;
			  FileOutput<<"Maxdistance      : "<<param.maxdistance<<endl;       
			  FileOutput<<"Maxnode          : "<<param.maxnode<<endl;       	  
			  FileOutput<<"Maxtime (s)      : "<<(param.maxtime/1000)<<endl;      
			  if (param.max_individual_path!=100) FileOutput<<"Maxpath          : "<<param.max_individual_path<<endl; 
			  FileOutput<<"\n-=[Miscellaneous]=-"<<endl;
			  if (param.seed<0) {
				  FileOutput<<"Seed            : clock time"<<endl;      
			  } else {
				  FileOutput<<"Seed             : "<<param.seed<<endl;      
			  }  
			  if (param.verbose) 
				  FileOutput<<"verbose          : "<<param.outputfile<<endl;  
			  FileOutput<<"\n===============================================================================\n";	
	  }
			
	
}

void help(){
	printf("\n%s",description);	
	printf("\nUsage :\nBRIDES -X=[filename] -Y=[filename] -outputfile=[filename]\n");
	printf("       -maxdistance=[1..n-1] -maxnodes=[1..n-1]\n");
	//printf("       -kmax=[0... n-1] -replicate=[1..n] -inputtype=[0,1]\n");

         printf("\nParameters :");
        printf("\n-X=file           [filename for network X]");
        printf("\n-Y=file           [filename for network Y]");
		printf("\n-attributes=file  [filename for node attributes]");
		printf("\n-usedist          [Use edge distances found in network files]");
		//printf("\n-invdist          [Use the inverse of edge distances found in network files]");
		printf("\nK=B               [attributes to considers as K separated by comma e.g. A,B,C]");
        printf("\nNK=A              [attributes to considers as non-K]");
		printf("\n-random=XXX       [sample XXX random pathways]");
		printf("\n-first=1          [first group of path to process]");
        printf("\n-last=n           [last group of path to process]");
		printf("\n-size=1000        [group size, default                         : 1000]");      
        printf("\n-maxdistance=100  [maximum path length to search, default      : 100]");
		printf("\n-maxnodes=100     [maximum augmented node (K) to search,default: 100]");
		printf("\n-maxtime=1        [maximum time for each path search, default  : 1 second]");
        printf("\n-output=file      [output to file each path information: taxa, distance, etc.]");
		printf("\n-seed=999         [set the random seed generator to a specific seed]");
		printf("\n\nExample : \n: ./BRIDES -X=sample_g1.txt -Y=sample_g2.txt\n\n");
}


int ExtraireDonnees(const char * chaine, char *champs, char * contenu){

	int cpt=0,i;
	int egale=false;
	int tailleChaine;

	if(chaine[0] != '-')
		return false;
	tailleChaine = (int)strlen(chaine);

	for(i=1;i<tailleChaine;i++){

		if (chaine[i] == '='){
			egale = true;
			champs[cpt] = '\0';
			cpt=0;
			continue;
		}
		if(egale)
			contenu[cpt++] = chaine[i];
		else
			champs[cpt++] = chaine[i];
	}
	
	if (!egale) {
	    contenu[0] = '\0';
		champs[cpt] = '\0';	
	} else {
		contenu[cpt] = '\0';
	}

	return true;
}

//===================================================================================
//= REAWD PARAMETERS FROM COMMAND LINE
//===================================================================================
int readParameters(Parameters *param, char **argv, int nargc){

	char champs[100];
	char contenu[100];
	int i;	
	//default 
        if (nargc<2) {               
               return -1;
           }
       // sprintf((*param).inputfile,"input.txt");
        (*param).maxdistance=100.0f; //--maxdistance of K nodes
		(*param).maxnode=100; //--Maxnode K to investigate
		(*param).maxtime=10000; //--Maxtime in millsecond (10s per path)
		(*param).size=1000; //-- default group size
		(*param).max_individual_path=100; //Since we keep all shortest path, this is the number of shortest path investigated for each path
        (*param).seed=-1;
		(*param).first=0;
		(*param).last=0;
		(*param).K.clear();
		(*param).nonK.clear();
		(*param).random=-1;		
         sprintf((*param).graph1,"%s","");
         sprintf((*param).graph2,"%s","");
		 sprintf((*param).attributes,"%s","");
         sprintf((*param).outputfile,"%s","");		 
		 sprintf((*param).debugfile,"%s","");
         (*param).info=false;
		 (*param).directed=false;         
		 (*param).verbose=false;         
		 (*param).found_g1=false;
         (*param).found_g2=false;
         (*param).found_attributes=false;
		 (*param).use_dist=false;
		 (*param).inv_dist=false;
         //(*param).use_g2=false;
         //sprintf((*param).outputfile,"%s","output.txt");
		/*(*param).distance_method=0;  //default cosine distance
        (*param).type_method=0; //default k-means
        (*param).inputtype=1; // default row are object - column variable
  		(*param).use_weight=false;
        (*param).kmax=0;
        (*param).replicate=1; //Default 1 replicate*/
        
	for(i=1;i<nargc;i++){

		if(ExtraireDonnees(argv[i],champs,contenu)){
 
			//--debug printf("\n : [%s] -> champs = %s , contenu = %s",argv[i],champs,contenu);

			if(strcmp("help",champs) == 0) {
                            return -1;
            }
			//======== debug ==============
			else if(strcmp("verbose",champs) == 0||strcmp("output",champs) == 0){
							sprintf((*param).outputfile,"%s",contenu);
							(*param).verbose=true;							
            }  
            //======= maxnode ==============       
			 else if(strcmp("maxnode",champs) == 0){                                
				  (*param).maxnode = atoi(contenu);
			}              
            //======= maxdistance ==============       
			 else if(strcmp("maxdistance",champs) == 0){                                
				  (*param).maxdistance = atof(contenu);
			}
			  //======= size ==============       
                         else if(strcmp("size",champs) == 0){                                
                              (*param).size = atoi(contenu);
			}
            //======= start ==============       
                         else if(strcmp("first",champs) == 0){                                
                              (*param).first = atoi(contenu);
			}
            //======= end ==============       
                         else if(strcmp("last",champs) == 0){                                
                              (*param).last = atoi(contenu);
			}
			//======= K ================
			 else if(strcmp("K",champs) == 0){    
			  std::stringstream   linestream(contenu);
             std::string token;
                while(std::getline(linestream, token, ',')) (*param).K.push_back(token); 
			}
			//======= non-K ================
			 else if(strcmp("nK",champs) == 0||strcmp("NK",champs) == 0){    
			  std::stringstream   linestream(contenu);
             std::string token;
                while(std::getline(linestream, token, ',')) (*param).nonK.push_back(token); 
			}
             //======= usedist ==============       
            else if(strcmp("usedist",champs) == 0){                                
                              (*param).use_dist = true;
			}
			 //======= invdist ==============       
            else if(strcmp("invdist",champs) == 0){                                
                              (*param).inv_dist = true;
							  (*param).use_dist = true;
			}
			//======= random path ==============       
            else if(strcmp("random",champs) == 0){
				(*param).random = atoi(contenu);
				if (strstr( contenu, "%" )!=NULL) (*param).random/=100;
			}
			 //======= SEED ==============       
            else if(strcmp("seed",champs) == 0){                                
                              (*param).seed = atoi(contenu);
			}
			 //======= MAXTIME ==============       
            else if(strcmp("maxtime",champs) == 0){                                
                              (*param).maxtime = atoi(contenu)*1000;
			}
			//======= MAXPATH ==============       
            else if(strcmp("maxpath",champs) == 0){                                
                              (*param).max_individual_path = atoi(contenu);
			}
			//======== input file ==============
			else if(strcmp("g1",champs) == 0||strcmp("X",champs) == 0||strcmp("x",champs) == 0){
                                sprintf((*param).graph1,"%s",contenu);	
                                FILE *Input1;
                                if ((Input1 = fopen((*param).graph1,"r"))!=0) {(*param).found_g1=true;}
                                fclose(Input1);
			}
            else if(strcmp("g2",champs) == 0||strcmp("Y",champs) == 0||strcmp("y",champs) == 0){
                                sprintf((*param).graph2,"%s",contenu);	
                               // (*param).use_g2=true;
                                 FILE *Input1;
                                if ((Input1 = fopen((*param).graph2,"r"))!=0) {(*param).found_g2=true;}
                                fclose(Input1);
			}
			 else if(strcmp("attributes",champs) == 0||strcmp("attr",champs) == 0||strcmp("A",champs) == 0){
                                sprintf((*param).attributes,"%s",contenu);	
                               // (*param).use_g2=true;
                                 FILE *Input1;
                                if ((Input1 = fopen((*param).attributes,"r"))!=0) {(*param).found_attributes=true;}
                                fclose(Input1);
			}
           //======== input file ==============
			/*else if(strcmp("weight",champs) == 0){
                            sprintf((*param).inputweightfile,"%s",contenu);				
                                (*param).use_weight=true;
			}*/
			//============ output file ===============
			// else if(strcmp("outputfile",champs) == 0){
				// sprintf((*param).outputfile,"%s",contenu);	
                                // //printf("%s\n",(*param).outputfile );
			// }
			//============= version ===============
			else if(strcmp("directed",champs) == 0){
					(*param).directed = true;
			}
			//============= info   ===============
			else if(strcmp("info",champs) == 0){
					(*param).info = true;
			}
			//============= version ===============
			else if(strcmp("version",champs) == 0){
                            return -1;
			}
			else {
				printf("\nBRIDES : invalid syntax for parameter : %s (needed -parameter=value)\n" , argv[i]);
							return -1;
			}
             }  //--end ExtraireDonnees         
        } //--end for
	
	
	if(strcmp((*param).graph1,"")==0){
		printf("\nerror : there is no input graph file");
		return -1;
	}
        
        if((*param).verbose&&strcmp((*param).outputfile,"")==0){
		 sprintf((*param).outputfile,"%s%s",(*param).graph1,".output.txt");		
	   }
	   //--Default debug
         if(strcmp((*param).debugfile,"")==0) {
			sprintf((*param).debugfile,"%s%s",(*param).graph1,".debug.txt");		 
		 } 
	return 0;
}