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
   #define omp_set_dynamic(i) 
   #define omp_set_num_threads(i)
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
                          "| This program implements a fast algorithm for characterizing evolving        |\n"
						  "| similarity networks using Breakthroughs, Roadblocks, Impasses, Detours,     |\n"
						  "| Equals and Shortcuts.                                                       |\n"
                          "===============================================================================\n";
		
const char *startMessage = "BRIDES V.1.0 - (2016) by Etienne Lord, V. Makarenkov, F-J. Lapointe\n"; 

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
    edge(int f, int t) {dist=1.0f;from=f;to=t;}
	
    bool operator<(const edge& rhs) const
    {
        return dist < rhs.dist;
    }
	
	
};
   
struct Parameters {
		char graph1[1024];       // filename g1, X
	char graph2[1024];       // filename g2,Y
	char attributes[1024];   // filename attributes
	char graphinfo[1024];    // Node appartenance
	char outputfile[1024];   // outputfile	
	char debugfile[1024];    // debugfile
	vector<std::string> K;   // attributes (string)to use for K nodes
	vector<std::string> nonK;// attributes (string) to use for non K nodes 


		bool found_g1;         //-- found g1 (networkX)file?
        bool found_g2;         //-- found g2 (networkY) file?
        bool found_attributes; //-- found Annotations file?
        bool info;             //-- Output debug information

		bool use_dist;   //--by default false, (dist set to 1, unweighted)
		bool inv_dist;  //-- use the inverse of the distance
        bool directed;  //-- Network are directed 
		bool removeK_from_X; //Remove k node from networkX
		
		int seed;       //-- Random seed
		float random;   //-- Number of random path to classify
		int size;       //-- Size of group 

        int first;      //-- Start group of vertex (defined by size)
        int last;       //-- End group of vertex
		int strategy;   //-- K Node ordering strategy (either 1- max(d(i,k), d(j,k)) or 2- (d(i,k)+d(j,k)))

        float maxdistance;       //--default search length (default=100);		
		int maxtime;             //--Max time to search for Detour or Equal path (default=10 s)
        float maxnode;           //--Max k node to investigate (default=100)


		int max_individual_path; //--Max path number (default=100)
		int maxthread;           //--Default: system 
		bool verbose;            //-- Output more result to file              

		int heuristic; //--BRIDES =1, BRIDES_YEN=2, BRIDES_YC=3, BRIDES_EC=4. DFS=5, BRIDES RELAXED=6
		int weight_method; //All weight as equal, inverse, default
		bool use_multiple_color; // TO DO
		int start_time_code;   //Time code for including vertex (Unix timestamps)
		int end_time_code;     //Time code for including vertex (Unix timestamps)
		bool output_graph;      
		bool use_time_code;
      
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
void dijkstra (int S, vector< vector<int> > adj, vector< std::map<int,float> > adj_dist,vector<float> &dist, vector< map<int,int> > &prev,  map<int,int> ignored_nodes, map<edge,int> ignored_edges);

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
bool good_path(vector<int> path, int S, int T);
void output_graph();
 
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
map<int,int> intersect2(vector<int> v1, vector<int> v2);

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
int total_nonk_node=0; //--should be equal to g1 normally
	
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

std::istream& getline(std::istream& is, std::string& t)
{
    t.clear();
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    for(;;) {
        int a = sb->sbumpc();
        switch (a) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:         
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)a;
        }
    }
}

void debug_info() {
	// Calculate some informations
	// Link type 
	cout<<"name\tid\tgroup\tattr\tnode1_cnt\tnode2_cnt"<<endl;
	for (int i=0; i<total_n;i++) {
		cout<<node_name[i]<<"\t"<<i<<" \t"<<Nsb[i]<<"\t"<<attributes[i]<<"\t"<<node_id_g1.count(i)<<"\t"<<node_id_g2.count(i)<<endl;
	}
	cout<<" NetworkX:"<<endl;
	for (int i=0; i<edge_list_g1.size();i++) {
			if (edge_list_g1[i].to!=-1) {
				cout<<"|"<<node_name[edge_list_g1[i].from]<<"|-|"<<node_name[edge_list_g1[i].to]<<"|"<<endl;
			} else {
				cout<<"|"<<node_name[edge_list_g1[i].from]<<"|- "<<endl;
			}
	}
	cout<<" NetworkY:"<<endl;
	for (int i=0; i<edge_list_g2.size();i++) {
			if (edge_list_g2[i].to!=-1) {
				cout<<"|"<<node_name[edge_list_g2[i].from]<<"|-|"<<node_name[edge_list_g2[i].to]<<"|"<<endl;
			} else {
				cout<<"|"<<node_name[edge_list_g2[i].from]<<"|- "<<endl;
			}
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
bool found_edges(int S, int T, map<edge,int> edges) {
	
	for (std::multimap<edge,int>::iterator it=edges.begin(); it!=edges.end();++it) {
		if (((*it).first.from==S&&(*it).first.to==T)||(!param.directed&&(*it).first.to==T&&(*it).first.from==S)) return true;
	}
	return false;
}
 
void dijkstra (int S, vector< vector<int> > adj, vector< std::map<int,float> > adj_dist,vector<float> &dist, vector< map<int,int> > &prev,  map<int,int> ignored_nodes, map<edge,int> ignored_edges) {
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
							 //cout<<u<<"-"<<v<<endl;
							 
							 if (found_edges(u,v,ignored_edges)) {
								continue;
								
							 }
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

 // use DFS to get BRIDES
     void dfs_BRIDES(int s, int t, vector< vector<int> > adj, vector< vector<int> > &paths, vector<int> &path, vector<bool> &onPath ) {		         
		if (paths.size()>0&&path.size()>paths[0].size()) return;
		if (path.size()>param.maxdistance) return;
	
		// found path from s to t
        if (s == t) {    			
			vector<int> tmp=path;
			if (good_path(path,path[0],t)) {			
				//We kept only the smaller one...
				if (paths.size()>0&&tmp.size()<paths[0].size()) {
					paths.clear();
				}	
				paths.push_back(tmp);			
			}
			return;
        } 
		  for(int i=0; i< adj[s].size();i++) {
			 int w=adj[s][i]; 	
			vector<int> tmp=path;			 
			if (!onPath[w]) {
				path.push_back(w);
				onPath[w] = true;
				dfs_BRIDES(w,t,adj, paths, path, onPath);
				path.pop_back();
				onPath[w] = false; 
			}
		  }		           	
    }

	
vector<int> get_path_dfs_BRIDES(int S, int T, vector< vector<int> > adj) {
	vector< vector<int> > paths;
	vector<bool> onPath(adj.size(), false);
	vector<int> path;
	path.push_back(S);
	onPath[S] = true;
	dfs_BRIDES(S,T,adj,paths, path,onPath);		
	if (paths.size()>0) return paths[0];
	path.clear();
	return path;
}	
	 
 // use DFS to get Path from the Dijkstra prev matrix
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

 // use DFS for YEN
     void dfs_pathYEN(int s, int t,vector< map<int,int> > prev, vector< vector<int> > &paths, vector<int> &path, vector<bool> &onPath ) {		 
        if (paths.size()>1) return;
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
				  dfs_pathYEN(s, w, prev,paths, path, onPath);
				  path.pop_back();
				onPath[w] = false; 
			}
		  }		           	
    }

//--This will output the path for the Dijkstra version of YenSPK
vector<int> get_pathYen(int S, int T, vector< map<int,int> > prev) {
	vector< vector<int> > paths;	
	vector<bool> onPath(prev.size(), false);
	onPath[T]=true;
	vector<int> path;	
	path.push_back(T);
	dfs_pathYEN(S,T, prev,paths, path,onPath);		
	
	for (int i=0; i<paths.size();i++ ) {
		std::reverse(paths[i].begin(), paths[i].end());
	}	
	if (paths.size()>0) return paths[0];
	path.clear();
	return path;
}
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
			len+=adj_dist[src][dest];
			src=dest;
		}
	return len;
}	 
	 
/**
 * This verify that a path is good i.e. contains one k-node
  * If the algo is relaxed (param.heuristic=6)
 */
bool good_path(vector<int> path, int S, int T) {
		bool found=false;
		bool found_S=false;
		bool found_T=false;
		map<int,int> ht;
		if (param.heuristic==6) {
			found=true;
			for (int i=0; i<path.size();i++) {
				if (ht[path[i]]!=0) return false; //Duplicate node
				ht[path[i]]=1;
				if (path[i]==T) found_T=true;
				if (path[i]==S) found_S=true;
			}
		
		} else {
			if (path.size()<3) return false;
			// any duplicate and path contains one k nodes
			for (int i=0; i<path.size();i++) {
				if (Nsb[path[i]]) found=true;
				if (ht[path[i]]!=0) return false;
				ht[path[i]]=1;
				if (path[i]==T) found_T=true;
				if (path[i]==S) found_S=true;
			}
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

//--Verify that a path doesn't contains the excluded node
bool good_path3(vector< vector<int> > paths, int exclude_node) {
	vector<int> tmp;
	for (int i=0; i<paths.size();i++) {
		for (int j=0; j<paths[i].size();j++) {
		 if (paths[i][j]==exclude_node) return false; 
		}
	}
	return true;
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

vector<int> join2(vector<int> v1, vector<int> v2) {
	vector<int> tmp;
		
	for (int i=0; i<v1.size();i++) tmp.push_back(v1[i]);
	for (int i=1; i<v2.size();i++) tmp.push_back(v2[i]);
	return tmp;
}
map<int,int> intersect(vector<int> v1, vector<int> v2,int k,int i2) {
	map<int,int> tmp;
	map<int,int> found;
	for (int i=0;i<v1.size();i++) {
		if (v1[i]!=k&&v1[i]!=i2) found.insert(std::pair<int,int>(v1[i],1));
	}
	for (int i=0;i<v2.size();i++) {
		if (found[v2[i]]>0) {			
		     tmp.insert(std::pair<int,int>(v2[i],1));
		} 
	}
	return tmp;
}

map<int,int> intersect2(vector<int> v1, vector<int> v2) {
	map<int,int> tmp;
	map<int,int> found;
	for (int i=0;i<v1.size();i++) {
		found.insert(std::pair<int,int>(v1[i],1));
	}
	for (int i=0;i<v2.size();i++) {
		if (found[v2[i]]>0) {			
		     tmp.insert(std::pair<int,int>(v2[i],1));
		} 
	}
	return tmp;
}

vector<int> extract(vector<int> path, int idx1, int idx2) {
	vector<int> tmp;
	for (int i=idx1; i<=idx2;i++)  tmp.push_back(path[i]);
	return tmp;
}

bool equal_vector(vector<int> path1, vector<int> path2) {
	if (path1.size()!=path2.size()) return false;
	return std::equal ( path1.begin(), path1.end(), path2.begin() );
}


vector< vector<int> > YenKSP(int S, int T,int K,vector< vector<int> > adj, vector< std::map<int,float> > adj_dist, vector<float> &dist) {

		    int len=adj.size();
			
			vector< map<int,int> > prev(len);
			vector<float> tmp_dist(len);
			prev.clear();			
			for (int i=0; i<len; i++) {
				map<int,int> tmp;
				prev.push_back(tmp);
			}
			std::fill(tmp_dist.begin(), tmp_dist.end(), 1e14);
		   
			dijkstra(S,adj,adj_dist,tmp_dist,prev);
			multimap<float,vector<int> > B;
			
			vector< vector<int> > A;
			vector<float> A_dist;
			A_dist.push_back(tmp_dist[T]);
			if (tmp_dist[T]>=Inf) {
				dist=A_dist;
				return A;
			}
			A.push_back(get_pathYen(S,T,prev));
			if (A.size()==0) { 
				
				dist.clear();
				dist.push_back(Inf);
				return A;
			}
			int k=1;
			
			while (k<=K) {
				map<edge,int> ignored_edges;
				map<int,int> ignored_nodes;
				for (int i=0; i<A[k-1].size()-1;i++) {
					int spurNode=A[k-1][i];
					vector<int> rootPath=extract(A[k-1],0,i);
					 ignored_edges.clear();
					 ignored_nodes.clear();
					 //--Remove edges
					for (int pi=0; pi<A.size();pi++) {
						//A[pi].size()>i
						if (A[pi].size()>i&&equal_vector(rootPath,extract(A[pi],0,i))) {
							edge e(A[pi][i],A[pi][i+1]);
							ignored_edges.insert(std::pair<edge,int>(e,1));
						}
					}
					//--Remove node
					for (int i=0; i<rootPath.size();i++) {
						int node=rootPath[i];
						if (node!=spurNode) ignored_nodes.insert(std::pair<int,int>(node,1));
					}
					//--Calculate the new spurPath
					prev.clear();			
					for (int i=0; i<len; i++) {
						map<int,int> tmp;
						prev.push_back(tmp);
					}
					vector<float> tmp_dist2(len);
					 std::fill(tmp_dist2.begin(), tmp_dist2.end(), 1e14);
					dijkstra(spurNode,adj,adj_dist,tmp_dist2,prev,ignored_nodes, ignored_edges);
					
					vector<int> spurPath=get_pathYen(spurNode,T,prev);
					if (intersect2(rootPath,spurPath).size()==1) {
						
						vector<int> newpath=join2(rootPath,spurPath);
						
						float newdist=tmp_dist[spurNode]+tmp_dist2[T];
						bool found=false;
						for (int i=0;i<A.size();i++) {
							if (equal_vector(A[i], newpath)) found=true;
						}
						for (std::multimap<float,vector <int> >::iterator it=B.begin(); it!=B.end(); ++it) {
							if (equal_vector((*it).second, newpath)) {
								found=true; 
							}
						}
						if (!found) B.insert(multimap< float,vector<int> >::value_type(newdist, newpath));
					}
				}
				
				if (B.size()==0) break;
				if ((B.begin())->second.size()>0) {
					A.push_back((B.begin())->second);
					A_dist.push_back((B.begin())->first);
					B.erase(B.begin());
				}
				
				k++;
				
				
						
			}
			dist=A_dist;
	return A;

}

vector< vector<int> > YenKSP_Stop(int S, int T,vector< vector<int> > adj, vector< std::map<int,float> > adj_dist, vector<float> &dist) {

		    int len=adj.size();
			int K=param.max_individual_path;
			vector< map<int,int> > prev(len);
			vector<float> tmp_dist(len);
			prev.clear();			
			for (int i=0; i<len; i++) {
				map<int,int> tmp;
				prev.push_back(tmp);
			}
			std::fill(tmp_dist.begin(), tmp_dist.end(), 1e14);
		   
			dijkstra(S,adj,adj_dist,tmp_dist,prev);
			multimap<float,vector<int> > B;
			
			vector< vector<int> > A;
			vector<float> A_dist;
			A_dist.push_back(tmp_dist[T]);
			if (tmp_dist[T]>=Inf) {
				dist=A_dist;
				return A;
			}
			A.push_back(get_pathYen(S,T,prev));
			if (A.size()==0) { 				
				dist.clear();
				dist.push_back(Inf);
				return A;
			}
			int k=1;
			
			while (k<=K) {
				map<edge,int> ignored_edges;
				map<int,int> ignored_nodes;
				for (int i=0; i<A[k-1].size()-1;i++) {
					int spurNode=A[k-1][i];
					vector<int> rootPath=extract(A[k-1],0,i);
					 ignored_edges.clear();
					 ignored_nodes.clear();
					 //--Remove edges
					for (int pi=0; pi<A.size();pi++) {
						if (A[pi].size()>i&&equal_vector(rootPath,extract(A[pi],0,i))) {
							edge e(A[pi][i],A[pi][i+1]);
							ignored_edges.insert(std::pair<edge,int>(e,1));
						}
					}
					//--Remove node
					for (int i=0; i<rootPath.size();i++) {
						int node=rootPath[i];
						if (node!=spurNode) ignored_nodes.insert(std::pair<int,int>(node,1));
					}
					//--Calculate the new spurPath
					prev.clear();			
					for (int i=0; i<len; i++) {
						map<int,int> tmp;
						prev.push_back(tmp);
					}
					vector<float> tmp_dist2(len);
					 std::fill(tmp_dist2.begin(), tmp_dist2.end(), 1e14);
					dijkstra(spurNode,adj,adj_dist,tmp_dist2,prev,ignored_nodes, ignored_edges);
					
					vector<int> spurPath=get_pathYen(spurNode,T,prev);
					if (intersect2(rootPath,spurPath).size()==1) {
						
						vector<int> newpath=join2(rootPath,spurPath);
						
						float newdist=tmp_dist[spurNode]+tmp_dist2[T];
						bool found=false;
						for (int i=0;i<A.size();i++) {
							if (equal_vector(A[i], newpath)) found=true;
						}
						for (std::multimap<float,vector <int> >::iterator it=B.begin(); it!=B.end(); ++it) {
							if (equal_vector((*it).second, newpath)) {
								found=true; 
							}
						}
						if (!found) B.insert(multimap< float,vector<int> >::value_type(newdist, newpath));
					}
				}
				
				if (B.size()==0) break;
				if ((B.begin())->second.size()>0) {
					A.push_back((B.begin())->second);
					A_dist.push_back((B.begin())->first);
					//STOP if we have a knodes
					if (good_path(A[k],S,T)) {
						return A;
					}
					B.erase(B.begin());
				}
				
				k++;
				
				
						
			}
			dist=A_dist;
	return A;

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
	 int XE=0; //--Number with elapsed time
	 int XD=0; //--Number with suppl. dijkstra
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
		 
		 while (more) {
			std::clock_t start_time = std::clock();
			bool el=false; //--Elapsed?
			 bool xd=false; //--Extra Dijkstra
			int  j=random_path[2*ii+1];  
			 int path_type=-1;
			 float real_dist_g2=0.0f;
 			//--Distance to j from i (directed or not)
			 vector< vector<int> > path=get_path(i,j,prev_g2);
			 vector<int> gp=good_path2(path,i,j); 
			 real_dist_g2=dist_g2[j];
			 //cout<<i<<" "<<j<<" "<<dist_g1[j]<<" "<<dist_g2[j]<<endl;
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
				} else if (param.heuristic==1||param.heuristic==4) {
					
					// Main big loop
					// test for k accessible nodes
					vector < map<int,int> > prev_g2_w_j(0); //--Reserve 10 ;
					vector < map<int,int> > prev_g2_w_i(0); //--Reserve 10 ;
					
					vector<float> dist_g2_w_j(total_n);
					vector<float> dist_g2_w_i(total_n);
					//--Note: for directed network, we need to test k to j and not j to k
					//--That mean that we cannot order the k to j part, unless we 
					//--For ordering the k in directed network, we use Floyd-Warshall 
					//--To be sure of the distance for i -> k and k -> j
					if (!param.directed) dijkstra(j,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_i, prev_g2_w_i,i);
					
					dijkstra(i,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_j, prev_g2_w_j,j);
					 std::vector<int> k_nodes;     //--K nodes (that we need to include in path)
					
					// Ttest for distance and order the k by this distance
						multimap<float,int> k_nodes_dist;
						for (int l=0; l<total_n;l++) {
							if (Nsb[l]) {
								float mdist=0.0;
								float sdist=0.0;
								if (!param.directed) {
									mdist=std::max(dist_g2_w_i[l],dist_g2_w_j[l]);
									sdist=dist_g2_w_i[l]+dist_g2_w_j[l];
								} else {
									sdist=mdist=dist_g2_w_j[l];
								}
								
								if (!param.directed&&dist_g2_w_i[l]<Inf&&dist_g2_w_j[l]<Inf&&mdist<=param.maxdistance) k_nodes_dist.insert(multimap<float,int>::value_type((param.strategy==1?mdist:sdist),l));
								if (param.directed&&dist_g2_w_j[l]<Inf&&mdist<=param.maxdistance) k_nodes_dist.insert(multimap<float,int>::value_type((param.strategy==1?mdist:sdist),l));
							}
						}
						
						int knodes_size=0;
						float last_distance=0.0f;					
						for (std::multimap<float,int>::iterator it=k_nodes_dist.begin(); it!=k_nodes_dist.end(); ++it)
						{															
							if (knodes_size<param.maxnode) {								
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
									int k=k_nodes[l];
									int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));									
									if (param.heuristic!=4&&(elapsed>param.maxtime)) {
										if (elapsed>param.maxtime) el=true;
										break;
									}
									
									
									//cout<<k<<" [k_nodes.size():"<<k_nodes.size()<<"]\n";
									vector< vector<int> > path_ik=get_path(i,k,prev_g2_w_j);
									//cout<<path_ik.size();
									vector< vector<int> > path_jk;
									//--Note, if directed, we would need to check the path from all k to all j
									if (param.directed) {
										vector < map<int,int> > prev_g2_k(0); //--Reserve 10 ;					
										vector<float> dist_g2_k(total_n);
										dijkstra(k,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_k, prev_g2_k,i);
										path_jk=get_path(k,j,prev_g2_k);
									} else {
										path_jk=get_path(j,k,prev_g2_w_i);																											
									}
									// cout<<node_name[i]<<"-"<<node_name[k]<<"|\n";cout_paths(path_ik);
									 // cout<<node_name[j]<<"-"<<node_name[k]<<"|\n";cout_paths(path_jk);
									for (int pi=0; pi<path_ik.size();pi++) {
										for (int pj=0;pj<path_jk.size();pj++) {											 
											 vector<int> tmp;
											 if (!param.directed) {
												 tmp=join(path_ik[pi],path_jk[pj]);
											 } else {
												 tmp=join2(path_ik[pi],path_jk[pj]);
											 }
											 if (good_path(tmp,i,j)) {
												 if (gp.size()==0||gp.size()>tmp.size()) {
													 gp=tmp; 
												 }
												 nodeadend=true;												
											 } 
										}	
									}									
									//--2. Try again with opt. but costly when we remove nodes duplicated in path i to k and j to k
									if (!nodeadend||param.heuristic==4){
										//gp.clear();
										for (int pi=0; pi<path_ik.size();pi++) {
												int elapsed2 = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
												if (elapsed2>param.maxtime&&param.heuristic!=4) {																										
													el=true;
													break;
												}
												for (int pj=0;pj<path_jk.size();pj++) {											 
													 map<int,int> bad_nodes=intersect(path_ik[pi], path_jk[pj],k,i);
													 // cout<<"** "<<bad_nodes.size()<<" "<<endl;
													 // cout_map(bad_nodes);
													 if (bad_nodes.size()>0) {														
														 vector < map<int,int> > prev_g2_w_i2(total_n); //--Reserve 10 ;	
														 vector<float> dist_g2_w_i2(total_n);														 
														 vector< vector<int> > path_jk2;
														 if (!param.directed) {
															 dijkstra(j,k,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_i2, prev_g2_w_i2,bad_nodes);														 
															 path_jk2=get_path(j,k,prev_g2_w_i2);
														 } else {
															dijkstra(k,j,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_i2, prev_g2_w_i2,bad_nodes);														 
															path_jk2=get_path(k,j,prev_g2_w_i2);
														 }
														 for (int pk=0; pk<path_jk2.size();pk++) {
															 vector<int> tmp2;
															 if (!param.directed) {
																tmp2=join(path_ik[pi],path_jk2[pk]);													 
															 } else {
																 tmp2=join2(path_ik[pi],path_jk2[pk]);													 
															 }
															 if (good_path(tmp2,i,j)) {
																if (gp.size()==0||gp.size()>tmp2.size()) {
																	gp=tmp2;																	
																}
																 nodeadend=true;
																 if (param.heuristic!=4) break;
															 } 
														 } //--End for 	pk	
														 elapsed2 = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
														 if (elapsed2>param.maxtime&&param.heuristic!=4) {																										
															el=true;
															break;
														}
															vector < map<int,int> > prev_g2_w_j2(total_n); //--Reserve 10 ;	
															 vector<float> dist_g2_w_j2(total_n);														 
															 dijkstra(i,k,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_w_j2, prev_g2_w_j2,bad_nodes);
															 vector< vector<int> > path_ik2=get_path(j,k,prev_g2_w_j2);													
															  //cout<<node_name[j]<<"-"<<node_name[k]<<"|\n";cout_paths(path_jk2);
															 for (int pk=0; pk<path_ik2.size();pk++) {
																 vector<int> tmp2;
																 if (!param.directed) {
																	tmp2=join(path_ik2[pk],path_jk[pj]);													 
																 } else {
																	 tmp2=join2(path_ik2[pk],path_jk[pj]);													 
																 }
																 if (good_path(tmp2,i,j)) {
																	if (gp.size()==0||gp.size()>tmp2.size()) {
																		gp=tmp2;																	
																	}
																	 nodeadend=true;
																	 if (param.heuristic!=4) break;
																 } 
															 } //--End for 	pk	
														 
													 } //--Endif badnodes
													if(nodeadend&&param.heuristic!=4) break;													
												} //--End for pj	
												if(nodeadend&&param.heuristic!=4) break;
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
									//cout<<gp.size()<<endl;
								}
								
											
					}
				} else if (param.heuristic==2) {					
					//BRIDES_Y
					vector<float> dist_yen(total_n);			
					
					vector< vector<int> > yen_paths=YenKSP_Stop(i,j,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_yen);					
					bool found=false;
					for (int pi=0; pi<yen_paths.size();pi++) {
						real_dist_g2=dist_yen[pi];
						if (good_path(yen_paths[pi],i,j)) {
							gp=yen_paths[pi];
							found=true;
							break;
						}
					}
					if (found&&real_dist_g2==dist_g1[j]) {
						path_type=c_equal;
					} else if (found) {
						path_type=c_detour;
					} else {
						path_type=c_roadblock;
					}
				} else if (param.heuristic==3) {
					//BRIDES_YC
					   std::vector<int> k_nodes;     //--K nodes (that we need to include in path)
					   map<int, vector< vector<int> > >Pik;
					   map<int, vector< vector<int> > >Pjk;
					   // Order the K using Dijkstra to ik and jk 
						vector < map<int,int> > prev_g2_j(0); //--Reserve 10 ;
						vector<float> dist_g2_j(total_n);					
						dijkstra(j,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_g2_j, prev_g2_j);
						multimap<float,int> k_nodes_dist;
						for (int l=0; l<total_n;l++) {
							if (Nsb[l]) {
								int k=l;
								float ik_dist=dist_g2[k];
								float jk_dist=dist_g2_j[k];
								float mdist=std::max(ik_dist,jk_dist);
								float sdist=ik_dist+jk_dist;
								if (ik_dist<Inf&&jk_dist<Inf&&mdist<=param.maxdistance) {
									k_nodes_dist.insert(multimap<float,int>::value_type((param.strategy==1?mdist:sdist),l));
								}
							}
						}
						int knodes_size=0;
						float last_distance=0;
						for (std::multimap<float,int>::iterator it=k_nodes_dist.begin(); it!=k_nodes_dist.end(); ++it)
						{															
							if (knodes_size<param.maxnode) {								
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
								bool nodeadend=false;
								gp.clear();
								for (int l=0; l<k_nodes.size();l++) {
									int k=k_nodes[l];
									int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
									if (param.heuristic!=4&&(elapsed>param.maxtime)) {
										if (elapsed>param.maxtime) el=true;
										break;
									}
									vector<float> dist_yen_i(total_n);
									vector<float> dist_yen_j(total_n);
									vector< vector<int> > Pik=YenKSP(i,k,param.max_individual_path,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_yen_i);
									vector< vector<int> > Pjk=YenKSP(j,k,param.max_individual_path,undirected_adjlist_g2,undirected_adjlist_g2_dist, dist_yen_j);
									for (int pi=0; pi<Pik.size();pi++) {
										for (int pj=0;pj<Pjk.size();pj++) {											 
											 vector<int> tmp=join(Pik[pi],Pjk[pj]);
											 if (good_path(tmp,i,j)) {
												 if (gp.size()==0||gp.size()>tmp.size()) {
													 gp=tmp; 
												 }
												 nodeadend=true;
											 } 
										}	
									}
									if(nodeadend) break;
							   }
							   if (nodeadend) {
									real_dist_g2=dist_path(gp,undirected_adjlist_g2_dist);
									if (real_dist_g2==dist_g1[j]) {
										path_type=c_equal;
									} else {
										path_type=c_detour;
									}
								} else {
									path_type=c_roadblock;
								} 
					}
				} else if (param.heuristic==5) {
					//heuristic DFS
					gp=get_path_dfs_BRIDES(i,j,undirected_adjlist_g2);
					
					if (good_path(gp,i,j)) {
						real_dist_g2=gp.size()-1; //Not the real distance
						if (real_dist_g2==dist_g1[j]) {
								path_type=c_equal;
							} else {
								path_type=c_detour;
							}					
					} else {
						path_type=c_roadblock;
					}
				}
			 } //--End else
			#pragma omp atomic
			XE+=(el?1:0);
			#pragma omp atomic
			XD+=(xd?1:0);
			#pragma omp atomic
			total++;			 		
			int elapsed = ((std::clock() - start_time) / (double) (CLOCKS_PER_SEC / 1000));
			#pragma omp atomic
			total_time+=elapsed;
			add_statistics(i,j,path_type,gp);
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
	results.push_back(XE);
	results.push_back(XD);
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
	  
      if (!param.found_g2&&!param.found_attributes&&!param.use_time_code) {
          cout<<"Error. Unable to locate the g2 or attributes file g2:"<<param.graph2<<" attributes:"<<param.attributes<<endl;
		  destructor();
          exit(-1);
      }
	  if (param.use_time_code) {
		  cout<<"Saving to "<<param.outputfile<<" the nodes from :"<<param.start_time_code<<" to "<<param.end_time_code<<endl;
		  FileOutput.open(param.outputfile);	  
		   if(!FileOutput.is_open()) {
			  printf("\n%s: verbose result file open failed...",param.outputfile);
			  destructor();
			  exit(1);
			   
		   }
		   output_graph();
		   
		   FileOutput.close();
		  exit(0);
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
		debug_info(); //--for debug
		
		vector<float> A;
		
		for (int i=0; i<total_n;i++) {
			for (int j=0;j<total_n;j++) {
				if (i<j) {
					cout<<"["<<node_name[i]<<"-"<<node_name[j]<<"]"<<endl;
					//vector< vector<int> > pa=YenKSP(i,j,undirected_adjlist_g2,undirected_adjlist_g2_dist, A);		
					vector<int> pa=get_path_dfs_BRIDES(i,j,undirected_adjlist_g2);
					//for (int p=0; p<pa.size();p++) {
						cout_paths(pa);
					//}
				}
			}
		}
		
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
	int XE=0;
	int XD=0;
	cout<<"\n============================ PARTIAL RESULTS ==================================\n";	
	cout<<" Group\tB\tR\tI\tD\tE\tS\tTotal\tCPU time (ms)"<<endl;
	if (param.verbose) {
		  FileOutput<<"\n============================ PARTIAL RESULTS ==================================\n";	
		  FileOutput<<"src\tdest\tDist_X\tDist_Y\tBRIDES\tCPU time (ms)\tPath\tTaxa"<<endl;
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
		XE+=r[8];
		XD+=r[9];
		#pragma omp critical 
		{		
			cout<<" "<<grp<<"\t"<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\t"<<r[4]<<"\t"<<r[5]<<"\t"<<r[6]<<"\t"<<r[7]<<"\t"<<endl;		
		}
	}
	#pragma omp critical 
	{		
		time_t endtime = time(0);
		double ttime=difftime(endtime, starttime);
		cout<<"\n================================= INFO =======================================\n";
		cout<<" (B) Breakthrough : path impossible in network X but possible in network Y\n";
		cout<<" (R) Roadblock    : path possible in network X but impossible in network Y\n";
		cout<<" (I) Impasse      : path impossible in both networks, X and Y\n";
		cout<<" (D) Detour       : path shorter in network X than in network Y\n";
		cout<<" (E) Equal        : path of same length in networks X and Y\n";
		cout<<" (S) Shortcut     : path longer in network X than in network Y\n";
		// cout<<"\n============================== STATISTICS ====================================\n";
		// cout<<"NK\tB\tR\tI\tD\tE\tS\tInside\tAttribute"<<endl;
		// for (int i=0; i<total_n;i++) {
				// if (Nsu[i]&&!Nsb[i]) {
					 // cout<<node_name[i]<<"\t"<<statistics[i][0]<<"\t"<<statistics[i][1]<<"\t"<<statistics[i][2]<<"\t"<<statistics[i][3]<<"\t"<<statistics[i][4]<<"\t"<<statistics[i][5]<<"\t"<<statistics[i][6]<<"\t"<<attributes[i]<<endl;
				// }
			// }
			// cout<<endl<<"K\tInside\tAttribute"<<endl;
			
			// for (int i=0; i<total_n;i++) {
				// if (Nsb[i]) {
					// cout<<node_name[i]<<"\t"<<statistics[i][6]<<"\t"<<attributes[i]<<endl;
				// }
			// }
		cout<<"\n================================ RESULTS ======================================\n";
		cout<<"\tB\tR\tI\tD\tE\tS\tTotal\tTime (s)"<<endl;	
		cout<<"\t"<<B<<"\t"<<R<<"\t"<<I<<"\t"<<D<<"\t"<<E<<"\t"<<S<<"\t"<<total<<"\t"<<ttime<<endl;	   
		cout<<"===============================================================================\n";	
		if (param.verbose) {
		    FileOutput<<"\n================================= INFO =======================================\n";
			FileOutput<<" (B) Breakthrough : path impossible in network X but possible in network Y\n";
			FileOutput<<" (R) Roadblock    : path possible in network X but impossible in network Y\n";
			FileOutput<<" (I) Impasse      : path impossible in both networks, X and Y\n";
			FileOutput<<" (D) Detour       : path shorter in network X than in network Y\n";
			FileOutput<<" (E) Equal        : path of same length in networks X and Y\n";
			FileOutput<<" (S) Shortcut     : path longer in network X than in network Y\n";
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
			FileOutput<<"*Inside: number of time a node belongs to another path;\n K: added nodes;NK: original nodes\n";
			FileOutput<<"\n================================ RESULTS ======================================\n";
			if (XE>0) FileOutput<<"Notice: "<<XE<<" pathway(s) classification took more than "<<(param.maxtime/1000)<<" seconds (maxtime).\n        Consider increasing the maxtime parameter."<<endl;
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
		//DO we allow colored nod in networkX?
		//if (e.to!=-1&&!Nsb[e.to]&&!Nsb[e.from]&&Nsu[e.to]&&Nsu[e.from])
		
		if (e.to!=-1&&Nsu[e.to]&&Nsu[e.from]) {
				undirected_adjlist_g1[e.from].push_back(e.to);
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
        while(getline(file, line))
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
				} 				
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
        float dist;
	float total_line=0;
   try {
        while(getline(file, line))
        {
            std::stringstream   linestream(line);
            // std::getline(linestream, data, '\t');
            std::vector<std::string> tokens;
            std::string token;
			//line.erase(line.find_last_not_of(" \n\r\t")+1);
            if (line.find("#")!=0&&line.find("%")!=0&&!line.empty()) {
                total_line++;
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
                if (tokens.size()>2) dist=atof(tokens[2].c_str());
				if (dist<=0.0f) dist=1.0f; //--Don't permit negative distance
                 if (!param.use_time_code) {
					 edge t;             
					 t.from=add_node(from);
					 t.to=add_node(to);
					 t.dist=dist; 
					 local_node.insert ( std::pair<int,int>(t.from,0));
					 if (t.to!=-1) local_node.insert ( std::pair<int,int>(t.to,0));				 
					 if (t.to!=t.from)  edge_list.push_back(t);
				 } else {
					 if (dist>=param.start_time_code&&dist<=param.end_time_code) {
						 edge t;             
						 t.from=add_node(from);
						 t.to=add_node(to);
						 t.dist=dist; 
						 local_node.insert ( std::pair<int,int>(t.from,0));
						 if (t.to!=-1) local_node.insert ( std::pair<int,int>(t.to,0));				 
						 if (t.to!=t.from)  edge_list.push_back(t);
					 }
				 }
            }
        }
   } catch(const std::exception& e) {}   
   cout<<"Done loading "<<filename<<"... ("<<local_node.size()<<" nodes, "<<total_line<<" edges)"<<endl;
   file.close();
  return(local_node);
} 

int randomInt (int min, int max){
    int n = rand()%(max-min)+min; 
    return n;
}

/**
 * Save to output the selected network
 */
void output_graph() {
	for (int i=0; i<edge_list_g1.size();i++) {
			if (edge_list_g1[i].to!=-1) {
				FileOutput<<node_name[edge_list_g1[i].from]<<"\t"<<node_name[edge_list_g1[i].to]<<"\t"<<edge_list_g1[i].dist<<endl;
			} else {
				FileOutput<<node_name[edge_list_g1[i].from]<<"\t"<<edge_list_g1[i].dist<<endl;
			}
	}
	
} 

/**
 * output some information to the file
 * Format will be similar to Nexus
 */
void output_header(char** argv) {
 
    //cout<<startMessage; 
	
	cout<<"\n=============================== PARAMETERS ====================================\n";		
	
	//-- OpenMP
	  cout<<"\n -=[Multicores support]=-"<<endl;
	  
	  cout<<" Maximum threads               : "<<omp_get_max_threads()<<endl;
	  cout<<" Number of cores               : "<<omp_get_num_procs()<<endl;
	 //--Loaded network 
	  cout<<"\n -=[Input files]=-"<<endl;
	                       					  
	  if (param.directed) cout<<" Networks                      : directed"<<endl;      
	  if (!param.directed) cout<<" Networks                      : undirected"<<endl;      
	  
	  cout<<" Network X                     : "<<param.graph1<<endl;
						  					  
	  if (param.found_g2) cout<<" Network Y                     : "<<param.graph2<<endl;
							      
	 
	  cout<<" Nodes in network X            : "<<total_nonk_node<<endl;   
      cout<<" Nodes in network Y            : "<<total_n_g2<<endl; 
	  //cout<<"Edges in networkX: "<<edge_list_g1.size()<<endl;   
      //cout<<"Edges in networkY: "<<edge_list_g2.size()<<endl; 
	  
      if (param.found_attributes)  {
		 
		  cout<<" Node's attributes|count       : "<<uniques_attributes.size()<<endl;   
		  for(map<string,int>::iterator it = uniques_attributes.begin(); it != uniques_attributes.end(); ++it) {
			  cout<<"\t"<<it->first<<"|"<<it->second<<endl;
		  }
	  }	   
    	                       
	  if (!param.nonK.empty()) cout<<" Original nodes (NK) attributes: "<<param.nonK<<endl;
	  cout<<" Total of added nodes          : "<<total_k_node<<endl; 
	  if (!param.K.empty()) cout<<" Added nodes (K) attributes    : "<<param.K<<endl;
	  
	  
	  cout<<" Total paths                   : "<<total_paths<<endl; 
      cout<<"\n -=[Run parameters]=-"<<endl;
	  if (param.random!=-1) {		
	      
		cout<<" Running mode                  : random"<<endl;  
		if (param.random<1.0) {
			cout<<" Number of investigated paths  : "<<(param.random*100)<<"% ("<<(int)(total_paths*param.random)<<")"<<endl; 
		} else {
			cout<<" Number of investigated paths  : "<<param.random<<endl;  
		}
	  } else {
	      cout<<" Running mode                  : normal"<<endl;
	  }
	  if (param.use_dist&&param.inv_dist) {
		cout<<" Edge distance                 : inverse weights"<<endl;      	  
	  } else if (param.use_dist) {
		cout<<" Edge distance                 : edge weights"<<endl;    
	  } else {
		cout<<" Edge distance                 : unweighted"<<endl;     
	  }
	  
	  cout<<" Group size                    : "<<param.size<<endl;      
	  cout<<" Number of groups              : "<<total_group<<endl;      
	  cout<<" First group                   : "<<param.first<<endl;      
      cout<<" Last group                    : "<<param.last<<endl;
      cout<<" Maxdistanc                    : "<<param.maxdistance<<endl;       
      cout<<" Maxnode                       : "<<param.maxnode<<endl;       	  
	  cout<<" Maxtime (s)                   : "<<(param.maxtime/1000)<<endl;    
	  cout<<" Maxpathnumber                 : "<<param.max_individual_path<<endl; 	                                             
	  cout<<"\n -=[Miscellaneous]=-"<<endl;
 	  if (param.seed<0) {
		 
		  cout<<" Seed                          : clock time"<<endl;      
	  } else {
		  cout<<" Seed                          : "<<param.seed<<endl;      
	  }  
	      cout<<" Added nodes ordering          : strategy "<<param.strategy<<endl;      
	  switch(param.heuristic) {
		  case 1:cout<<" Algorithm                     : BRIDES     (1)"<<endl;break;
		  case 2:cout<<" Algorithm                     : BRIDES_Y   (2)"<<endl;break;
		  case 3:cout<<" Algorithm                     : BRIDES_YC  (3)"<<endl;break;
		  case 4:cout<<" Algorithm                     : BRIDES_EC  (4)"<<endl;break;
		  case 5:cout<<" Algorithm                     : DFS        (5)"<<endl;break;
		  case 6:cout<<" Algorithm                     : B. RELAXED (6)"<<endl;break;
	  }
	  
	  
	  if (param.maxthread!=0) {
		  cout<<" Maxthread                     : "<<param.maxthread<<endl;      
	  }
	  if (param.verbose) 
		  cout<<" Output file                   : "<<param.outputfile<<endl;      
	  cout<<"\n===============================================================================\n";	
      
	  if (param.verbose) {
		  FileOutput<<"\n=============================== PARAMETERS ====================================\n";		
	
			//-- OpenMP
			  FileOutput<<"\n -=[Multicores support]=-"<<endl;			 
			  FileOutput<<" Maximum threads               : "<<omp_get_max_threads()<<endl;
			  FileOutput<<" Number of cores               : "<<omp_get_num_procs()<<endl;
			 //--Loaded network 
			  FileOutput<<"\n -=[Input files]=-"<<endl;
			  if (param.directed)  FileOutput<<" Networks                      : directed"<<endl;      
			  if (!param.directed) FileOutput<<" Networks                      : undirected"<<endl;      
			  	    
			  FileOutput<<" Network X                     : "<<param.graph1<<endl;
			                            
			  if (param.found_g2) FileOutput<<" Network Y                     : "<<param.graph2<<endl;
			           
			  FileOutput<<" Nodes in network X            : "<<total_nonk_node<<endl;   
			  FileOutput<<" Nodes in network Y            : "<<total_n_g2<<endl; 			  
			  if (param.found_attributes)  {
				        
				  FileOutput<<" Node's attributes|count       : "<<uniques_attributes.size()<<endl;   
				  for(map<string,int>::iterator it = uniques_attributes.begin(); it != uniques_attributes.end(); ++it) {
					  FileOutput<<"\t"<<it->first<<"|"<<it->second<<endl;
				  }
			  }	   
											 
			  if (!param.nonK.empty()) FileOutput<<" Non-added nodes(NK) attributes: "<<param.nonK<<endl;
         			
			  FileOutput<<" Total of added nodes          : "<<total_k_node<<endl; 
			  if (!param.K.empty()) FileOutput<<" Added nodes (K) attributes    : "<<param.K<<endl;
			  FileOutput<<" Total paths                   : "<<total_paths<<endl; 
      
			  //--Running parameters
			  FileOutput<<"\n -=[Run parameters]=-"<<endl;
			  if (param.random!=-1) {		
			    
				FileOutput<<" Running mode                  : random"<<endl;  
				if (param.random<1.0) {
					FileOutput<<" Number of investigated paths  : "<<(param.random*100)<<"% ("<<(int)(total_paths*param.random)<<")"<<endl; 
				} else {
					FileOutput<<" Number of investigated paths  : "<<param.random<<endl;  
				}
			  } else { 
			   
			    FileOutput<<" Running mode                  : normal"<<endl;
			  }
			  if (param.use_dist&&param.inv_dist) {
				FileOutput<<" Edge distance                 : inverse weights"<<endl;      	  
			  } else if (param.use_dist) {
				FileOutput<<" Edge distance                 : edge weights"<<endl;    
			  } else {
				FileOutput<<" Edge distance                 : unweighted"<<endl;     
			  }
			  FileOutput<<" Group size                    : "<<param.size<<endl;      
			  FileOutput<<" Number of groups              : "<<total_group<<endl;      
			  FileOutput<<" First group                   : "<<param.first<<endl;      
			  FileOutput<<" Last group                    : "<<param.last<<endl;
			  FileOutput<<" Maxdistance                   : "<<param.maxdistance<<endl;       
			  FileOutput<<" Maxnode                       : "<<param.maxnode<<endl;       	  
			  FileOutput<<" Maxtime (s)                   : "<<(param.maxtime/1000)<<endl;      
			  FileOutput<<" Maxpathnumber                 : "<<param.max_individual_path<<endl; 	                                             
			  
			  FileOutput<<"\n -=[Miscellaneous]=-"<<endl;
			  if (param.seed<0) {
				  FileOutput<<" Seed                          : clock time"<<endl;				  
			  } else {
				  FileOutput<<" Seed                          : "<<param.seed<<endl;      
			  }  
			       
			FileOutput<<" Added nodes ordering          : strategy "<<param.strategy<<endl;      
			  switch(param.heuristic) {
				  case 1:FileOutput<<" Algorithm                     : BRIDES     (1)"<<endl;break;
				  case 2:FileOutput<<" Algorithm                     : BRIDES_Y   (2)"<<endl;break;
				  case 3:FileOutput<<" Algorithm                     : BRIDES_YC  (3)"<<endl;break;
				  case 4:FileOutput<<" Algorithm                     : BRIDES_EC  (4)"<<endl;break;
				  case 5:FileOutput<<" Algorithm                     : DFS        (5)"<<endl;break;
				  case 6:FileOutput<<" Algorithm                     : B. RELAX   (6)"<<endl;break;
			  }
			  if (param.maxthread!=0) {
				  FileOutput<<" Maxthread                     : "<<param.maxthread<<endl;      
			  }
			  if (param.verbose)      
				  FileOutput<<" Output file                   : "<<param.outputfile<<endl;  
			  FileOutput<<"\n===============================================================================\n";	
	  }
			
	
}

void help(){
	printf("\n%s",description);	
	printf("\nUsage :\nBRIDES -X=[filename] -Y=[filename] -outfile=[filename]\n");
	printf("       -maxdistance=[1..n-1] -maxnodes=[1..n-1]\n");
         printf("\nParameters :");
        printf("\n-X=file           [filename for original network X]");
        printf("\n-Y=file           [filename for augmented network Y]");
		printf("\n-A=file           [filename for node attributes]");
		printf("\n-usedist          [use edge distances found in network files]");
		printf("\n-invdist          [use the inverse of edge distances]");
		printf("\n-directed         [specifies that the networks are directed]");
		printf("\n-K=B,C            [attributes to consider as added nodes (K), e.g. B,C]");
        printf("\n-NK=A             [attributes to considers as original nodes (NK), e.g. A]");
		printf("\n-random=XX        [sample XX random paths]");
		printf("\n-first=1          [first group of path to process]");
        printf("\n-last=n           [last group of path to process]");
		printf("\n-size=1000        [group size, default: 1000]");      
        printf("\n-maxdistance=100  [maximum distance to an added node (K), default: 100]");
		printf("\n-maxnodes=100     [maximum number of added nodes to examine, default: 100]");
		printf("\n-maxtime=10       [maximum time for each path search, default: 10 second]");
		printf("\n-maxpathnumber=100[maximum number of shortest-path return by Dijkstra: 100]");
		printf("\n-maxthread=XXX    [maximum number of OpenMP threads to use, default unlimited]");
        printf("\n-outfile=file     [output path information into file: taxa, distance, etc.]");
		printf("\n-seed=999         [set the random seed generator to a specific seed]");
		printf("\n-strategy=1       [set the K nodes ordering strategy: (1) maxdist.,(2) sum.]");
		printf("\n-algo=1           [1-BRIDES, 2-BRIDES_Y, 3-BRIDES_YC, 4-BRIDES_EC,5-DFS]");
		printf("\n                  [6-RELAXED]");
		printf("\n\nExample : \n: ./BRIDES -X=networkX.txt -Y=networkY.txt\n\n");
}

//--Get each parameters
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
		(*param).maxtime=100000; //--Maxtime in millsecond (100s per path)
		(*param).size=1000; //-- default group size
		(*param).max_individual_path=100; //Since we keep all shortest path, this is the number of shortest path investigated for each path
		(*param).maxthread=0; 
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
         (*param).heuristic=1;
		 (*param).strategy=1;
		 (*param).info=false;
		 (*param).directed=false;         
		 (*param).verbose=false;         
		 (*param).found_g1=false;
         (*param).found_g2=false;
         (*param).found_attributes=false;
		 (*param).use_dist=false;
		 (*param).inv_dist=false;
		 (*param).removeK_from_X=false;
		 (*param).output_graph=false;
		 (*param).use_time_code=false;
		 (*param).start_time_code=0;   //Time code for including vertex (Unix timestamps)
		 (*param).end_time_code=1;
		
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
			else if(strcmp("verbose",champs) == 0||strcmp("output",champs) == 0||strcmp("outfile",champs) == 0){
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
			//======= start timecode ==============       
                         else if(strcmp("tstart",champs) == 0){                                
                              (*param).start_time_code = atoi(contenu);
							  (*param).use_time_code=true;
			}
            //======= end ==============       
                         else if(strcmp("tend",champs) == 0){                                
                              (*param).end_time_code = atoi(contenu);
							  (*param).use_time_code=true;
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
            else if(strcmp("maxpathnumber",champs) == 0){                                
                              (*param).max_individual_path = atoi(contenu);
			}
			//======== input file ==============
			else if(strcmp("g1",champs) == 0||strcmp("X",champs) == 0||strcmp("x",champs) == 0){
                                sprintf((*param).graph1,"%s",contenu);	
                                FILE *Input1;
                                if ((Input1 = fopen((*param).graph1,"r"))!=0) {(*param).found_g1=true;fclose(Input1);}
			}
            else if(strcmp("g2",champs) == 0||strcmp("Y",champs) == 0||strcmp("y",champs) == 0){
                                sprintf((*param).graph2,"%s",contenu);	
                               // (*param).use_g2=true;
                                 FILE *Input1;
                                if ((Input1 = fopen((*param).graph2,"r"))!=0) {(*param).found_g2=true;fclose(Input1);}
             }
			 else if(strcmp("attributes",champs) == 0||strcmp("attr",champs) == 0||strcmp("A",champs) == 0){
                                sprintf((*param).attributes,"%s",contenu);	
                               // (*param).use_g2=true;
                                 FILE *Input1;
                                if ((Input1 = fopen((*param).attributes,"r"))!=0) {(*param).found_attributes=true;fclose(Input1);}
			}
           //======== heuristic ================
			else if(strcmp("heuristic",champs) == 0||strcmp("algo",champs) == 0){
                  (*param).heuristic = atoi(contenu);
			}
			//============ output file ===============
			// else if(strcmp("outputfile",champs) == 0){
				// sprintf((*param).outputfile,"%s",contenu);	
                                // //printf("%s\n",(*param).outputfile );
			// }
			else if(strcmp("removeK_from_X",champs) == 0){
					(*param).removeK_from_X = true;
			}
			//============= version ===============
			else if(strcmp("directed",champs) == 0){
					(*param).directed = true;
			}
			 //======== heuristic ================
			else if(strcmp("strategy",champs) == 0){
                  (*param).strategy = atoi(contenu);
			}
			//============= info   ===============
			else if(strcmp("info",champs) == 0||strcmp("debug",champs) == 0){
					(*param).info = true;
			}
			//============= thread   ===============
			else if(strcmp("maxthread",champs) == 0){
					(*param).maxthread =  atoi(contenu);
					omp_set_dynamic(0);     
					omp_set_num_threads((*param).maxthread); 
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
