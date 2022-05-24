#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>

#include "slurm/slurm.h"

#include "src/common/slurm_topology.h"
#include "src/common/switch.h"
#include "src/common/node_conf.h"
#include "src/common/xmalloc.h"
#include "src/common/slurm_step_layout.c"


#include "src/slurmctld/job_scheduler.h"
#include "src/slurmctld/proc_req.h"
#include "src/slurmctld/slurmctld.h"
#include "src/slurmctld/calc_hops.h"
#include "src/plugins/select/linear/select_linear.h"

#include <metis.h>

extern int switch_levels;
extern struct switch_record *switch_record_table;
extern int switch_record_cnt;
uint32_t* node_cnt;

#define MAX_PARTITION_LEVEL 5

#define SWITCH_ORDER_SIZE 100
extern struct table *alloc_node_table;
extern struct table *switch_idx_table;

extern int* switch2node;
//define a structure to maintain (switchindex, freenodes(switchindex)) pair
typedef struct {
		int switch_idx; 
		uint32_t free_nodes;
} node_cnt_struct;
node_cnt_struct* node_cnt_array;

typedef struct {
		int no;
		int size;
} part_struct;


uint32_t hashCode(struct table *t,uint32_t key){
    if(key<0)
        return -(key%t->size);
    return key%t->size;
}

void insert(struct table *t,uint32_t key,int* val, int n){
    uint32_t pos = hashCode(t,key);
    struct node *list = t->list[pos];
    struct node *newNode = (struct node*)malloc(sizeof(struct node));
    newNode->key = key;
    for(int i=0;i<n;i++)
        newNode->val[i] = val[i];
    newNode->n = n;
    newNode->next = list;
    t->list[pos] = newNode;
}

// Arr is the array where the lookup result is stored if key is found
int lookup(struct table *t,uint32_t key, int*arr, int *n){ 
    uint32_t pos = hashCode(t,key);
    //debug("JobId=%d Pos=%d",key,pos);
    struct node *list = t->list[pos];
    struct node *temp = list;
    while(temp){
        if(temp->key==key){
            *n = temp->n;
            for(int i=0;i<*n;i++)
                arr[i] = temp->val[i];
            return 0;
        }
        temp = temp->next;
    }
    return -1;
}

// To sort array in descending order
int desc_cmp(const void *a, const void *b) {
        int idxa = *(int *)a;
        int idxb = *(int *)b; 
        if (node_cnt[idxa] != node_cnt[idxb])
                return node_cnt[idxa] > node_cnt[idxb] ? -1 : 1;
        else
                return switch_record_table[idxa].comm_jobs < switch_record_table[idxb].comm_jobs ? -1 : 1;
}

// To sort array in increasing order 
int inc_cmp(const void *a, const void *b) {
	int idxa = *(int *)a;
	int idxb = *(int *)b;
	if (node_cnt[idxa] == 0)
		return 1;
	else if (node_cnt[idxb] == 0)
		return -1;
	else
		if(node_cnt[idxa] != node_cnt[idxb])
			return node_cnt[idxa] < node_cnt[idxb] ? -1 : 1;
		else
			return switch_record_table[idxa].comm_jobs > switch_record_table[idxb].comm_jobs ? -1 : 1;
}

// To sort array in descending order
int desc_cmp_struct(const void *a, const void *b) {
		node_cnt_struct val1 = *(node_cnt_struct*)a;
		node_cnt_struct val2 = *(node_cnt_struct*)b;

		if (val1.free_nodes != val2.free_nodes)
				return val1.free_nodes > val2.free_nodes ? -1 : 1;
		else
				return switch_record_table[val1.switch_idx].comm_jobs < switch_record_table[val2.switch_idx].comm_jobs ? -1 : 1;
}

// To sort array in increasing order 
int inc_cmp_struct(const void *a, const void *b) {  
		node_cnt_struct val1 = *(node_cnt_struct*)a;
		node_cnt_struct val2 = *(node_cnt_struct*)b;
		
		if (val1.free_nodes == 0)
				return 1;
		else if (val2.free_nodes == 0)
				return -1;
		else
				if (val1.free_nodes != val2.free_nodes)
						return val1.free_nodes < val2.free_nodes ? -1 : 1;
				else
						return switch_record_table[val1.switch_idx].comm_jobs > switch_record_table[val2.switch_idx].comm_jobs ? -1 : 1;
}

int comp(const int * a, const int * b) {
		return *b - *a;
}

int dec_cmp_part(const void *a, const void *b) {
		part_struct val1 = *(part_struct*)a;
		part_struct val2 = *(part_struct*)b;
		
		return val2.size - val1.size;
}

//calculated hop-bytes based on switches array (which gives which node is present on which switch) for the entire matrix
long double calc_hops_comm_matrix(int switches[], int size, struct job_record *job_ptr) {
        //file pointer to read communication pattern from communication file path
        FILE* fcomm = fopen(job_ptr->comment + 2, "r");         
        int nodes = job_ptr->node_cnt; //number of nodes required for job
        
        //scanning communication pattern from file (path provided in job comment parameter) and storing in commpattern 2D array
        int node1, node2;
        long double val;
        long double comm_hops_max = 0; //stores max hops encountered 
        long double comm_hops_local;
        long double comm_hops_total = 0; //total hop bytes

        for (node1 = 0; node1 < nodes; node1++) {
                for (node2 = 0; node2 < nodes; node2++) {
                        fscanf(fcomm, "%Lf", &val);
						//if mat[node1][node2] > 0
                        if (val > 0) {
                                float c = 0, c1 = 0, c2 = 0, c3 = 0;
                                if (switches[node1] == switches[node2]) {
										//find contention factor
                                        c = (switch_record_table[switches[node1]].comm_jobs) / ((float)switch_record_table[switches[node1]].num_nodes);
                                        comm_hops_local = 2 * val * (1 + c); 
                                }
                                else {
										//find contention on both switches
                                        c1 = (switch_record_table[switches[node1]].comm_jobs) / ((float)switch_record_table[switches[node1]].num_nodes);
                                        c2 = (switch_record_table[switches[node2]].comm_jobs) / ((float)switch_record_table[switches[node2]].num_nodes);
                                        c3 = (switch_record_table[switches[node1]].comm_jobs + switch_record_table[switches[node2]].comm_jobs) / 
                                                ((float)switch_record_table[switches[node1]].num_nodes + (float)switch_record_table[switches[node2]].num_nodes);
                                        c = c1 + c2 + c3/2; //final contention
                                        comm_hops_local = (2 * (switch_levels+1)) * val * (1 + c);
                                }
								//update max hop-bytes
                                if (comm_hops_local > comm_hops_max) {
                                        comm_hops_max = comm_hops_local;
                                }
                                comm_hops_total += comm_hops_local;
                        }
                }
        }
        fclose(fcomm); //close file
        return comm_hops_total / 2; //since c[i][j] and c[j][i] should only be counted once for finding hop bytes
		//return comm_hops_max;
}

//takes as input the xadj, adjncy, adjwgt, vsize parameters as required by METIS
//nodes gives total number of nodes, npart gives total number of required partitions
//totalv is a boolean which denotes whether to use totalv or edgecut (default) as objective function
idx_t* graph_to_partition(idx_t* xadj, idx_t* adjncy, idx_t* adjwgt, int nodes, int npart, int totalv, float part_ufactor) {
        // debug("part array allocating");
        // debug("nodes: %d", nodes);
        idx_t* part = (idx_t*)malloc(nodes * sizeof(idx_t)); //partition array
        // debug("part array allocated");
        idx_t nvtxs = nodes;
        idx_t ncon = 1; //The number of balancing constraints. It should be at least 1.
        idx_t objval; //stores the edge-cut or the total communication volume of the partitioning solution.
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_UFACTOR] = 1000;
        //METIS_OPTION_UFACTOR 1.5 denotes the size of largest partition is 1.5 times the average, the factor in this case should be 500
        float ufactor = (part_ufactor - 1) * 1000; //factor = (1 - reqfactor) * 1000
        if (ufactor == 0) {
                ufactor = 30; //default value is 30. i.e. factor of 1.03
        }
        options[METIS_OPTION_UFACTOR] = ufactor;
        debug("ufactor for partitioning: %f", ufactor);

		//use totalv as objective (reduce the total communication volume)
		//it makes use of the vsize array where vsize[i] denotes the data sent to all other nodes by node i
        if (totalv == 1) {
                options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
				// debug("Calling kway with obj=totalv");
				//giving vsize here gives segmentation fault for 512 nodes (does not give for < 512 nodes)
                // int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, vsize, 
                //         NULL, &npart, NULL, NULL, options, &objval, part);
                int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, 
                        adjwgt, &npart, NULL, NULL, options, &objval, part);
                debug("partitioning ret:%d (totalvol) : %"PRIu32"", ret, objval);
        }
		//use edgecut as objective (default)
		//does not make use of vsize array
        else {
                int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, 
                        adjwgt, &npart, NULL, NULL, options, &objval, part);        
                debug("partitioning ret:%d (edgecut) : %"PRIu32"", ret, objval);
        }

        // debug("graphtopartition complete");
		//part[i] gives the partition for node i
        return part;
}

//takes as input the 2D comm matrix and generates the partitioning based on it
idx_t* matrix_to_partition(int nodes, int edges, uint64_t* commpattern, int npart, float part_ufactor) {
        idx_t* xadj = (idx_t*)malloc((nodes+1) * sizeof(idx_t)); //initialize xadj array
        idx_t* adjncy = (idx_t*)malloc((edges) * sizeof(idx_t)); //initialize adjncy array
        idx_t* adjwgt = (idx_t*)malloc((edges) * sizeof(idx_t)); //initialize adjwgt array
        //idx_t* vsize = (idx_t*)malloc((nodes) * sizeof(idx_t)); //vsize array (to be used when obj = minimize total comm volume)

        xadj[0] = 0; //start for vertex 0
        int k = 0; //to populate adjncy and adjwgt arrays

        for (int i = 0; i < nodes; i++) {
                //vsize[i] = 0;
                int adj_edges = 0; //number of edges adjacent to vertex i
                for (int j = 0; j < nodes; j++) {
                        if (commpattern[i*nodes + j] / 4096 > 0) {
                                adjncy[k] = j; //j is adjacent to i
                                adjwgt[k] = commpattern[i*nodes + j] / 4096; //weight of edge i-j
                                // debug("%"PRIu32"", adjwgt[k]);
                                //vsize[i] += commpattern[i*nodes + j]; //vertex i sends to j
                                k += 1;
                                adj_edges += 1;
                        }
                }
                xadj[i+1] = xadj[i] + adj_edges;
        }

        // for (int i = 0; i < nodes; i++) {
        //         debug("vsize[%d]: %d", i, vsize[i]);
        // }
        debug("adjacency matrix for graph loaded");

        // debug("Nodes: %d, Edges: %d", nodes, edges);
        // for (i = 0; i <= nodes; i++) {
        //         debug("xadj[%d]: %d", i, xadj[i]);
        // }
        // for (i = 0; i < edges; i++) {
        //         debug("adjncy[%d]: %d,  adjwgt[%d]: %d", i, adjncy[i], i, adjwgt[i]);
        // }

        idx_t* part = graph_to_partition(xadj, adjncy, adjwgt, nodes, npart, 0, part_ufactor);
        //idx_t* part2 = graph_to_partition(xadj, adjncy, adjwgt, vsize, nodes, npart, 1);
        //free(part2);

        // for (int i = 0; i < nodes; i++)
        // {
	//         debug ("part[%d]: %d", i, part[i]);
        // }
        free(xadj);
        free(adjncy);
        free(adjwgt);
        //free(vsize);

        return part;
} 

//map nodes from partition pno to switch switchno in the switch2node array
void map_nodes_to_switch(int* switch2node, int curr_alloc_nodes, int* node_alloc, int pno, int switchno, idx_t* part, int nodes, int* groups, int partition_level, int jobnodes) {
        //iterate over all nodes in partition level = partition_level
        for (int i = 0; i < nodes; i++) {
                //if node i is on partition pno
                if (part[i] == pno) {
                        int real_rank = i;
                        int curr_level = partition_level;
                        //find real rank for node i (if partition_level > 0)
                        while (curr_level > 0) {
                                real_rank = groups[(curr_level-1)*jobnodes + real_rank]; //groups[level-1][current_rank]
                                curr_level -= 1; //decrement current level
                        }
                        switch2node[switchno*jobnodes + curr_alloc_nodes] = real_rank; //switch2node[switchno][i] = node
                        curr_alloc_nodes += 1; //increment alloc_nodes on switch
                        node_alloc[real_rank] = 1; //set node to allocated
                }
        }
}

//switch_node_cnt contains number of free nodes in every switch
//want_nodes contains number of nodes required for the job
//switch_alloc_node will contain number of nodes allocated on each switch, after the function has completed
//switch2node will contain the nodes that are to be placed on a particular switch, after the function has completed
//2D array switch2node[#switches][#nodesrequired] is flattened to 1D
void combal_alloc(struct job_record *job_ptr, uint32_t* switch_node_cnt, int* switch_idx, 
                        uint32_t want_nodes, int* switch_alloc_nodes, int* switch2node) {
        uint32_t curr_size = want_nodes; //total number of required nodes
        uint32_t rem_nodes = want_nodes; //number of remaining nodes to allocate
        int i, j, nalloc;
        
        //create new array of structures with elements = #switches
        node_cnt_array = (node_cnt_struct*)malloc(switch_record_cnt * sizeof(node_cnt_struct));

        //populate node_cnt_array
        for (i = 0; i < switch_record_cnt; i++) {
                node_cnt_array[i].switch_idx = i; 
                node_cnt_array[i].free_nodes = switch_node_cnt[i]; 
                //debug("switch:%d freenodes:%lu, switch_node_cnt:%lu", node_cnt_array[i].switch_idx, node_cnt_array[i].free_nodes, switch_node_cnt[i]);
        }

        //sort switch_idx array (containing switch numbers) according to node_cnt array (containing free numbers)
        if (job_ptr->comment && strncmp(job_ptr->comment, "1", 1) == 0)
       		qsort(node_cnt_array, switch_record_cnt, sizeof(*node_cnt_array), desc_cmp_struct);
		else 
			qsort(node_cnt_array, switch_record_cnt, sizeof(*node_cnt_array), inc_cmp_struct);

        debug("Switches with free nodes after sorting::");
        for (i = 0; i < switch_record_cnt; i++) {
                if (node_cnt_array[i].free_nodes > 0)
                        debug("switch:%d freenodes:%lu", node_cnt_array[i].switch_idx, node_cnt_array[i].free_nodes);
        }
        //debug("%d", job_ptr->node_cnt); //gives 0 nodes (it is not populated as of yet)
        int nparts; //number of partitions
        idx_t* part; //array to store partition
        int jobnodes = want_nodes; //number of nodes required to job

        // debug("%d %d", sizeof(uint64_t), sizeof(long double));

        FILE* fcomm = fopen(job_ptr->comment + 2, "r"); //open the file containing comm matrix
        uint64_t* commpattern = (uint64_t*)malloc(jobnodes * jobnodes * sizeof(uint64_t)); //allocate memory for n x n comm matrix
        int edges = 0; //number of edges in graph (each edge is counted twice for undirected graph)
        
        //scanning communication pattern from file (path provided in job comment parameter) and storing in commpattern 2D array
        for (i = 0; i < jobnodes; i++) {
                for (j = 0; j < jobnodes; j++) {
                        fscanf(fcomm, "%ld", &commpattern[i*jobnodes + j]);
                        if (commpattern[i*jobnodes + j] > 0) 
                                edges += 1; //increment number of edges if non-zero communication
                }
        }
        fclose(fcomm); //close file
        debug("base matrix loaded with %d nodes and %d edges", jobnodes, edges);
        
        int* node_alloc = (int*)malloc(jobnodes * sizeof(int)); //boolean array to check which nodes have been allocated to switches
        for (i = 0; i < jobnodes; i++) node_alloc[i] = 0;

        //Communication-intensive job
        if (job_ptr->comment && strncmp(job_ptr->comment, "1", 1) == 0) {
                int partition_level = 0; //to maintain level of partitioning
                
                //if all nodes can be allocated on switch with highest number of free nodes
                if (want_nodes && want_nodes <= node_cnt_array[0].free_nodes) {
                        int switchno = node_cnt_array[0].switch_idx;
                        for (i = 0; i < jobnodes; i++) {
                                //switch2node[switchno][i] = node (=i)
                                switch2node[switchno*jobnodes + i] = i; 
                                node_alloc[i] = 1; //current node is now allocated
                        }
                        switch_alloc_nodes[switchno] = want_nodes; //allocate required nodes from switch with highest free nodes
                        debug("%s: (allinone) found switch:%d for allocation- nodes:%d "
                      		"allocated:%u ", __func__, switchno, node_cnt_array[switchno].free_nodes, switch_alloc_nodes[switchno]);
                        want_nodes = 0;    
                        
                        free(node_cnt_array);
                        free(node_alloc);
                        free(commpattern);
                        return;
                }

                int* groups = (int*)malloc(MAX_PARTITION_LEVEL * jobnodes * sizeof(int)); //groups array to maintain original ranks (will be used when partitioning level >= 1)

                int median_switch = switch_record_cnt - 1; //initialize as switch with least free nodes (should be a non-leaf switch)
                //decrement until we find a leaf switch with free nodes
                while (node_cnt_array[median_switch].free_nodes == 0) median_switch -= 1;
                median_switch /= 2; //get the median switch as (0 + median_switch) / 2
                //median_switch now contains the index of the median switch in the switch_record_cnt array

                //while there are required nodes AND required nodes are greater than median number of free nodes, then partition the graph
                while (want_nodes && want_nodes > node_cnt_array[median_switch].free_nodes && partition_level < MAX_PARTITION_LEVEL) {
                        //if #free_nodes in median_switch=0, find the next switch containing free nodes
                        int median_free_nodes = node_cnt_array[median_switch].free_nodes;
                        
                        // if (median_free_nodes == 0) break;

                        //number of partitions = ceil(total remaining nodes / median number of free nodes)
                        debug("want_nodes:%d median_switch:%d median_free_nodes:%d", want_nodes, node_cnt_array[median_switch].switch_idx, median_free_nodes);
                        nparts = ceil(want_nodes / (float)median_free_nodes);
                        //nparts = 2; //just for testing
                        debug("#partitions: %d", nparts);
                        
                        //find the partitions from the 2D comm matrix
                        // debug("pattern: %s", job_ptr->comment);
                        float part_ufactor = node_cnt_array[0].free_nodes / (float)median_free_nodes; //unbalanced factor for partitions
                        part = matrix_to_partition(want_nodes, edges, commpattern, nparts, part_ufactor); //part[i] gives partition number for node i
                        debug("Partitioning level %d completed", partition_level);

                        part_struct* partsize = (part_struct*)malloc(nparts * sizeof(part_struct));
                        for (int i = 0; i < nparts; i++) {
                                partsize[i].no = i;
                                partsize[i].size = 0;
                        }

                        //for (int i = 0; i < nparts; i++) debug("%d %d", partsize[i].no, partsize[i].size);

                        //get number of nodes in every partition
                        for (int i = 0; i < want_nodes; i++) {
                                partsize[part[i]].size += 1;
                        }
                        

                        // for (int i = 0; i < nparts; i++) debug("%d %d", partsize[i].no, partsize[i].size);
                        qsort(partsize, nparts, sizeof(part_struct), dec_cmp_part); //sort partsize in dec order

                        for (int i = 0; i < nparts; i++) {
                                debug("partition%d size: %d", partsize[i].no, partsize[i].size);
                        }

                        int k = 0; //to iterate switches
                        int* part_alloc = (int*)malloc(nparts * sizeof(int)); //to track if the nodes in the partition were allocated to a switch
                        for (i = 0; i < nparts; i++) part_alloc[i] = 0;

                        //iterate partitions one by one in dec order
                        for (i = 0; i < nparts; i++) {
                                if (partsize[i].size == 0) {
                                        part_alloc[partsize[i].no] = 1;
                                        continue;
                                }
                                //if kth switch (in order of maximum free nodes) can accomodate current partsize, allocate nodes from it
                                if (node_cnt_array[k].free_nodes >= partsize[i].size) {
                                        int switchno = node_cnt_array[k].switch_idx;
                                        //all nodes on partition i are allocated to switch switchno
                                        map_nodes_to_switch(switch2node, switch_alloc_nodes[switchno], node_alloc, partsize[i].no, switchno, part, want_nodes, groups, partition_level, jobnodes); 
                                        switch_alloc_nodes[switchno] += partsize[i].size;
                                        debug("%s: found switch:%d for partition %d, allocation- nodes:%d "
                      		                "allocated:%u ", __func__, switchno, partsize[i].no, node_cnt_array[k].free_nodes, partsize[i].size);
                                        node_cnt_array[k].free_nodes -= partsize[i].size;
                                        part_alloc[partsize[i].no] = 1; //partition allocation complete
                	                k += 1; //go to next switch
                                }
                        }

                        int unalloc_nodes = 0; //number of nodes from partitions that were unallocated
                        for (i = 0; i < nparts; i++) {
                                //if current partition was unallocated
                                if (part_alloc[partsize[i].no] == 0) {
                                        unalloc_nodes += partsize[i].size; //add all nodes from unallocated partition
                                }
                        }
                        //if there are unallocated partitions
                        if (unalloc_nodes > 0) {
                                debug("%d nodes are unallocated", unalloc_nodes);
                                int t = 0;
                                //populate groups array
                                for (int i = 0; i < want_nodes; i++) {
                                        //if current node belongs to unallocated partition
                                        if (part_alloc[part[i]] == 0) {
                                                groups[partition_level*jobnodes + t] = i; //groups[partition_level][t] = i
                                                t += 1; //previous rank i has new rank t now (rank=i in partition_level-1, rank=t in partition_level)
                                        }
                                }
                                // debug("groups array created");

                                //allocate memory for new comm matrix
                                uint64_t* commpattern_new = (uint64_t*)malloc(unalloc_nodes * unalloc_nodes * sizeof(uint64_t)); 
                                //Construct new communication matrix consisting of nodes from unallocated partitions
                                edges = 0;
                                for (int i = 0; i < unalloc_nodes; i++) {
                                        for (int j = 0; j < unalloc_nodes; j++) {
                                                int rank1 = groups[partition_level*jobnodes + i]; //previous rank of new rank i
                                                int rank2 = groups[partition_level*jobnodes + j]; //previous rank of new rank j
                                                //debug("rank1:%d rank2:%d", rank1, rank2);
                                                //commnew[i][j] = comm[rank1][rank2]
                                                commpattern_new[i*unalloc_nodes + j] = commpattern[rank1*want_nodes + rank2]; //populate new comm matrix
                                                if (commpattern_new[i*unalloc_nodes + j] > 0) edges += 1;
                                        }
                                }

                                debug("new commpattern created");
                                free(commpattern); //free old communication matrix
                                commpattern = commpattern_new; //new communication matrix becomes present communication matrix

                                //sort the array of structures containing switch_idx and free nodes again as free_nodes are now updated on each switch
                                qsort(node_cnt_array, switch_record_cnt, sizeof(*node_cnt_array), desc_cmp_struct);
                        }
                        want_nodes = unalloc_nodes; //update want_nodes (remaining nodes to be allocated)
                        partition_level += 1; //increment partition level
                        free(partsize);
                        free(part);
                        free(part_alloc);

                        median_switch = switch_record_cnt - 1; //initialize as switch with least free nodes (should be a non-leaf switch)
                        //decrement until we find a leaf switch with free nodes
                        while (node_cnt_array[median_switch].free_nodes == 0) median_switch -= 1;
                        median_switch /= 2; //get the median switch as (0 + median_switch) / 2
                }
                free(groups); //free groups array
                
                //if paritioning until maximum level still could not allocate all nodes, allocate them accordingly

                median_switch = switch_record_cnt - 1; //initialize as switch with least free nodes (should be a non-leaf switch)
                //decrement until we find a leaf switch with free nodes
                while (node_cnt_array[median_switch].free_nodes == 0) median_switch -= 1;
                median_switch /= 2; //get the median switch as (0 + median_switch) / 2
                        
                //if some nodes are still unallocated AND all the required nodes can be allocated from the switch with median number of free nodes
                if (want_nodes && want_nodes <= node_cnt_array[median_switch].free_nodes) {
                        int median_free_nodes = node_cnt_array[median_switch].free_nodes;
                        int switchno = node_cnt_array[median_switch].switch_idx; //switchno for median switch
                        int t = 0;
                        for (i = 0; i < jobnodes; i++) {
                                //if current node is unallocated
                                if (node_alloc[i] == 0) {
                                        //switch2node[switchno][i] = node (=i)
                                        switch2node[switchno*jobnodes + switch_alloc_nodes[switchno] + t] = i; 
                                        node_alloc[i] = 1;
                                        t += 1;
                                }
                        }
                        switch_alloc_nodes[switchno] += want_nodes; //allocate required nodes from switch with median free nodes
                        
                        debug("%s: (remaining) found switch:%d for allocation- nodes:%d "
                      		"allocated:%u ", __func__, switchno, node_cnt_array[switchno].free_nodes, switch_alloc_nodes[switchno]);
                        want_nodes = 0;
                }

                // for (i = 0; i < switch_record_cnt; i++) {
                //         if (node_cnt_array[i].free_nodes > 0)
                //                 debug("switch:%d freenodes:%lu", node_cnt_array[i].switch_idx, node_cnt_array[i].free_nodes);
                // }
                // debug("remaining nodes: %d", want_nodes);
                //if still some nodes are unallocated, allocate them in increasing order of free nodes sequentially
                i = switch_record_cnt - 1;
                while (want_nodes > 0 && i >= 0) { 
                        //calculate nodes to allocate to present switch
                	nalloc = (node_cnt_array[i].free_nodes < want_nodes) ? node_cnt_array[i].free_nodes : want_nodes; 
                        int switchno = node_cnt_array[i].switch_idx;
                        if (nalloc == 0) {
                                i -= 1;
                                continue;
                        }

                        int t = 0; 
                        for (int j = 0; j < jobnodes; j++) {
                                if (node_alloc[j] == 0) {
                                        //switch2node[switchno][i] = node (=i)
                                        switch2node[switchno*jobnodes + switch_alloc_nodes[switchno] + t] = j; 
                                        node_alloc[j] = 1;
                                        t += 1;
                                }
                                //if number of nodes to allocate on this switch are over, break out of the loop
                                if (t == nalloc) break;   
                        }
                	debug("%s: (sequentially) found switch:%d for allocation- nodes:%d "
                      		"allocated:%u ", __func__, switchno, node_cnt_array[i].free_nodes, nalloc);
                	switch_alloc_nodes[switchno] += nalloc;
                	node_cnt_array[i].free_nodes -= nalloc;
                	want_nodes -= nalloc;
			i--;
        	}
		}
		else {
			//Compute intensive job
			//allocate nodes sequentially starting from switch in ascending order of free nodes
					for(i = 0; (i < switch_record_cnt && want_nodes && node_cnt_array[i].free_nodes); i++) {
							//number of nodes to allocate on current switch
							nalloc = (node_cnt_array[i].free_nodes < want_nodes) ? node_cnt_array[i].free_nodes : want_nodes;
							int switchno = node_cnt_array[i].switch_idx;
							if (nalloc == 0) continue;

							int t = 0; 
							for (j = 0; j < jobnodes; j++) {
									if (node_alloc[j] == 0) {
											//switch2node[switchno][i] = node (=i)
											switch2node[switchno*jobnodes + switch_alloc_nodes[switchno] + t] = j; 
											node_alloc[j] = 1;
											t += 1;
									}
									//if number of nodes to allocate on this switch are over, break out of the loop
									if (t == nalloc) break;   
							}

				debug("%s: (compute) found switch:%d for allocation- nodes:%d "
									"allocated:%u ", __func__, node_cnt_array[i].switch_idx, node_cnt_array[i].free_nodes, nalloc);
				switch_alloc_nodes[node_cnt_array[i].switch_idx] = nalloc;
				node_cnt_array[i].free_nodes -= nalloc;
				want_nodes -= nalloc;
		}
	}

	debug("Comm Balanced allocation complete");
	free(node_cnt_array);
	free(node_alloc);
	free(commpattern);

/*
	for (i = 0; i < switch_record_cnt; i++) {
			if (switch_alloc_nodes[i] == 0) continue;
			debug("SWITCH %d ALLOCATED NODES:-", i);
			for (j = 0; j < switch_alloc_nodes[i]; j++) {
					debug("%d", switch2node[i*jobnodes + j]);
			}
	}
*/
	//debug("returning from combal()");
	return;
}



void hop(struct job_record *job_ptr)
{
	FILE *f;
	char* dir = getenv("SLURM_COST_DIR"); //environment variable SLURM_COST_DIR gives the path where to store the cost files
	//debug("getenv successful with value %s", dir);
	char path_hops[100];
	char path_debug[100];
	strcpy(path_hops, dir);
	strcat(path_hops, "hops.txt");
	strcpy(path_debug, dir);
	strcat(path_debug, "debug.txt");
       
	f = fopen (path_hops, "a");
	int i,j, begin, end,k=0;
	int size = job_ptr->node_cnt;
	int switches[size];
	int index = 0;
	struct node_record *node_ptr;
	int *switch_idx; //To store the result of switch_idx_table lookup
	int *switch_alloc_nodes; //To store the result of alloc_node_table lookup
	int switch_result; //To see if switch_array lookup was successful
	int node_result; // To see if node lookup was successful
	int *n = (int*)malloc(sizeof(int));
	int nunique=0; //No of unique switches (lookup does not give the exact number of switches in array)	
        //debug("Original size:%d, switch levels:%d",size,switch_levels);

	switch_idx = xcalloc(switch_record_cnt, sizeof(int)); // Ordered array of switches selected
        switch_alloc_nodes = xcalloc(switch_record_cnt, sizeof(int)); //Ordered array of nodes allocated on each switch 
	
	switch_result = lookup(switch_idx_table,job_ptr->job_id,switch_idx,n);
	node_result = lookup(alloc_node_table,job_ptr->job_id,switch_alloc_nodes,n);	
	
	if(switch_result == 0 && node_result == 0){
		if (*n != switch_record_cnt)
			debug("Arrays are of inconsistent size");
		debug("Generating arrays");
	}
	else {
		debug("THIS SHOULD NEVER HAPPEN: KEY WAS NOT FOUND");
		begin = bit_ffs(job_ptr->node_bitmap);
        	if (begin >=0 )
                	end = bit_fls(job_ptr->node_bitmap);
        	else
                	end = -1;
       		for (i=begin; i<=end; i++){
                	if(!bit_test(job_ptr->node_bitmap, i))
                        	continue;
                	node_ptr = node_record_table_ptr + i;
                	switches[index]= node_ptr->leaf_switch;
                	//debug("Node name = %s , switches[%d]=%d",
                	//node_ptr->name,index,switches[index]);
                	index+=1;
        	}
	}
	
	// Now add these to strings 
	char job_info[1000]; //For jobname, id, comment, and nunique
	char switch_info[2000]; //For switches
	char alloc_info[2000]; //For allocated nodes
	char comm_info[2000]; //For comm nodes
	char total_info[2000]; //For total nodes
	
	int i_switch=0;
	int i_alloc=0;
	int i_comm=0;
	int i_total=0;
	int i_nodes = 0;
	char node_info[2000];

	for (i = 0; i < switch_record_cnt; i++) {
			//if switch has nodes allocated to it and is a leaf switch
			if (switch_alloc_nodes[i] != 0 && switch_record_table[i].level == 0) {
					nunique += 1;
					i_switch += sprintf(&switch_info[i_switch], "%d ", i); //write switchno
					i_alloc += sprintf(&alloc_info[i_alloc], "%d ", switch_alloc_nodes[i]); //#allocated nodes on switch
					i_comm += sprintf(&comm_info[i_comm], "%d ", switch_record_table[i].comm_jobs); //#commjobs on switch
					i_total += sprintf(&total_info[i_total], "%d ", switch_record_table[i].num_nodes); //#nodes on switch

					//find nodes of the job matrix allocated to the switch
					for (j = 0; j < switch_alloc_nodes[i]; j++) {
							switches[switch2node[i*size + j]] = i;
							i_nodes += sprintf(&node_info[i_nodes], "%d ", switch2node[i*size + j]); //job node present on switch
					}
					i_nodes += sprintf(&node_info[i_nodes], " | ");
			}
	}

	//print newly created switches array
	// for (int i = 0; i < size; i++) {
	// 	debug("Node %d on switch %d", i, switches[i]);
	// }

	// Although for experiments comment should always be given
	if (job_ptr->comment)
                sprintf(job_info,"%s %"PRIu32" %s %d",job_ptr->name,job_ptr->job_id,job_ptr->comment,nunique);
        else
                sprintf(job_info,"%s %"PRIu32" 0 %d",job_ptr->name,job_ptr->job_id,nunique);

	// Add this information to a file
	FILE *info;
	info = fopen(path_debug,"a");
	fputs(job_info,info); // Append jobinfo 
    fprintf(info,"\n");
	
	fputs(switch_info,info);
    fprintf(info,"\n");
	
	fputs(alloc_info,info);
    fprintf(info,"\n");
		
	fputs(node_info, info);
    fprintf(info,"\n");

	fputs(comm_info,info);
    fprintf(info,"\n");

	fputs(total_info,info);
    fprintf(info,"\n");
	fclose(info);

	/**** Writing to debug file over **/	

	char temp[2000];
	long double comm_hops_max = 0;
        comm_hops_max = calc_hops_comm_matrix(switches, size, job_ptr);

        sprintf(temp,"%s %"PRIu32" %s %Lf",job_ptr->name, job_ptr->job_id, job_ptr->comment, comm_hops_max);
        debug("Hops calculated according to communication pattern: %Lf", comm_hops_max);

	fputs(temp,f);
	fprintf(f,"\n");
	fclose(f);
	free(n);
	xfree(switch_idx);
	xfree(switch_alloc_nodes);
}
