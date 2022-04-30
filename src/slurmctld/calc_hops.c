#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>

#include "slurm/slurm.h"

#include "src/common/slurm_topology.h"
#include "src/common/switch.h"
#include "src/common/node_conf.h"

#include "src/slurmctld/job_scheduler.h"
#include "src/slurmctld/proc_req.h"
#include "src/slurmctld/slurmctld.h"

extern int switch_levels;
extern struct switch_record *switch_record_table;
//extern int nodes_per_switch;

/* cnt is total node count */

float fatrecursive(int arr[], int size, int start, int cnt){
	float hops = 0;
	float max_hops = 0;
	int i = 0;
	float c=0, c1=0, c2=0, c3=0;
	for (i =start; i<start +(size/2); i++){
		if ( i+(size/2) < cnt ){
			if (arr[i] == arr [i + (size/2)]){
				c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
				hops=2 + 2*c;
				/*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
					i,i+(size/2),switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
			}
			else{
				c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
				c2 = (switch_record_table[arr[i+(size/2)]].comm_jobs)/((float)switch_record_table[arr[i+(size/2)]].num_nodes);
				c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[i+(size/2)]].comm_jobs)/
					((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[i+(size/2)]].num_nodes);
                                c = c1+c2+c3/2;
				hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
					i,i+(size/2),switch_record_table[arr[i]].comm_jobs,
					switch_record_table[arr[i+(size/2)]].comm_jobs,
					arr[i],arr[i+(size/2)],c,hops);*/
			}
		}
		else
			continue;
		if (hops > max_hops)
			max_hops = hops;
	}
	return max_hops;
}
float treerecursive(int arr[], int size, int start, int cnt){
        float hops = 0;
        float max_hops =0;
        int i=0;
        float c=0, c1=0, c2=0, c3=0;
	for (i =start; i<start +(size/2); i++){
                if ( i+(size/2) < cnt ){
                        if (arr[i] == arr [i + (size/2)]){
			  	c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
				hops=2 + 2*c;
                                /*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
                                        i,i+(size/2),switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
                        }
                        else{
                                c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
                                c2 = (switch_record_table[arr[i+(size/2)]].comm_jobs)/((float)switch_record_table[arr[i+(size/2)]].num_nodes);
                                c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[i+(size/2)]].comm_jobs)/
                                        ((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[i+(size/2)]].num_nodes);
                                c = c1+c2+c3;
				hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
                                        i,i+(size/2),switch_record_table[arr[i]].comm_jobs,
                                        switch_record_table[arr[i+(size/2)]].comm_jobs,
                                        arr[i],arr[i+(size/2)],c,hops);*/
                        }
                }
                else    
                        continue;
                if (hops > max_hops)
                        max_hops = hops;
        }
        return max_hops;
}

float fatreduce(int arr[], int size, int start, int cnt){
	float hops = 0;
	float max_hops = 0;
	int i = 0;
	float c=0, c1=0, c2=0, c3=0;
	for (i=start; i < cnt-size; i+=2*size){
		if (arr[i] == arr[i+size]){
                                c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
		       		hops=2 + 2*c;
                                /*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
		}
		else{
				c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
				c2 = (switch_record_table[arr[i+size]].comm_jobs)/((float)switch_record_table[arr[i+size]].num_nodes);
				c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[i+size]].comm_jobs)/
					((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[i+size]].num_nodes);
				c = c1+c2+c3/2;
		       		hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,
                                        switch_record_table[arr[i+size]].comm_jobs,
                                        arr[i],arr[i+size],c,hops);*/
		}
		
		if (hops > max_hops)
			max_hops = hops;
	}
	return max_hops;
}

float treereduce(int arr[], int size, int start, int cnt){
        float hops = 0;
        float max_hops = 0;
        int i = 0;
        float c=0, c1=0, c2=0, c3=0;
        for (i=start; i < cnt-size; i+=2*size){
                if (arr[i] == arr[i+size]){
                                c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
				hops=2 + 2*c;
                                /*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
                }
                else{
				c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
                                c2 = (switch_record_table[arr[i+size]].comm_jobs)/((float)switch_record_table[arr[i+size]].num_nodes);
                                c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[i+size]].comm_jobs)/
                                        ((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[i+size]].num_nodes);
                                c = c1+c2+c3;
                                hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,
                                        switch_record_table[arr[i+size]].comm_jobs,
                                        arr[i],arr[i+size],c,hops);*/
                }

                if (hops > max_hops)
                        max_hops = hops;
        }
        return max_hops;
}

float binomial(int arr[], int size, int cnt){
        float hops = 0;
        float max_hops = 0;
        int i = 0;
        float c=0, c1=0, c2=0, c3=0;
        for (i = 0; i<size; i++){
                if (arr[i] == arr[i+size]){
                                c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
                                hops=2 + 2*c;
                                /*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
                }
                else{
                                c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
                                c2 = (switch_record_table[arr[i+size]].comm_jobs)/((float)switch_record_table[arr[i+size]].num_nodes);
                                c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[i+size]].comm_jobs)/
                                        ((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[i+size]].num_nodes);
                                c = c1+c2+c3/2;
                                hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
                                        i,i+size,switch_record_table[arr[i]].comm_jobs,
                                        switch_record_table[arr[i+size]].comm_jobs,
                                        arr[i],arr[i+size],c,hops);*/
                }

                if (hops > max_hops)
                        max_hops = hops;
        }
        return max_hops;
}

float ring(int arr[], int cnt){
        float hops = 0;
        float max_hops = 0;
        int i = 0;
        float c=0, c1=0, c2=0, c3=0;
        for (i=0; i<cnt; i++){
                int j = (i+1)%cnt;
                if (arr[i] == arr[j]){
                                c = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes) ;
                                hops=2 + 2*c;
                                /*debug("%d<->%d : Comm_Jobs=%d Contention=%f Hops=%f Switch =%d",
                                        i,j,switch_record_table[arr[i]].comm_jobs,c,hops,arr[i]);*/
                }
                else{
                                c1 = (switch_record_table[arr[i]].comm_jobs)/((float)switch_record_table[arr[i]].num_nodes);
                                c2 = (switch_record_table[arr[j]].comm_jobs)/((float)switch_record_table[arr[j]].num_nodes);
                                c3 = (switch_record_table[arr[i]].comm_jobs + switch_record_table[arr[j]].comm_jobs)/
                                        ((float)switch_record_table[arr[i]].num_nodes + (float)switch_record_table[arr[j]].num_nodes);
                                c = c1+c2+c3/2;
                                hops=2*(switch_levels+1) + 2*(switch_levels+1)*c ;
                                /*debug("%d<->%d : Comm_jobs=%d,%d Switch =%d,%d Contention=%f Hops=%f",
                                        i,j,switch_record_table[arr[i]].comm_jobs,
                                        switch_record_table[arr[j]].comm_jobs,
                                        arr[i],arr[j],c,hops);*/
                }

                if (hops > max_hops)
                        max_hops = hops;
        }
        return (max_hops*(cnt-1));
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
}

void hop(struct job_record *job_ptr)
{
	FILE *f;
	f = fopen ("/home/gagandeep/slurmcost/hops.txt", "a");
	int i, begin, end;
	int size = job_ptr->node_cnt;
	int switches[size];
	int index = 0;
	struct node_record *node_ptr;
	int nunique = 0; // Find number of unique leaf switches
	int prev = -1; // Keep track of previous switch, switches come in ordered form

        size = pow(2,ceil(log(size)/log(2)));
        //debug("Original size:%d, switch levels:%d",size,switch_levels);


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
		if (prev != node_ptr->leaf_switch){
			// Leaf switch has changed
			nunique++;
			prev = node_ptr->leaf_switch;
		}
		/*debug("Node name = %s , switches[%d]=%d",
			node_ptr->name,index,switches[index]);*/
		index+=1;
	}
	/**** Writing to the debug file ***/
	// Now create the consolidated arrays
	int switch_arr[nunique];
	int alloc_arr[nunique];
	int comm_arr[nunique];
	int total_nodes[nunique];

	int j=-1; // index for switch info and other arrays
	prev=-1; // Keep track of previous switch
	for (i=0;i<size;i++){
		if (switches[i]!=prev){
			// New leaf switch
			j++; //Move to next entry
			prev = switches[i];
			switch_arr[j] = switches[i];
			comm_arr[j] = switch_record_table[switches[i]].comm_jobs;
			total_nodes[j] = switch_record_table[switches[i]].num_nodes;
			alloc_arr[j] = 1;// Switch has occured for the first time	
		}
		else
			alloc_arr[j]++;	
	}
	// Now add these to strings 
	char job_info[100]; //For jobname, id, comment, and nunique
	char switch_info[2000]; //For switches
	char alloc_info[2000]; //For allocated nodes
	char comm_info[2000]; //For comm nodes
	char total_info[2000]; //For total nodes
	
	// Although for experiments comment should always be given
	if (job_ptr->comment)
                sprintf(job_info,"%s %"PRIu32" %s %d",job_ptr->name,job_ptr->job_id,job_ptr->comment,nunique);
        else
                sprintf(job_info,"%s %"PRIu32" 0 %d",job_ptr->name,job_ptr->job_id,nunique);
	int i_switch=0;
	int i_alloc=0;
	int i_comm=0;
	int i_total=0;

	for (i=0;i<nunique;i++){
		i_switch += sprintf(&switch_info[i_switch], "%d ", switch_arr[i]);
		i_alloc += sprintf(&alloc_info[i_alloc], "%d ", alloc_arr[i]);
		i_comm += sprintf(&comm_info[i_comm], "%d ", comm_arr[i]);
		i_total += sprintf(&total_info[i_total], "%d ", total_nodes[i]);
	}
	//debug("%s",job_info);
	//debug("%s", switch_info);
	//debug("%s", alloc_info);
	//debug("%s", comm_info);
	//debug("%s", total_info);	

	// Add this information to a file
	FILE *info;
	info = fopen("/home/gagandeep/slurmcost/debug.txt","a");
	fputs(job_info,info); // Append jobinfo 
        fprintf(info,"\n");
	
	fputs(switch_info,info);
        fprintf(info,"\n");
	
	fputs(alloc_info,info);
        fprintf(info,"\n");

	fputs(comm_info,info);
        fprintf(info,"\n");

	fputs(total_info,info);
        fprintf(info,"\n");
	fclose(info);

	/**** Writing to debug file over **/

	float hops = 0;
	float max_hops =0;

// Calculate Hops for recursive halving
	float rec_fathops =0;
	int rec_size = size;
	int msize =1; // Message size for recursive halving calculations
	//debug("Calculating fat tree recursive hops");
	while(rec_size > 1){
		max_hops = 0;
		for (i=0; i<job_ptr->node_cnt; i+= rec_size){
			hops = fatrecursive(switches,rec_size,i,job_ptr->node_cnt);
			if (hops > max_hops)
				max_hops = hops;
		}
		//debug(" rec_fathops = %d x %f ",msize,max_hops);
		rec_fathops += msize * max_hops;
		msize = msize * 2;
		rec_size = rec_size /2;
	}

	float rec_treehops =0;
        rec_size = size;
        msize =1; // Message size for recursive halving calculations
        //debug("Calculating tree recursive hops");
        while(rec_size > 1){
                max_hops = 0;
                for (i=0; i<job_ptr->node_cnt; i+= rec_size){
                        hops = treerecursive(switches,rec_size,i,job_ptr->node_cnt);
                        if (hops > max_hops)
                                max_hops = hops;
                }
                //debug(" rec_treehops = %d x %f ",msize,max_hops);
                rec_treehops += msize * max_hops;
                msize = msize * 2;
                rec_size = rec_size /2;
        }


// Calculate Hops for reduce
	float red_fathops =0;
	int red_size = 1;
	//debug("Calculating fat tree reduce hops");
	while(red_size < size){
		red_fathops += fatreduce(switches,red_size,0,job_ptr->node_cnt);
		red_size *=2;
	}
        float red_treehops =0;
        red_size = 1;
        //debug("Calculating tree reduce hops");
        while(red_size < size){
                red_treehops += treereduce(switches,red_size,0,job_ptr->node_cnt);
                red_size *=2;
        }
// Calculate binomial hops
        float bin_hops = 0;
        int bin_size = 1;
        //debug("Calculating binomial hops");
        while (bin_size < size){
                bin_hops += binomial(switches, bin_size, job_ptr->node_cnt);
                bin_size *=2;
        }
// Ring Hops
        //debug("Calculating ring hops");
        float ring_hops = ring(switches, job_ptr->node_cnt);

        char temp[1000];

		//if communication matrix path is given
		if (strlen(job_ptr->comment) > 3) {
			long double comm_hops = calc_hops_comm_matrix(switches, size, job_ptr);
			sprintf(temp,"%s %"PRIu32" %s %Lf", job_ptr->name, job_ptr->job_id, job_ptr->comment, comm_hops);
			debug("Hops calculated according to communication matrix: %Lf", comm_hops);
		}
		//comment is either 0/1/1:x
		else {
			if (job_ptr->comment)
					sprintf(temp,"%s %"PRIu32" %s %f %f %f %f %f %f",job_ptr->name,job_ptr->job_id,job_ptr->comment,rec_fathops,rec_treehops,red_fathops,red_treehops,bin_hops,ring_hops);
			else
					sprintf(temp,"%s %"PRIu32" 0 %f %f %f %f %f %f",job_ptr->name,job_ptr->job_id,rec_fathops,rec_treehops,red_fathops,red_treehops,bin_hops,ring_hops);

			debug("Recursive FatHops:%f TreeHops:%f | Reduce FatHops = %f Treehops =%f | Binomial:%f Ring:%f | temp: %s",rec_fathops,rec_treehops,red_fathops,red_treehops,bin_hops,ring_hops,temp);
		}

	fputs(temp,f);
	fprintf(f,"\n");
	fclose(f);
}
