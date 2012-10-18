#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#include "yaxt.h"

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))
//#define DEBUG

struct config {

   int x, y, z; // global data size (default 1024 * 512 * 1)

   struct {
      int x_leng, y_leng, z_leng;        //local data size
      unsigned local_size[3];            
      unsigned local_start[3];
   } src, dst; //!< source and target decomposion
};

void print_index(Xt_idxlist idxsection, int rank, int i){
  #ifdef DEBUG
  for (int p = 0; p < xt_idxlist_get_num_indices(idxsection); p ++){
      Xt_idx value;
      xt_idxlist_get_index (idxsection, p, &value);
      
      if (i == 0)
        printf ("rank: %d , INDEX_SRC: %d \n", rank, value);
      else 
        printf ("rank: %d , INDEX_DST: %d \n", rank, value);
    }
  #endif
}

void calc_blocks (int nprocs, int dim2, int dim1, int *dim2_block, int *dim1_block){
  // global data proportion
  double prop = 1;
  // temp num of blocks in x direction 
  int tmp_dim1 = -1;
  // temp prop to compare last and actual iteration
  double tmp_prop = -1.0; 
  
  int snprocs = sqrt(nprocs);

  // calculate global data proportion by division of bigger number by smaller number, so prop <= 1
  prop = (dim2 < dim1) ? (double)dim2/(double)dim1 : (double)dim1/(double)dim2;
    
  //
  for (int f = snprocs; f > 0; f-- ){
    if (nprocs % f == 0){    
      *dim1_block = nprocs / f;
      *dim2_block = f;
        
      double actu = fabs (prop - (double)*dim2_block/(double)*dim1_block);
      double past = fabs  (prop - tmp_prop) ;
      
      //compare if actual prop better then last. If, make it to actual and do next iteration, otherwise take last and break
      if (actu < past ){
        tmp_prop = (double)*dim2_block / (double)*dim1_block;
        tmp_dim1 = *dim1_block;
      }
      else {
        *dim1_block = tmp_dim1;
        *dim2_block = nprocs / tmp_dim1;
        break;
      }
    }
  }
  
  // adapt proportion to global data
  if(dim2 > dim1){
    *dim1_block = *dim2_block;
    *dim2_block = tmp_dim1;   
  }
}

void init_data(int *data_array, int array_leng, int j){
  // 0 = source
  if (j == 0)
    for (int i = 0; i < array_leng; i++) {
      data_array[i] = 100 + i;
    }
  // 1 = dest
  else
    for (int i = 0; i < array_leng; i++) {
      data_array[i] = -1;
    }
}

void calc_config (struct config *conf, int nprocs){
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  
  // --- SRC
  // source num blocks in global structure
  int src_y_block = 1;
  int src_x_block = 1;
  
  // calculate num of blocks regarding global data proportion
  calc_blocks (nprocs, conf->y, conf->x, &src_y_block, &src_x_block);
  
  // source local data size in x, y and z direction - num of elements

  conf->src.z_leng = conf->z;
  conf->src.y_leng = conf->y / src_y_block + ((rank / src_x_block) < (conf->y % src_y_block));
  conf->src.x_leng = conf->x / src_x_block + ((rank % src_x_block) < (conf->x % src_x_block));
  
  // --- DEST
  //destination num blocks in global structure
  int dst_z_block = 1;
  int dst_y_block = 1;
  //int dst_x_block = 1;
    
  calc_blocks (nprocs, conf->z, conf->y, &dst_z_block, &dst_y_block);
    
  //calculate actually length of local destination data array 
  conf->dst.z_leng = conf->z / dst_z_block + ((rank / dst_y_block) < (conf->z % dst_z_block));
  conf->dst.y_leng = conf->y / dst_y_block + ((rank % dst_y_block) < (conf->y % dst_y_block));
  conf->dst.x_leng = conf->x;

  // calculate parameters for source section
  conf->src.local_size[0] = conf->z;
  conf->src.local_size[1] = conf->src.y_leng;
  conf->src.local_size[2] = conf->src.x_leng;
  int src_block_pos_y = rank / src_x_block;
  int src_block_pos_x = rank % src_x_block;
     
  // how many blocks get an extra row or column (source_dim_blocks_bonus)
  int src_y_bb = conf->y % src_y_block;
  int src_x_bb = conf->x % src_x_block;
  
  conf->src.local_start[0] = 0;
  conf->src.local_start[1] = (conf->y/ src_y_block) * src_block_pos_y + MIN(src_block_pos_y, src_y_bb);
  conf->src.local_start[2] = (conf->x/ src_x_block) * src_block_pos_x + MIN(src_block_pos_x, src_x_bb);

  //dst
  conf->dst.local_size[0] = conf->dst.z_leng;
  conf->dst.local_size[1] = conf->dst.y_leng;
  conf->dst.local_size[2] = conf->x;

  int dst_block_pos_z = rank / dst_y_block;
  int dst_block_pos_y = rank % dst_y_block;    
    
  // how many blocks get an extra row or column (source_dim_blocks_bonus)
  int dst_z_bb = conf->z % dst_z_block;   
  int dst_y_bb = conf->y % dst_y_block;

  conf->dst.local_start[0] = (conf->z/ dst_z_block) * dst_block_pos_z + MIN(dst_block_pos_z, dst_z_bb);
  conf->dst.local_start[1] = (conf->y/ dst_y_block) * dst_block_pos_y + MIN(dst_block_pos_y, dst_y_bb);
  conf->dst.local_start[2] = 0;                              
}


void get_config(struct config *conf, int argc, char** argv, int nprocs){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int map = 0;
  int par = 0;
  
  while(( par = getopt(argc, argv, "x:y:z:")) != -1) {
      switch(par) {
          case 'x':
              conf->x = atoi(optarg);
              map |= 1;
              break;
          case 'y':
              conf->y = atoi(optarg);
              map |= 2;
              break;
          case 'z':
             conf->z = atoi(optarg);
             map |= 4;
              break;
          case '?': 
             //printf("No arguments \n");
             break;
      }
  }

  if (argc == 1){
    if(rank == 0)
      printf("No arguments, set default: x=1024, y=512, z=1 \n");  
    conf->x = 1024;
    conf->y = 512;
    conf->z = 1;
  }else if((map != 7) || (optind != 7 )) {
    if (rank == 0)
      printf(" Option fail \n Example: ./tr_gp2fs -x 10 -y 20 -z 100 \n Abort \n");
    exit(EXIT_FAILURE);
  }
  
  calc_config (conf, nprocs);
}

int main(int argc, char** argv) {
  
  //init mpi
 
  MPI_Init(NULL, NULL);
 
  xt_initialize (MPI_COMM_WORLD);
  
  int nprocs, rank;
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  
  struct config conf;  

  {
    get_config(&conf, argc, argv, nprocs);
    
    // --- local data array
    int *src_array, *dst_array; 
    
    //allocate memory for source
    int src_leng = conf.src.z_leng * conf.src.y_leng * conf.src.x_leng;
    src_array = (int*)malloc(sizeof(*src_array)* src_leng);
  
    // allocate memory of actuall local data size
    int dst_leng = conf.dst.z_leng * conf.dst.y_leng * conf.dst.x_leng;
    dst_array = (int*)malloc(sizeof(dst_array) * dst_leng);  
   
    // init arrays with fake data
    init_data(src_array, src_leng, 0);
    init_data(dst_array, dst_leng, 1);   
  
    // global param
    unsigned num_dims = 3;
    unsigned global_size[3] = {conf.z, conf.y, conf.x};   
    
    // -- SOURCE index list by section
    Xt_idxlist src_idxsection;
    Xt_idx start = 0;

    //create source index section
    src_idxsection = xt_idxsection_new(start, num_dims, global_size, conf.src.local_size, conf.src.local_start);
    
    #ifdef DEBUG
    // 0 = soruce
    print_index(src_idxsection, rank, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
    // --- DEST-ination index list by section
    Xt_idxlist dst_idxsection;

    dst_idxsection = xt_idxsection_new(start, num_dims, global_size, conf.dst.local_size, conf.dst.local_start);

    // 1 = dest
     print_index(dst_idxsection, rank, 1);
    
    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxsection, dst_idxsection, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_INTEGER);

    //array poiter, especially necessary for data array number > 1
    int* src_array_p = src_array;
    int* dst_array_p = dst_array;

    //Exchange
    xt_redist_s_exchange(redist, (void**)(&src_array_p), 1, (void**)(&dst_array_p), 1);
        
    //clean up
    free (src_array);
    free (dst_array);
    
    xt_idxlist_delete(dst_idxsection);
    xt_idxlist_delete(src_idxsection);
    xt_xmap_delete(xmap);
    xt_redist_delete(redist);
  }
  
  MPI_Finalize();

  return 0;
}

