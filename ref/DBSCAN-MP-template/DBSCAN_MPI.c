/* Copyright 2015 Gagarine Yaikhom (MIT License) */
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "point.h"
#include "constants.h"
#include "list.h"
#include "mpi.h"
#include <math.h>

index_list *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon);
index_list *get_epsilon_neighbours2(
    double coX, double coY, double coZ,
    point_t *local,
    unsigned int num,
    double epsilon);
void print_epsilon_neighbours(
    point_t *points,
    index_list *en);
void destroy_epsilon_neighbours(index_list *en);
index_list *dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts);
int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts);
int spread(
    unsigned int index,
    index_list *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts);
int **get_pairs(int step, int size);

void print_epsilon_neighbours(
    point_t *points,
    index_list *en)
{
    if (en) {
        node_t *h = en->head;
        while (h) {
            printf("(%lfm, %lf, %lf)\n",
                   points[h->index].x,
                   points[h->index].y,
                   points[h->index].z);
            h = h->next;
        }
    }
}

index_list *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon)
{
    index_list *en = (index_list *)
        calloc(1, sizeof(index_list));
    if (en == NULL)
    {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }
    for (int i = 0; i < num_points; ++i)
    {
        if (i == index)
            continue;
        if (dist(points[index], points[i]) > epsilon)
            continue;
        else
        {
            if (append_at_end(i, en) == FAILURE)
            {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}
//redefinition of function. Ability to provide random coordinates is given
index_list *get_epsilon_neighbours2(
    double coX, double coY, double coZ,
    point_t *local,
    unsigned int num,
    double epsilon)
{
    index_list *en = (index_list *)
        calloc(1, sizeof(index_list));
    if (en == NULL)
    {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }
    for (int i = 0; i < num; ++i)
    {
        if (dist2(coX, coY, coZ, local[i].x, local[i].y, local[i].z) > epsilon)
            continue;
        else
        {
            if (append_at_end(i, en) == FAILURE)
            {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}

void destroy_epsilon_neighbours(index_list *en)
{
    if (en)
    {
        node_t *t, *h = en->head;
        while (h)
        {
            t = h->next;
            free(h);
            h = t;
        }
        free(en);
    }
}

index_list *dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts)
{
    index_list *cp = (index_list *)
        calloc(1, sizeof(index_list));
    unsigned int i, cluster_id = 0;
    for (i = 0; i < num_points; ++i)
    {
        if (points[i].cluster_id == UNCLASSIFIED)
        {
            int is_core_point = expand(i, cluster_id, points,
                                       num_points, epsilon, minpts);
            if (is_core_point == CORE_POINT)
            {
                ++cluster_id;
                append_at_end(i, cp);
            }
        }
    }
    return cp;
}

int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts)
{
    int return_value = NOT_CORE_POINT;
    index_list *seeds =
        get_epsilon_neighbours(index, points,
                               num_points, epsilon);
    if (seeds == NULL)
        return FAILURE;

    if (seeds->num_members < minpts)
        points[index].cluster_id = NOISE;
    else
    {
        points[index].cluster_id = cluster_id;
        node_t *h = seeds->head;
        while (h)
        {
            points[h->index].cluster_id = cluster_id;
            h = h->next;
        }

        h = seeds->head;
        while (h)
        {
            spread(h->index, seeds, cluster_id, points,
                   num_points, epsilon, minpts);
            h = h->next;
        }

        return_value = CORE_POINT;
    }
    destroy_epsilon_neighbours(seeds);
    return return_value;
}

int spread(
    unsigned int index,
    index_list *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts)
{
    index_list *spread =
        get_epsilon_neighbours(index, points,
                               num_points, epsilon);
    if (spread == NULL)
        return FAILURE;
    if (spread->num_members >= minpts)
    {
        node_t *n = spread->head;
        point_t *d;
        while (n)
        {
            d = &points[n->index];
            if (d->cluster_id == NOISE ||
                d->cluster_id == UNCLASSIFIED)
            {
                if (d->cluster_id == UNCLASSIFIED)
                {
                    if (append_at_end(n->index, seeds) == FAILURE)
                    {
                        destroy_epsilon_neighbours(spread);
                        return FAILURE;
                    }
                }
                d->cluster_id = cluster_id;
            }
            n = n->next;
        }
    }
    destroy_epsilon_neighbours(spread);
    return SUCCESS;
}

//function to get radius of each cluster
void get_core_radius(point_t *points, unsigned int num_of_points, index_list *cp)
{
    node_t *temp_core_point = cp->head;
    unsigned int i;
    while (temp_core_point)
    {
        unsigned int cp_i = temp_core_point->index;
        double rad = 0;
        for (i = 0; i < num_of_points; i++)
        {
            if (points[i].cluster_id == points[cp_i].cluster_id)
            {
                double distance = dist(points[i], points[cp_i]);
                if (rad < distance)
                    rad = distance;
            }
        }
        temp_core_point->radius = rad;
        temp_core_point = temp_core_point->next;
    }
}

int main(int argc, char *argv[])
{
    int rank, size;
    int tag0 = 40;
    int tag1 = 50;
    int tag2 = 60;
    double tStart, tStop;
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int blocklengths[4] = {1, 1, 1, 1};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Datatype mpi_point_t_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(point_t, x);
    offsets[1] = offsetof(point_t, y);
    offsets[2] = offsetof(point_t, z);
    offsets[3] = offsetof(point_t, cluster_id);

    MPI_Type_create_struct(4, blocklengths, offsets, types, &mpi_point_t_type);
    MPI_Type_commit(&mpi_point_t_type);

    //Read file of points coordinates into directory
    FILE *in = fopen("points.dat", "r");

    point_t *points;
    double epsilon;
    unsigned int minpts, i;
    unsigned int num_points =
        parse_input(in, &points, &epsilon, &minpts);
    fclose(in);
    if(rank == 0){
        tStart = clock();
    }
    unsigned int num = num_points/size;
    point_t *myLocal;
    myLocal = (point_t *)calloc(num_points,sizeof(point_t));
    
    //Scatter elements across MPI_COMM_WORLD deviding points set into equal size subsets
    MPI_Scatter(points, num, mpi_point_t_type, myLocal, num, mpi_point_t_type, 0, MPI_COMM_WORLD);

    if (num)
    {
        //Classify points of subset while finding the center of each cluster. Then save it into cp array.
        index_list *cp = (index_list *)calloc(num, sizeof(index_list));
        cp = dbscan(myLocal, num, epsilon,
                            minpts);

        //add radius of each cluster into cp array
        get_core_radius(myLocal, num, cp);

        /*printf("Epsilon: %lf\n", epsilon);
        printf("Minimum points: %u\n", minpts);
        print_points(myLocal, num);
        print_list(points, cp);*/

        //convert cp array to double array and add cluster_id of the cluster to which it refers to.
        double (*myArr)[5];
        myArr = calloc(cp->num_members*5, sizeof(double*));
        myArr[cp->num_members][5];
        unsigned int i = 0;
        node_t *h = cp->head;
        while(h){
            myArr[i][0] = myLocal[h->index].x;
            myArr[i][1] = myLocal[h->index].y;
            myArr[i][2] = myLocal[h->index].z;
            myArr[i][3] = h->radius;
            myArr[i][4] = myLocal[h->index].cluster_id;
            i++;
            h = h->next;
        }
        
        //preallocation and predefinition of arrays and variables respectively.
        double (*newArr)[5];      //array to be produced after cluster comparison between processes
        point_t *receivedLocal;   //received array to be copared with myArr
        double (*receivedArr)[5];
        int newCountOfCores;
        int numOfCores;

        //for loop to count steps with respect to the number of processors which are being involved
        for(int step=0; step<(int)log2(size); step++){
            //step 0 is executed seperately from the rest steps
            if(step==0){
                int senderRank, recverRank;
                int n = size/2;
                for(int i=0; i<n; i++){        //Characterizing the processes as sender and receiver
                    senderRank = 2*i;
                    recverRank = 2*i+1;
                    if(rank == senderRank){
                        int countArr0 = cp->num_members*5;
                        //Send myArr (with clusters centroids radiuses and ids) as well as the local subset of points to be classified
                        MPI_Send(&myArr[0][0], countArr0, MPI_DOUBLE, recverRank, tag1, MPI_COMM_WORLD);
                        MPI_Send(myLocal, num, mpi_point_t_type, recverRank, tag2, MPI_COMM_WORLD);
                    }else if(rank == recverRank){
                        //Count the number of elements that are being received
                        int count;
                        MPI_Status status;
                        MPI_Probe(senderRank, tag1, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_DOUBLE, &count);
                        
                        //allocation for received arrays
                        receivedArr = calloc(count, sizeof(double*));
                        receivedLocal = (point_t *)calloc(num,sizeof(point_t));
                        //receiving arrays
                        MPI_Recv(receivedArr, count, MPI_DOUBLE, senderRank, tag1, MPI_COMM_WORLD, &status);
                        MPI_Recv(receivedLocal, num, mpi_point_t_type, senderRank, tag2, MPI_COMM_WORLD, &status);
                        
                        //Finde maxId of clusters classification whithin the receiver process
                        int maxId = 0;
                        for(int j=0; j<num; j++){
                            if(myLocal[j].cluster_id > maxId){
                                maxId = myLocal[j].cluster_id;
                            }
                        }
                        for(int j=0; j<num; j++){
                            receivedLocal[j].cluster_id = receivedLocal[j].cluster_id + maxId + 1;
                        }
                        for(int j=0; j<count/5; j++){
                            receivedArr[j][4] = receivedArr[j][4] + maxId + 1;
                        }
                        maxId = 0;
                        for(int j=0; j<num; j++){
                            if(receivedLocal[j].cluster_id > maxId){
                                maxId = receivedLocal[j].cluster_id;
                            }
                        }
                        
                        //allocation of merged array and begin the merging of two subsets of points
                        point_t *mergeLocal;
                        mergeLocal = (point_t *)calloc(2*num,sizeof(point_t));
                        
                        int i=0;                        
                        while(i<2*num){
                            if(i<num){
                                mergeLocal[i].x  = receivedLocal[i].x;
                                mergeLocal[i].y  = receivedLocal[i].y;
                                mergeLocal[i].z  = receivedLocal[i].z;
                                mergeLocal[i].cluster_id  = receivedLocal[i].cluster_id;
                            }else{
                                mergeLocal[i].x  = myLocal[i-num].x;
                                mergeLocal[i].y  = myLocal[i-num].y;
                                mergeLocal[i].z  = myLocal[i-num].z;
                                mergeLocal[i].cluster_id  = myLocal[i-num].cluster_id;
                            }
                            i++;
                        }

                        //allocation of merged array and begin the merging of two local arrays of clusters characteristics
                        newCountOfCores = cp->num_members + count/5;
                        newArr = calloc(newCountOfCores*5, sizeof(double*)); 
                        int j=0;
                        while(j<newCountOfCores){
                            for(int l=0; l<cp->num_members; l++){
                                newArr[j][0] = myArr[l][0];
                                newArr[j][1] = myArr[l][1];
                                newArr[j][2] = myArr[l][2];
                                newArr[j][3] = myArr[l][3];
                                newArr[j][4] = myArr[l][4];
                                j++;
                            }
                            for(int l=0; l<count/5; l++){
                                newArr[j][0] = receivedArr[l][0];
                                newArr[j][1] = receivedArr[l][1];
                                newArr[j][2] = receivedArr[l][2];
                                newArr[j][3] = receivedArr[l][3];
                                newArr[j][4] = receivedArr[l][4];
                                j++;
                            }
                        }  
                        //case 1: one lies entirely within the other;
                        for(i=0; i<cp->num_members; i++){
                            for(int r=0; r<count/5; r++){
                                double d = (double)sqrt(pow(myArr[i][0]-receivedArr[r][0], 2)+pow(myArr[i][1]- //d is the distance between centroids
                                                receivedArr[r][1], 2)+pow(myArr[i][2]-receivedArr[r][2], 2));

                                double radDif = (double)abs(myArr[i][3] - receivedArr[r][3]);                  //radDif is the difference between correspondant radiuses
                                //First comparison. 
                                if(d<=radDif){
                                    //Depending the result, one cluster will take the id of the other
                                    if(myArr[i][3]>receivedArr[r][3]){
                                        for(int j=0; j<num; j++){
                                            if(receivedLocal[j].cluster_id == (int)receivedArr[r][4]){
                                                receivedLocal[j].cluster_id = (int)myArr[i][4];
                                            }
                                        }
                                        for(int i=0; i<newCountOfCores; i++){
                                            for(int j=0; j<5; j++){
                                                for(int k=0; k<count/5; k++){
                                                    for(int p=0; p<5; p++){
                                                        if(newArr[i][j] == receivedArr[k][p]){
                                                            newArr[i][j] = 0;
                                                        }
                                                    }
                                                }
                                            }
                                        }                                            
                                    }else{
                                        for(int j=0; j<num; j++){
                                            if(myLocal[j].cluster_id == myArr[i][4]){
                                                myLocal[j].cluster_id = (int)receivedArr[r][4];
                                            }
                                        }
                                        for(int i=0; i<newCountOfCores; i++){
                                            for(int j=0; j<5; j++){
                                                for(int k=0; k<cp->num_members; k++){
                                                    for(int p=0; p<5; p++){
                                                        if(newArr[i][j] == myArr[k][p]){
                                                            newArr[i][j] = 0;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    for(int j=0; j<2*num; j++){
                                        for(int l=0; l<num; l++){
                                            if(mergeLocal[j].x == myLocal[l].x && mergeLocal[j].y == myLocal[l].y && mergeLocal[j].z == myLocal[l].z){
                                                mergeLocal[j].cluster_id = myLocal[l].cluster_id;
                                            }
                                        }
                                        for(int l=0; l<num; l++){
                                            if(mergeLocal[j].x == receivedLocal[l].x && 
                                                mergeLocal[j].y == receivedLocal[l].y && mergeLocal[j].z == receivedLocal[l].z){
                                                mergeLocal[j].cluster_id = receivedLocal[l].cluster_id;
                                            }
                                        }
                                    } 
                                }
                            }
                        }                        
                        
                        int n=1;
                        for(int i=0; i<cp->num_members; i++){
                            int counter = 1;
                            double (*collidedClusters)[5];
                            collidedClusters = calloc(20, sizeof(double*));
                            for(int r=0; r<count/5; r++){
                                //case 2: intersect in two points;
                                double d = (double)sqrt(pow(myArr[i][0]-receivedArr[r][0], 2)+
                                            pow(myArr[i][1]-receivedArr[r][1], 2)+pow(myArr[i][2]-receivedArr[r][2], 2));
                                double radConcat = myArr[i][3] + receivedArr[r][3]; //radConcat stands for the concatenation of correspondant radiuses
                                if(d<radConcat){
                                    int count1 = 0;
                                    int count2 = 0;
                                    int j = 0;
                                    int k = 0;
                                    int l = 0;
                                    while(j<num){
                                        if(myLocal[j].cluster_id == myArr[i][4]){
                                            count1++;
                                        }
                                        j++;
                                    }
                                    while(k<num){
                                        if(receivedLocal[k].cluster_id == (int)receivedArr[r][4]){
                                            count2++;
                                        }
                                        k++;
                                    }
                                    //verify if points of same cluster_id's exist in the correspondant clusters under consideration 
                                    if(count1 == 0 || count2 == 0){
                                        continue;
                                    }else{
                                        //save collided clusters into one array
                                        counter = counter+1;
                                        if(counter == 2){
                                            collidedClusters[0][0] = myArr[i][0];
                                            collidedClusters[0][1] = myArr[i][1];
                                            collidedClusters[0][2] = myArr[i][2];
                                            collidedClusters[0][3] = myArr[i][3];
                                            collidedClusters[0][4] = myArr[i][4];
                                            collidedClusters[1][0] = receivedArr[r][0];
                                            collidedClusters[1][1] = receivedArr[r][1];
                                            collidedClusters[1][2] = receivedArr[r][2];
                                            collidedClusters[1][3] = receivedArr[r][3];
                                            collidedClusters[1][4] = receivedArr[r][4];
                                        }else if(counter>2){
                                            for(int p=counter-1; p<counter; p++){
                                                collidedClusters[p][0] = receivedArr[r][0];
                                                collidedClusters[p][1] = receivedArr[r][1];
                                                collidedClusters[p][2] = receivedArr[r][2];
                                                collidedClusters[p][3] = receivedArr[r][3];
                                                collidedClusters[p][4] = receivedArr[r][4];
                                            }
                                        }
                                    }
                                }
                            }
                            //set as UNCLASSIFIED all points that exist into collided clusters
                            int count1=0, count2=0, k=0, j=0, l=0;
                            while(j<num){
                                for(int i=0; i<counter; i++){
                                    if(myLocal[j].cluster_id == (int)collidedClusters[i][4]){
                                        myLocal[j].cluster_id = UNCLASSIFIED;
                                        count1++;
                                    }
                                }
                                j++;
                            }
                            while(k<num){
                                for(int i=0; i<counter; i++){
                                    if(receivedLocal[k].cluster_id == (int)collidedClusters[i][4]){
                                        receivedLocal[k].cluster_id = UNCLASSIFIED;
                                        count2++;
                                    }
                                }
                                k++;
                            }
                            //mine all UNCLASSIFIED points from local and received arrays and save them into unclassLocal array for further manipulation
                            int total = count1+count2;
                            point_t *unclassLocal;
                            unclassLocal = (point_t *)calloc(total*4,sizeof(point_t));
                            while(l<count1){
                                for(int i=0; i<num; i++){
                                    if(myLocal[i].cluster_id == UNCLASSIFIED){
                                        unclassLocal[l].x = myLocal[i].x;
                                        unclassLocal[l].y = myLocal[i].y;
                                        unclassLocal[l].z = myLocal[i].z;
                                        unclassLocal[l].cluster_id = myLocal[i].cluster_id;
                                        l++;
                                    }
                                }
                            }
                            while(l<total){
                                for(int m=0; m<num; m++){
                                    if(receivedLocal[m].cluster_id == UNCLASSIFIED){
                                        unclassLocal[l].x = receivedLocal[m].x;
                                        unclassLocal[l].y = receivedLocal[m].y;
                                        unclassLocal[l].z = receivedLocal[m].z;
                                        unclassLocal[l].cluster_id = receivedLocal[m].cluster_id;
                                        l++;
                                    }
                                }
                            }
                            //find median coordinates of all correlated points  
                            double tempx =0, tempy=0, tempz=0;
                            for(k=0; k<total; k++){
                                tempx = tempx + unclassLocal[k].x;
                                tempy = tempy + unclassLocal[k].y;
                                tempz = tempz + unclassLocal[k].z;
                            }

                            double newCorePx = tempx/total;
                            double newCorePy = tempy/total;
                            double newCorePz = tempz/total;
                            double dist;
                            double maxdist = 0;
                            int ff = 0;
                            //measure distance and classify from scratch
                            for(k=0; k<total; k++){
                                dist = dist2(newCorePx, newCorePy, newCorePz, unclassLocal[k].x, unclassLocal[k].y, unclassLocal[k].z);
                                if(dist >= epsilon){
                                    unclassLocal[k].cluster_id = maxId+1;
                                    ff++;
                                    if(maxdist<dist){
                                        maxdist = dist;
                                    }
                                }else{
                                    unclassLocal[k].cluster_id = NOISE;
                                }
                            }
                            
                            
                            if(ff<minpts){
                                continue;
                            }else{
                                //reclassify points into mergeLocal comparing with the points of previous process of classification
                                for(k=0; k<2*num; k++){
                                    for(l=0; l<total; l++){
                                        if(mergeLocal[k].x == unclassLocal[l].x && mergeLocal[k].y == unclassLocal[l].y
                                        && mergeLocal[k].z == unclassLocal[l].z){
                                            mergeLocal[k].cluster_id = unclassLocal[l].cluster_id;
                                        }
                                    }
                                }
                                
                                for(k=0; k<total; k++){
                                    if(unclassLocal[k].cluster_id != NOISE && unclassLocal[k].cluster_id != UNCLASSIFIED){
                                        maxId = unclassLocal[k].cluster_id;
                                    }
                                }
                                
                                for(int i=0; i<newCountOfCores; i++){
                                    for(int j=0; j<counter; j++){
                                        if(newArr[i][0] == collidedClusters[j][0] && newArr[i][1] == collidedClusters[j][1]
                                        && newArr[i][2] == collidedClusters[j][2] && newArr[i][3] == collidedClusters[j][3]
                                        && newArr[i][4] == collidedClusters[j][4]){
                                            newArr[i][0] = 0;
                                            newArr[i][1] = 0;
                                            newArr[i][2] = 0;
                                            newArr[i][3] = 0;
                                            newArr[i][4] = 0;
                                        }
                                    }
                                }                    
                                
                                newArr[n+newCountOfCores][0] = newCorePx;
                                newArr[n+newCountOfCores][1] = newCorePy;
                                newArr[n+newCountOfCores][2] = newCorePz;
                                newArr[n+newCountOfCores][3] = dist;
                                newArr[n+newCountOfCores][4] = maxId;
                                n++;
                                
                                maxId = maxId + 1;
                            }
                        }
                        //rearrange newArr and myLocal and set them ready to send for further manipulation
                        int p = 0;
                        for(int i=0; i<n+newCountOfCores; i++){
                            if(newArr[i][0] == 0 && newArr[i][1] == 0 
                            && newArr[i][2] == 0 && newArr[i][3] == 0 && newArr[i][4] == 0){
                                p++;
                            }
                        }
                        
                        numOfCores = n+newCountOfCores - p;
                        myArr = calloc(numOfCores*5, sizeof(double*));
                        int f = 0;
                        for(int p=0; p<newCountOfCores+n; p++){
                            if(newArr[p][0] == 0 && newArr[p][1] == 0
                            && newArr[p][2] == 0 && newArr[p][3] == 0 && newArr[p][4] == 0){
                                continue;
                            }else{
                                myArr[f][0] = newArr[p][0];
                                myArr[f][1] = newArr[p][1];
                                myArr[f][2] = newArr[p][2];
                                myArr[f][3] = newArr[p][3];
                                myArr[f][4] = newArr[p][4];
                                f++;
                            }
                        }
                        
                        for(int i=0; i<2*num; i++){
                            myLocal[i].x = mergeLocal[i].x;
                            myLocal[i].y = mergeLocal[i].y;
                            myLocal[i].z = mergeLocal[i].z;
                            myLocal[i].cluster_id = mergeLocal[i].cluster_id;
                        }
                    }
                }
            }else{
                int senderRank, recverRank;
                int i = pow(2, (step+1)) - 1;
                int sizeAtEachStep = (int)pow(2, step);
                for(int j=0; j<size/(2 * pow(2, step)); j++){
                    recverRank = i;
                    senderRank = i - (int) pow(2, step);
                    i += pow(2, step+1);

                    if(rank == senderRank & rank != size-1){

                        MPI_Send(myLocal, sizeAtEachStep*num, mpi_point_t_type, recverRank, tag1, MPI_COMM_WORLD);
                        MPI_Send(&myArr[0][0], numOfCores*5, MPI_DOUBLE, recverRank, tag2, MPI_COMM_WORLD);

                    }else if(rank == recverRank){

                        MPI_Status status;

                        int count1;
                        MPI_Probe(senderRank, tag1, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, mpi_point_t_type, &count1);
                        receivedLocal = (point_t *)calloc(count1,sizeof(point_t));
                        MPI_Recv(receivedLocal, count1, mpi_point_t_type, senderRank, tag1, MPI_COMM_WORLD, &status);
                        
                        int count2;
                        MPI_Probe(senderRank, tag2, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_DOUBLE, &count2);
                        receivedArr = calloc(count2, sizeof(double*));
                        MPI_Recv(receivedArr, count2, MPI_DOUBLE, senderRank, tag2, MPI_COMM_WORLD, &status);

                        int maxId1 = 0;
                        for(int j=0; j<count1; j++){
                            if(myLocal[j].cluster_id > maxId1){
                                maxId1 = myLocal[j].cluster_id;
                            }
                        }
                        
                        int maxId2 = 0;
                        for(int j=0; j<count1; j++){
                            if(receivedLocal[j].cluster_id > maxId2){
                                maxId2 = receivedLocal[j].cluster_id;
                            }
                        }

                        int maxId;
                        if(maxId1>maxId2){
                            maxId = maxId1;
                            for(int j=0; j<count1; j++){
                                receivedLocal[j].cluster_id = receivedLocal[j].cluster_id + maxId1 + 1;
                            }
                        }else{
                            maxId = maxId2;
                            for(int j=0; j<count1; j++){
                                myLocal[j].cluster_id = myLocal[j].cluster_id + maxId2 + 1;
                            }
                        }

                        point_t *mergeLocal;
                        mergeLocal = (point_t *)calloc(2*count1,sizeof(point_t));
                        
                        int i=0;                        
                        while(i<2*count1){
                            if(i<count1){
                                mergeLocal[i].x  = receivedLocal[i].x;
                                mergeLocal[i].y  = receivedLocal[i].y;
                                mergeLocal[i].z  = receivedLocal[i].z;
                                mergeLocal[i].cluster_id  = receivedLocal[i].cluster_id;
                            }else{
                                mergeLocal[i].x  = myLocal[i-count1].x;
                                mergeLocal[i].y  = myLocal[i-count1].y;
                                mergeLocal[i].z  = myLocal[i-count1].z;
                                mergeLocal[i].cluster_id  = myLocal[i-count1].cluster_id;
                            }
                            i++;
                        }

                        newCountOfCores = numOfCores + count2/5;
                        newArr = calloc(newCountOfCores*5, sizeof(double*)); 
                        int j=0;
                        while(j<newCountOfCores){
                            for(int l=0; l<numOfCores; l++){
                                newArr[j][0] = myArr[l][0];
                                newArr[j][1] = myArr[l][1];
                                newArr[j][2] = myArr[l][2];
                                newArr[j][3] = myArr[l][3];
                                newArr[j][4] = myArr[l][4];
                                j++;
                            }
                            for(int l=0; l<count2/5; l++){
                                newArr[j][0] = receivedArr[l][0];
                                newArr[j][1] = receivedArr[l][1];
                                newArr[j][2] = receivedArr[l][2];
                                newArr[j][3] = receivedArr[l][3];
                                newArr[j][4] = receivedArr[l][4];
                                j++;
                            }
                        }

                        //case 1: one lies entirely within the other;
                        for(i=0; i<numOfCores; i++){
                            for(int r=0; r<count2/5; r++){
                                double d = (double)sqrt(pow(myArr[i][0]-receivedArr[r][0], 2)+pow(myArr[i][1]-receivedArr[r][1], 2)+pow(myArr[i][2]-receivedArr[r][2], 2));
                                double radDif = (double)abs(myArr[i][3] - receivedArr[r][3]);
                                
                                if(d<=radDif){
                                    if(myArr[i][3]>receivedArr[r][3]){
                                        for(int j=0; j<count1; j++){
                                            if(receivedLocal[j].cluster_id == (int)receivedArr[r][4]){
                                                receivedLocal[j].cluster_id = (int)myArr[i][4];
                                            }
                                        }
                                        for(int i=0; i<newCountOfCores; i++){
                                            for(int j=0; j<5; j++){
                                                for(int k=0; k<count2/5; k++){
                                                    for(int p=0; p<5; p++){
                                                        if(newArr[i][j] == receivedArr[k][p]){
                                                            newArr[i][j] = 0;
                                                        }
                                                    }
                                                }
                                            }
                                        }                                            
                                    }else{
                                        for(int j=0; j<count1; j++){
                                            if(myLocal[j].cluster_id == myArr[i][4]){
                                                myLocal[j].cluster_id = (int)receivedArr[r][4];
                                            }
                                        }
                                        for(int i=0; i<newCountOfCores; i++){
                                            for(int j=0; j<5; j++){
                                                for(int k=0; k<numOfCores; k++){
                                                    for(int p=0; p<5; p++){
                                                        if(newArr[i][j] == myArr[k][p]){
                                                            newArr[i][j] = 0;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    for(int j=0; j<2*count1; j++){
                                        for(int l=0; l<count1; l++){
                                            if(mergeLocal[j].x == myLocal[l].x && mergeLocal[j].y == myLocal[l].y && mergeLocal[j].z == myLocal[l].z){
                                                mergeLocal[j].cluster_id = myLocal[l].cluster_id;
                                            }
                                        }
                                        for(int l=0; l<count1; l++){
                                            if(mergeLocal[j].x == receivedLocal[l].x && mergeLocal[j].y == receivedLocal[l].y && mergeLocal[j].z == receivedLocal[l].z){
                                                mergeLocal[j].cluster_id = receivedLocal[l].cluster_id;
                                            }
                                        }
                                    } 
                                }
                            }
                        }
                        int n=1;
                        for(int i=0; i<numOfCores; i++){
                            int counter = 1;
                            double (*collidedClusters)[5];
                            collidedClusters = calloc(50, sizeof(double*));
                            for(int r=0; r<count2/5; r++){
                                //case 2: intersect in two points;
                                double d = (double)sqrt(pow(myArr[i][0]-receivedArr[r][0], 2)+pow(myArr[i][1]-receivedArr[r][1], 2)+pow(myArr[i][2]-receivedArr[r][2], 2));
                                double radConcat = myArr[i][3] + receivedArr[r][3];
                                if(d<radConcat){
                                    int count3 = 0;
                                    int count4 = 0;
                                    int m = 0;
                                    int k = 0;
                                    int l = 0;
                                    while(m<count1){
                                        if(myLocal[m].cluster_id == myArr[i][4]){
                                            count3++;
                                        }
                                        m++;
                                    }
                                    while(k<count1){
                                        if(receivedLocal[k].cluster_id == (int)receivedArr[r][4]){
                                            count4++;
                                        }
                                        k++;
                                    }

                                    if(count3 == 0 || count4 == 0){
                                        continue;
                                    }else{
                                        counter = counter+1;
                                        if(counter == 2){
                                            collidedClusters[0][0] = myArr[i][0];
                                            collidedClusters[0][1] = myArr[i][1];
                                            collidedClusters[0][2] = myArr[i][2];
                                            collidedClusters[0][3] = myArr[i][3];
                                            collidedClusters[0][4] = myArr[i][4];
                                            collidedClusters[1][0] = receivedArr[r][0];
                                            collidedClusters[1][1] = receivedArr[r][1];
                                            collidedClusters[1][2] = receivedArr[r][2];
                                            collidedClusters[1][3] = receivedArr[r][3];
                                            collidedClusters[1][4] = receivedArr[r][4];
                                        }else if(counter>2){
                                            for(int p=counter-1; p<counter; p++){
                                                collidedClusters[p][0] = receivedArr[r][0];
                                                collidedClusters[p][1] = receivedArr[r][1];
                                                collidedClusters[p][2] = receivedArr[r][2];
                                                collidedClusters[p][3] = receivedArr[r][3];
                                                collidedClusters[p][4] = receivedArr[r][4];
                                            }
                                        }
                                    }
                                }
                            }
                            
                            int count3=0, count4=0, k=0, m=0, l=0;
                            while(m<count1){
                                for(int i=0; i<counter; i++){
                                    if(myLocal[m].cluster_id == (int)collidedClusters[i][4]){
                                        myLocal[m].cluster_id = UNCLASSIFIED;
                                        count3++;
                                    }
                                }
                                m++;
                            }
                            while(k<count1){
                                for(int i=0; i<counter; i++){
                                    if(receivedLocal[k].cluster_id == (int)collidedClusters[i][4]){
                                        receivedLocal[k].cluster_id = UNCLASSIFIED;
                                        count4++;
                                    }
                                }
                                k++;
                            }
                            int total = count3+count4;
                            point_t *unclassLocal;
                            unclassLocal = (point_t *)calloc(total*4,sizeof(point_t));
                            while(l<count3){
                                for(int i=0; i<count1; i++){
                                    if(myLocal[i].cluster_id == UNCLASSIFIED){
                                        unclassLocal[l].x = myLocal[i].x;
                                        unclassLocal[l].y = myLocal[i].y;
                                        unclassLocal[l].z = myLocal[i].z;
                                        unclassLocal[l].cluster_id = myLocal[i].cluster_id;
                                        l++;
                                    }
                                }
                            }
                            while(l<total){
                                for(int m=0; m<count1; m++){
                                    if(receivedLocal[m].cluster_id == UNCLASSIFIED){
                                        unclassLocal[l].x = receivedLocal[m].x;
                                        unclassLocal[l].y = receivedLocal[m].y;
                                        unclassLocal[l].z = receivedLocal[m].z;
                                        unclassLocal[l].cluster_id = receivedLocal[m].cluster_id;
                                        l++;
                                    }
                                }
                            }
                                
                            double tempx =0, tempy=0, tempz=0;
                            for(k=0; k<total; k++){
                                tempx = tempx + unclassLocal[k].x;
                                tempy = tempy + unclassLocal[k].y;
                                tempz = tempz + unclassLocal[k].z;
                            }

                            double newCorePx = tempx/total;
                            double newCorePy = tempy/total;
                            double newCorePz = tempz/total;
                            double dist;
                            double maxdist = 0;
                            int ff = 0;
                            for(k=0; k<total; k++){
                                dist = dist2(newCorePx, newCorePy, newCorePz, unclassLocal[k].x, unclassLocal[k].y, unclassLocal[k].z);
                                if(dist >= epsilon){
                                    unclassLocal[k].cluster_id = maxId+1;
                                    ff++;
                                    if(maxdist<dist){
                                        maxdist = dist;
                                    }
                                }else{
                                    unclassLocal[k].cluster_id = NOISE;
                                }
                            }

                            if(ff<minpts){
                                continue;
                            }else{
                                for(k=0; k<2*count1; k++){
                                    for(l=0; l<total; l++){
                                        if(mergeLocal[k].x == unclassLocal[l].x && mergeLocal[k].y == unclassLocal[l].y
                                        && mergeLocal[k].z == unclassLocal[l].z){
                                            mergeLocal[k].cluster_id = unclassLocal[l].cluster_id;
                                        }
                                    }
                                }
                                
                                for(k=0; k<total; k++){
                                    if(unclassLocal[k].cluster_id != NOISE && unclassLocal[k].cluster_id != UNCLASSIFIED){
                                        maxId = unclassLocal[k].cluster_id;
                                    }
                                }
                                
                                for(int i=0; i<newCountOfCores; i++){
                                    for(int j=0; j<counter; j++){
                                        if(newArr[i][0] == collidedClusters[j][0] && newArr[i][1] == collidedClusters[j][1]
                                        && newArr[i][2] == collidedClusters[j][2] && newArr[i][3] == collidedClusters[j][3]
                                        && newArr[i][4] == collidedClusters[j][4]){
                                            newArr[i][0] = 0;
                                            newArr[i][1] = 0;
                                            newArr[i][2] = 0;
                                            newArr[i][3] = 0;
                                            newArr[i][4] = 0;
                                        }
                                    }
                                }                    
                                
                                newArr[n+newCountOfCores][0] = newCorePx;
                                newArr[n+newCountOfCores][1] = newCorePy;
                                newArr[n+newCountOfCores][2] = newCorePz;
                                newArr[n+newCountOfCores][3] = dist;
                                newArr[n+newCountOfCores][4] = maxId;
                                n++;
                                
                                maxId = maxId + 1;
                            }
                        }
                        int p = 0;
                        for(int i=0; i<n+newCountOfCores; i++){
                            if(newArr[i][0] == 0 && newArr[i][1] == 0 
                            && newArr[i][2] == 0 && newArr[i][3] == 0 && newArr[i][4] == 0){
                                p++;
                            }
                        }
                        
                        numOfCores = n+newCountOfCores - p;
                        myArr = calloc(numOfCores*5, sizeof(double*));
                        int f = 0;
                        for(int p=0; p<newCountOfCores+n; p++){
                            if(newArr[p][0] == 0 && newArr[p][1] == 0
                            && newArr[p][2] == 0 && newArr[p][3] == 0 && newArr[p][4] == 0){
                                continue;
                            }else{
                                myArr[f][0] = newArr[p][0];
                                myArr[f][1] = newArr[p][1];
                                myArr[f][2] = newArr[p][2];
                                myArr[f][3] = newArr[p][3];
                                myArr[f][4] = newArr[p][4];
                                f++;
                            }
                        }
                        
                        for(int i=0; i<2*count1; i++){
                            myLocal[i].x = mergeLocal[i].x;
                            myLocal[i].y = mergeLocal[i].y;
                            myLocal[i].z = mergeLocal[i].z;
                            myLocal[i].cluster_id = mergeLocal[i].cluster_id;
                        }
                    }   
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        if(rank == size-1){
            int maxID = 0;
            int minID = 0;
            int countr = 0;
            int classValue = 0;
            
            for(int i=0; i<num_points; i++){
                if(maxID<myLocal[i].cluster_id){
                    maxID = myLocal[i].cluster_id;
                }
            }
            minID = maxID;
            for(int i=0; i<num_points; i++){
                if(minID>myLocal[i].cluster_id && myLocal[i].cluster_id > -1){
                    minID = myLocal[i].cluster_id;
                }
            }
            while(minID<=maxID){
                for(int i=0; i<num_points; i++){
                    if(minID == myLocal[i].cluster_id){
                        myLocal[i].cluster_id = classValue;
                        countr++;
                    }
                }
                if(countr>0){
                    classValue++;
                }
                minID++;
                countr = 0;                
            }
            //printf("Epsilon: %lf\n", epsilon);
            //printf("Minimum points: %u\n", minpts);
            //print_points(myLocal, num_points);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0){
        tStop = clock();
        double t = (tStop-tStart)/CLOCKS_PER_SEC;
        printf("t = %f\n", tStop-tStart);
    }
    MPI_Finalize();
    return 0;
}