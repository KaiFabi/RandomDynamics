#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <float.h>

/***************/
/* Definitions */
/***************/
#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )

#define DIM 2
#define ETA 0.99
#define DIST 0.05

#define HYPOT2(x, y) sqrt(sqr(x) + sqr(y))
#define sqr(x) ((x)*(x))

#define TRUE 1
#define FALSE 0

#define POWDIM 4

typedef enum{particle, pseudoParticle} nodetype;

/**************/
/* Structures */
/**************/
typedef struct{
    double x[DIM];
    double v[DIM];
    double vtmp[DIM];
    int moved;
    int todelete;
} Particle;

typedef struct{
    double lower[DIM];
    double upper[DIM];
} Box;

typedef struct TreeNode{
    Particle p;
    Box box;
    nodetype node;
    struct TreeNode *son[POWDIM];
} TreeNode;

/**************/
/* Prototypes */
/**************/
void outputResults(TreeNode *root, FILE *fp);
void freeTree(TreeNode *root);
void moveParticles(TreeNode *root);
void setFlags(TreeNode *t);
void moveLeaf(TreeNode *t, TreeNode *root);
int outsideBox(TreeNode *t);
void repairTree(TreeNode *t);
void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box);
void initData(TreeNode **root, Box *domain, int N);
void inputParameters(double *dt, double *Tmax, Box *box, int *N);
int sonNumber(Box *box, Box *sonbox, Particle *p);
void insertTree(Particle *p, TreeNode *t);
void genDist(Particle *p, int n);
void updateX(Particle *p, double dt);
void compX(TreeNode *root, TreeNode *t);
double urand(double, double);

/********/
/* Main */
/********/
int main(int argc, char *argv[]){
    TreeNode *root;
    Box box;
    double dt, Tmax;
    int N;
    inputParameters(&dt, &Tmax, &box, &N);
    initData(&root, &box, N);
    timeIntegration(0, dt, Tmax, root, box);
    freeTree(root);
    return 0;
}

/*************/
/* Functions */
/*************/

void outputResults(TreeNode *t, FILE *fp){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            outputResults(t->son[i], fp);
        }
        // Operation on *t
        if(t->node == particle){ 
            fprintf(fp, "%lf %lf %lf %lf\n",t->p.x[0],t->p.x[1],t->p.v[0],t->p.v[1]);
        }
    }
}

void moveParticles(TreeNode *root){
    setFlags(root);
    moveLeaf(root, root);
    repairTree(root);
}

void setFlags(TreeNode *t){
	if(t != NULL){
		for(int i=0; i<POWDIM; i++){
			setFlags(t->son[i]);
        }
        // Operation on *t
		t->p.moved = FALSE;
		t->p.todelete = FALSE;
	}
}

// Check if particle is still in the box
int outsideBox(TreeNode *t){ 
    for(int d=0; d<DIM; d++){
        if((t->p.x[d] < t->box.lower[d]) || (t->p.x[d] > t->box.upper[d])){
            return TRUE; 
        }
    }
    return FALSE;
}

// Move particels in tree if they are outside of their boxes
void moveLeaf(TreeNode *t, TreeNode *root){
	if(t != NULL){
		for(int i=0; i<POWDIM; i++){
			moveLeaf(t->son[i], root);
        }
        // Operation on *t
		if((t->node == particle)&&(t->p.moved == FALSE)){
			t->p.moved = TRUE;
			if(outsideBox(t) == TRUE){
                insertTree(&t->p, root);
				t->p.todelete = TRUE;
			}
		}
	}
}

// Repair tree
void repairTree(TreeNode *t){
	if(t != NULL){
		for(int i=0; i<POWDIM; i++){
			repairTree(t->son[i]);
        }
        // Operation on *t
		if(t->node != particle){
			int numberofsons = 0;
    		int d = 0;
			for(int i=0; i<POWDIM; i++) {
				if(t->son[i] != NULL) {
					if(t->son[i]->p.todelete == TRUE){
                        free(t->son[i]);
                        t->son[i] = NULL;
                    }
                    else {
			            numberofsons++;
					    d=i;
				    }
                }
			}

			if(0 == numberofsons) {
				t->p.todelete = TRUE;
			} else if (1 == numberofsons) {
				t->p = t->son[d]->p;
				//free(t->son[d]);
                //t->son[d] = NULL;
			}
		}
	}
}

void freeTree(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            freeTree(t->son[i]);
        }
        // Operation on *t
        free(t);
    }
}

// Update position
void updateX(Particle *p, double dt){
    for(int d=0; d<DIM; d++){
        p->x[d] += dt*(p->v[d]);
    }
}

// Sum up all neighbour particles' velocity
void adaptation(Particle *i, Particle *j){
    for(int d=0; d<DIM; d++){
        i->vtmp[d] += j->v[d];
    }
}

double distanceDomainBoxParticle2D(TreeNode *tl, TreeNode *t){
    double x_min = t->box.lower[0];
    double x_max = t->box.upper[0];
    double y_min = t->box.lower[1];
    double y_max = t->box.upper[1];
    double x = tl->p.x[0];
    double y = tl->p.x[1];

    if (x < x_min) {
        if (y < y_min)
            return HYPOT2(x_min-x, y_min-y);
        else if (y <= y_max)
            return x_min - x;
        else
            return HYPOT2(x_min-x, y_max-y);
    } else if (x <= x_max) {
        if (y < y_min)
            return y_min - y;
        else if (y <= y_max)
            return 0;
        else
            return y - y_max;
    } else {
        if (y < y_min)
            return HYPOT2(x_max-x, y_min-y);
        else if (y <= y_max)
            return x - x_max;
        else
            return HYPOT2(x_max-x, y_max-y);
    }
}

void compB(TreeNode *tl, TreeNode *t, int level){
    if((tl != t) && (t != NULL)){
        if((distanceDomainBoxParticle2D(tl, t) < DIST) || (level == 0)){ // check if distance to potential neighbour particles is small enough
            if(t->node == particle){
                double r = 0;
                for(int d=0; d<DIM; d++){
                    r += sqr(t->p.x[d] - tl->p.x[d]);
                }
                r = sqrt(r); 
                if(r < DIST){
                    adaptation(&tl->p, &t->p);
                }
            }else{
                for(int i=0; i<POWDIM; i++){
                    compB(tl, t->son[i], level+1);
                }
            }

        }
    }
}


void compX(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compX(t->son[i], root);
        }
        // Operationen on *t
        if(t->node == particle){
            for(int d=0; d<DIM; d++){
                t->p.vtmp[d] = 0;
            }
            compB(t, root, 0); // "t" represents current leaf; get all neighbour leaves from current position
        }
    }
}


void compR(TreeNode *t, double dt){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compR(t->son[i], dt);
        }
        // Operation on *t
        if(t->node == particle){
            updateX(&t->p, dt);
        }
    }
}

void compY(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compY(t->son[i]);
        }
        // Operation on *t
        if(t->node == particle){
            double vabs = 0;                //
            for(int d=0; d<DIM; d++){       //
                vabs += sqr(t->p.vtmp[d]);  // compute magnitude of velocity
            }                               //
            vabs = sqrt(vabs);              // 
            if(vabs > 0){
                // compute new velocity of particle
                for(int d=0; d<DIM; d++){
                    t->p.v[d] = ETA * t->p.v[d] + (1.-ETA) * t->p.vtmp[d] / vabs;
                }
            }
        }
    }
}


void get_file_name(int k, char* buffer, size_t buflen){
    snprintf(buffer, buflen, "logs_test/data_%d.log", k);
}

void outputResults2File(TreeNode *root, int n){
    const size_t BUFLEN = 50;
    char file_name[BUFLEN];
    get_file_name(n, file_name, BUFLEN);
    FILE *fp;
    fp = fopen(file_name, "w+");
    outputResults(root, fp);
    fclose(fp);
}

void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box){
    int n=0;
    while (t<Tmax){
        t += dt;
        outputResults2File(root, n);    // print some results
        compX(root, root);              // get all neighbour particles' velocities
        compY(root);            // compute new velocity 
        compR(root, dt);        // compute new position
        moveParticles(root);    // move particles in tree
        n++;
    }
}

void initData(TreeNode **root, Box *domain, int N){
    // Allocate memory for all particles
    Particle *p = (Particle*)malloc(N*sizeof(*p));

    // Generate data
    genDist(p, N);

    // Allocate memory for tree root node
    *root = (TreeNode*)calloc(1, sizeof(TreeNode));

    // Add particles to tree
    (*root)->p = p[0];
    (*root)->box = *domain; 
    for(int i=1; i<N; i++){
        insertTree(&p[i], *root);
    }
    free(p);
}

void inputParameters(double *dt, double *Tmax, Box *box, int *N){
    *dt = 0.0001;
    *Tmax = 0.01;
    *N = 20000; // # of agents
    // Max domain size
    for(int d=0; d<DIM; d++){
        box->lower[d] = -100; 
        box->upper[d] = 100; 
    }
}

double urand(double low, double high){
    return low+((double)rand() / (double)RAND_MAX)*(high-low);
}

void genDist(Particle *p, int n){ 
    time_t t;
    srand((unsigned) time(&t));
    double eta = 10E-4;
    for(int i=0; i<n; i++){
       p[i].x[0] = urand(-1,1);
       p[i].x[1] = urand(-1,1);       
       p[i].v[0] = eta * urand(-1,1);
       p[i].v[1] = eta * urand(-1,1);
   }
}

int sonNumber(Box *box, Box *sonbox, Particle *p){
    int b=0;
    for(int d=DIM-1; d>=0; d--){
        if(p->x[d] < .5*(box->upper[d] + box->lower[d])){
            b = 2*b;
            sonbox->lower[d] = box->lower[d];
            sonbox->upper[d] = .5*(box->upper[d] + box->lower[d]);
        }else{
            b = 2*b+1;
            sonbox->lower[d] = .5*(box->upper[d] + box->lower[d]);
            sonbox->upper[d] = box->upper[d];
        }
    } 
    return b;
}

void insertTree(Particle *p, TreeNode *t){
    Box sonbox;
    int b=sonNumber(&t->box, &sonbox, p);
    
    if(t->son[b] == NULL){
        if(t->node == particle){
            Particle p2 = t->p;
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
            t->node = pseudoParticle;

            insertTree(&p2, t);
        }else{
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
        }
    }else{
        insertTree(p, t->son[b]);
    }
}

