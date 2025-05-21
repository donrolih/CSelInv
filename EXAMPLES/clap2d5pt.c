#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include <unistd.h>

typedef struct { double r, i; } doublecomplex;

#define mesh(i,j) mesh[nx*((j)-1)+(i)-1]
#define nzvals(i) nzvals[(i)-1]
#define rowind(i) rowind[(i)-1]
#define colptr(i) colptr[(i)-1]
#define xsol(i)   xsol[(i)-1]
#define rhs(i)    rhs[(i)-1]
#define diag(i)   diag[(i)-1]
#define diag2(i)  diag2[(i)-1]
#define rowind_inva(i) rowind_inva[(i)-1]
#define colptr_inva(i) colptr_inva[(i)-1]
#define inva(i) inva[(i)-1]

extern int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
extern int ldlt_factor__(int *, int *, int *, doublecomplex *);
extern int ldlt_solve__(int *, doublecomplex *, doublecomplex *);
extern int ldlt_free__(int *);
extern int ldlt_blkselinv__(int *, int*, int*, doublecomplex *,int *);

extern int nd2d(int nx, int ny, int *mesh, int *perm);
#ifdef TIMING
extern double getime();
#endif

int main(int argc, char ** argv)
{
   int i, j, k, nnodes,nedges, ia, nnz, ibeg, iend ;
   int nx= -1, ny = -1, count, node;
   int *mesh, *perm;
   int *rowind, *colptr;
   doublecomplex *nzvals, *rhs, *xsol, *diag, *diag2;
   int token, order=-1;
   int Lnnz;
   double t0,t1, errmax; 
   long myclktck;
   doublecomplex dval;
   double errabs;
   int pbc = 0, chkerr = 0, printa = 0;
   int ierr = 0;
   int* rowind_inva;
   int* colptr_inva;
   doublecomplex* inva;
   int dumpL = 0;
 
   ia = 1;
   while (ia < argc) {
      if ( !strncmp(argv[ia],"-nx",3) ) {
 	 nx = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-ny",3) ) {
    	 ny = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-order",6) ) {
    	 order = atoi(&argv[ia][7]);
      }
      else if ( !strncmp(argv[ia],"-pbc",4) ) {
    	 pbc = atoi(&argv[ia][5]);
         if (pbc != 0) pbc = 1;
      }
      else if ( !strncmp(argv[ia],"-chkerr",7) ) {
    	 chkerr = atoi(&argv[ia][8]);
         if (chkerr != 0) chkerr = 1;
      }
      else if ( !strncmp(argv[ia],"-printa",7) ) {
    	 printa = atoi(&argv[ia][8]);
         if (printa != 0) printa = 1;
      }
      else {
 	 fprintf(stderr, "invalide argument!\n");
         fprintf(stderr, "usage: lap2d5pt -nx=<x> -ny=<y> -order=<n> -pbc=1 -chkerr=1 -printa=1\n");
         return 1;
      }
      ia++;
   }

   if (nx == -1 && ny > 0) nx = ny;
   if (ny == -1 && nx > 0) ny = nx;

   if (nx == -1) nx = 10;
   if (ny == -1) ny = 10;

   nnodes = nx*ny;
   mesh = (int*)malloc(nnodes*sizeof(int));

   for (i = 0; i<nnodes; i++) mesh[i]=i+1;
   
   /* first pass to count the number of edges */
   if (pbc) {
      /* Periodic BC */ 
      nedges = 4*nnodes;
   }
   else {
      /* Dirichlet BC */ 
      nedges = 0; 
      for (j = 1; j <= ny; j++) {
         for (i = 1; i <= nx; i++) {
           if (j < ny) nedges++;
           if (j > 1)  nedges++;
           if (i < nx) nedges++;
           if (i > 1)  nedges++;   
         }
      }
   }
   /* print the matrix dimension and number of nonzeros */
   nnz = nedges/2 + nnodes;
   printf("%d  %d\n", nnodes, nnz);

   colptr = (int*)malloc((nnodes+1)*sizeof(int));
   rowind = (int*)malloc(nnz*sizeof(int));
   nzvals = (doublecomplex*)malloc(nnz*sizeof(doublecomplex));
   colptr(1) = 1;
   count = 0;
   node  = 0;

   if (pbc) {
      dval.r = 10.0;
      dval.i = 0.0;
      printf(" Periodic boundary condition\n");
   }
   else {
      dval.r = 4.0;
      dval.i = 0.0;
      printf(" Dirichlet boundary condition\n");
   }

   for (j = 1; j <= ny; j++) {
      for (i = 1; i <= nx; i++) {
	 /* diagonal */
         if (printa)
            printf("%d %d  %8.2e\n", mesh(i,j), mesh(i,j), dval); 

         rowind[count] = mesh(i,j);
         nzvals[count].r = dval.r;
         nzvals[count].i = dval.i;
         count++;

         /* lower */
         if (i < nx) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i+1,j), mesh(i,j));

            rowind[count] = mesh(i+1,j);
            nzvals[count].r = -1.0;
            nzvals[count].i = 0.0;
            count++;
         }

         if (pbc) {
	    /* bottom of the mesh */
            if (i == 1) {
               if (printa) 
                  printf("%d %d -1.0\n", mesh(nx,j), mesh(i,j));

               rowind[count] = mesh(nx,j);
               nzvals[count].r = -1.0;
               nzvals[count].i = 0.0;
               count++;
            }
         } 

         /* right */
         if (j < ny) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i,j+1), mesh(i,j)); 

            rowind[count] = mesh(i,j+1);
            nzvals[count].r = -1.0;
            nzvals[count].i = 0.0;
            count++;
         }        

         if (pbc) {
	    /* right end of the mesh */
            if  (j==1) {
               if (printa)
                  printf("%d %d -1.0\n", mesh(i,ny), mesh(i,j));

               rowind[count] = mesh(i,ny);
               nzvals[count].r = -1.0;
               nzvals[count].i = 0.0;
               count++;
            }
         }
         node++;
         colptr(node+1) = count+1; 
      } 
   }
   if (count != nnz) {
       printf(" count = %d, nnz = %d\n", count, nnz);  
       return 1;
   }

/*
   for (j = 0; j <= nnodes; j++) {
      ibeg = colptr[j]-1;
      iend = colptr[j+1]-2;
      for (i = ibeg; i<= iend; i++) {
         printf("row = %d, col = %d, val = %11.3e\n", rowind[i], j+1, nzvals[i]);
      }
   }
*/
   token = 0;
  
   if (order == 0) {
      perm = (int*)malloc(nnodes*sizeof(int));
      ierr = nd2d(nx,ny,mesh,perm); 
      /* for (i = 0; i < nnodes; i++)
         printf("perm(%d) = %d\n", i+1, perm[i]); */
   }

   ldlt_preprocess__(&token, &nnodes, colptr, rowind, &Lnnz, &order, perm);   
   ldlt_fact__(&token, colptr, rowind, nzvals);

   rhs = (doublecomplex*)calloc(nnodes,sizeof(doublecomplex));
   xsol = (doublecomplex*)calloc(nnodes,sizeof(doublecomplex));
   diag = (doublecomplex*)malloc(nnodes*sizeof(doublecomplex));

   if (chkerr) {
      /* extract the diagonal the slow way*/
      /* myclktck = sysconf(_SC_CLK_TCK);
         printf("_SC_CLK_TCK = %d\n",(int)myclktck); */
#ifdef TIMING
      t0 = getime();
#endif
      for (i=1; i<=nnodes; i++) {
         rhs(i).r = 1.0;
         rhs(i).i = 0.0;
         ldlt_solve__(&token, xsol, rhs);
         diag(i).r = xsol(i).r;
         diag(i).i = xsol(i).i;
         rhs(i).r = 0.0;
         rhs(i).i = 0.0;
      } 
#ifdef TIMING
      t1 = getime();
#endif
      printf(" Time for slow extraction = %11.3e\n", t1-t0);
   }


   /* fast extraction */
   colptr_inva = (int*)malloc(sizeof(int)*(nnodes+1));
   rowind_inva = (int*)malloc(sizeof(int)*Lnnz);
   inva  = (doublecomplex*) malloc(sizeof(doublecomplex)*Lnnz);
   ldlt_blkselinv__(&token, colptr_inva, rowind_inva, inva, &dumpL);

   diag2 = (doublecomplex*)calloc(nnodes,sizeof(doublecomplex));
   for(j = 1; j < nnodes+1; j++){
     for(i = colptr_inva(j); i < colptr_inva(j+1); i++){
       if(rowind_inva(i) == j){
	 diag2(j) = inva(i);
       }
     }
   }


   if(0){
     for(j = 1; j < nnodes+1; j++){
       for(i = colptr_inva(j); i < colptr_inva(j+1); i++){
	 fprintf(stdout, "%4d %4d %25.10e %25.10e\n", 
		 rowind_inva(i), j, inva(i).r, inva(i).i);
       }
     }
   }

   // for(i = 0; i < nnodes;i++){
     // printf("diag [%4d] = (%15.5f,%15.5f)\n", i, diag(i).r, diag(i).i);
     // printf("diag2[%4d] = (%15.5f,%15.5f)\n", i, diag2(i).r, diag2(i).i);
   // }

   if (chkerr) {
      errmax = 0.0;
      for (i=0; i<nnodes; i++) {
          errabs = fabs(diag(i).r-diag2(i).r) + fabs(diag(i).i-diag2(i).i);
          if ( errabs > errmax )
              errmax = fabs(diag(i).r-diag2(i).r)
		     + fabs(diag(i).i-diag2(i).i) ; 
      }
      printf(" errmax = %11.3e\n", errmax);
   }

   ldlt_free__(&token); 

   if (order == 0) free(perm);
   free(mesh);
   free(colptr); 
   free(rowind); 
   free(nzvals); 
   free(rhs); 
   free(xsol); 
   free(diag); 
   free(diag2); 
   free(colptr_inva);
   free(rowind_inva);
   free(inva);
}
