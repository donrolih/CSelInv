#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MESHMIN 3

#define mesh(i,j)      mesh[nx*((j)-1)+(i)-1]
#define meshtop(i,j)   mesh[nxtop*((j)-1)+(i)-1]
#define meshbot(i,j)   mesh[nxbot*((j)-1)+(i)-1]
#define meshleft(i,j)  mesh[nx*((j)-1)+(i)-1]
#define meshright(i,j) mesh[nx*((j)-1)+(i)-1]
#define perm(i)        perm[(i)-1]
#define ptop(i)        ptop[(i)-1]
#define pbot(i)        pbot[(i)-1]
#define pleft(i)       pleft[(i)-1]
#define pright(i)      pright[(i)-1]

int nd2d(int nx, int ny, int *mesh, int *perm)
{
   int status = 0;
   int i, j, count;
   int *meshtop, *meshbot, *meshleft, *meshright;
   int *ptop, *pbot, *pleft, *pright; 
   int nxtop, nxbot, nyleft, nyright;
   int midcol, midrow;
   int ierr;

   if (nx <= MESHMIN && ny <= MESHMIN) {
      /* at a leaf, walk thru in column major */
      count = 0;
      for (j = 1; j <=ny; j++) 
         for (i=1; i<=nx; i++) {
            perm[count] = mesh(i,j);
            count++;
         }  
   } 
   else if ( nx > ny) {
      /* cut mesh in the middle row */
      midrow = (int)ceil((double)nx/2.0);

      /* order the top half */
      nxtop = midrow-1;
      meshtop = (int*)malloc((nxtop*ny)*sizeof(int)); 
      ptop = (int*)malloc(nxtop*ny*sizeof(int));
      for (j = 1; j<=ny; j++) 
         for (i=1; i<=nxtop; i++) {
	    /* meshtop(i,j) = mesh(i,j); */
	    meshtop[nxtop*((j)-1)+(i)-1] = mesh(i,j); 
         }

      ierr = nd2d(nxtop, ny, meshtop, ptop);
      for (i = 0; i<nxtop*ny; i++)
         perm[i] = ptop[i];

      free(ptop);
      free(meshtop);

      /* order the bottom half */
      nxbot = nx - nxtop - 1;
      meshbot = (int*)malloc(nxbot*ny*sizeof(int)); 
      for (j = 1; j<=ny; j++) 
	 for (i=1; i<=nxbot; i++) { 
	   /*meshbot(i,j) = mesh(nxtop+1+i,j);*/
	   meshbot[nxbot*((j)-1)+(i)-1] = mesh(nxtop+1+i,j);
         }

      pbot = (int*)malloc(nxbot*ny*sizeof(int));
      ierr = nd2d(nxbot, ny, meshbot, pbot);
      for (i = 0; i< nxbot*ny; i++)
         perm[nxtop*ny+i] = pbot[i];

      free(pbot); 
      free(meshbot);

      /* append the seperator */
      for (j = 0; j<ny; j++)
	 perm[(nxtop+nxbot)*ny+j]=mesh(midrow,j+1); 

   }
   else {
      /* cut mesh in the middle column */
      midcol = (int)ceil((double)ny/2);
      nyleft = midcol-1;
      meshleft = (int*)malloc(nx*nyleft*sizeof(int)); 
      for (j = 1; j<=nyleft; j++) {
         for (i=1; i<=nx; i++) 
	   /* meshleft(i,j) = mesh(i,j); */
	   meshleft[(j-1)*nx+i-1] = mesh(i,j); 
      }
      pleft = (int*)malloc(nx*nyleft*sizeof(int));
      ierr = nd2d(nx, nyleft, meshleft, pleft);
      for (i = 0; i<nyleft*nx; i++)
         perm[i] = pleft[i];

      free(pleft);
      free(meshleft);

      nyright = ny - nyleft - 1;
      meshright = (int*)malloc(nx*nyright*sizeof(int)); 
      for (j = 1; j<=nyright; j++) {
         for (i=1; i<=nx; i++) 
	   /* meshright(i,j) = mesh(i,nyleft+1+j); */
	   meshright[(j-1)*nx+i-1] = mesh(i,nyleft+1+j); 
      }
      pright = (int*)malloc(nx*nyright*sizeof(int));
      ierr = nd2d(nx, nyright, meshright, pright);
      for (i = 0; i<nyright*nx; i++)
         perm[nyleft*nx+i] = pright[i];
      free(pright);
      free(meshright);

      for (i = 0; i<nx; i++)
	perm[(nyleft+nyright)*nx+i] = mesh(i+1,midcol);

   } 
   return status;
}
