/* $Author: schaid $ */
/* $Date: 2004/02/16 21:04:09 $ */
/* $Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo_em_pin.c,v 1.13 2004/02/16 21:04:09 schaid Exp $ */
/* $Locker:  $ */
/*
 * $Log:
*/
/*
*License: 
*
*Copyright 2003 Mayo Foundation for Medical Education and Research. 
*
*This program is free software; you can redistribute it and/or modify it under the terms of 
*the GNU General Public License as published by the Free Software Foundation; either 
*version 2 of the License, or (at your option) any later version.
*
*This program is distributed in the hope that it will be useful, but WITHOUT ANY 
*WARRANTY; without even the implied warranty of MERCHANTABILITY or 
*FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
*more details.
*
*You should have received a copy of the GNU General Public License along with this 
*program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*Boston, MA 02111-1307 USA
*
*For other licensing arrangements, please contact Daniel J. Schaid.
*
*Daniel J. Schaid, Ph.D.
*Division of Biostatistics
*Harwick Building û Room 775
*Mayo Clinic
*200 First St., SW
*Rochester, MN 55905
*
*phone: 507-284-0639
*fax:      507-284-9542
*email: schaid@@mayo.edu 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <S.h> 
#include "haplo_em_pin.h"


/* Progressive insertion of loci into haplotypes with EM algorithm */


/*************** Global vars ******************************************************/

static long n_loci, *loci_used;               /* used for qsort functions         */   
static HAP **ret_hap_list;                    /* stored for later return to S+    */
static HAPUNIQUE **ret_u_hap_list;
static long ret_n_hap, ret_n_u_hap, ret_max_haps;

/**********************************************************************************/

void haplo_em_pin(
   long    *S_n_loci,             /* number of loci                                 */
   long    *n_subject,            /* number of subjects                             */
   double  *weight,               /* weight per subject                             */
   long    *geno_vec,             /* vector of genotypes,  col-major from           */
                                  /* n_subject x 2*n_loci matrix                    */
   long    *n_alleles,            /* vector of number alleles per locus,            */
				  /* length=n_loci                                  */
   long    *max_haps,             /* number of maximum haplotypes over all subjects */
				  /*  - CRITICAL for memory alloc                   */
   long    *max_iter,             /* max num. iters for each EM loop                */
   long    *loci_insert_order,    /* vector for order of insert of loci for         */
                                  /* progressive insertion; length = n_loci         */
   double  *min_prior,            /* trim haplo's with prior < min_prior            */
   double  *min_posterior,        /* trim subject's pair of haplos if               */
				  /* post < min_posterior                           */
   double  *tol,                  /* convergence tolerance for change in lnlike in  */ 
				  /*  EM loop                                       */
   long    *insert_batch_size,    /* number of loci to insert in a batch before     */
                                  /*  each EM loop; order of inserted               */
                                  /* loci determined by loci_insert_order           */
   long    *converge,             /* convergence indicator for EM                   */
   double  *S_lnlike,             /* lnlike from final EM                           */
   long    *S_n_u_hap,            /* number of unique haplotypes                    */
   long    *n_hap_pairs,          /* total number of pairs of haplotypes over all   */
                                  /* subjects                                       */
   long    *random_start,         /* indicator of random posteriors should be       */
				  /* generated at the start of each EM loop.        */
                                  /* 1 = Yes, 0 = No                                */
   long    *iseed1,               /* seeds for AS183 random unif                    */
   long    *iseed2,	
   long    *iseed3,
   long    *verbose)              /* indicator if verbose pringing during  run,     */
                                  /* for debugging. verbose=0 means no printing     */
                                  /* verbose=1 means lots of printing to screen     */
{


  long i, j, k, iter, n_iter, insert_loc;
  long is, ie, n_batch;
  long n_u_hap, n_hap, n_trim, pair_id, len_hap_list;
  long **geno;
  double lnlike, lnlike_old;
  double *prior;

  HAP **hap_list;     /* List of all haplotypes = array of pointers to hap structs */
  HAPUNIQUE **u_hap_list;   /* List of unique haplotypes */
  HAP *h1, *h2;


  /* convert from S vecs to  C structures */

  n_loci = *S_n_loci;
  geno = long_vec_to_mat(geno_vec, *n_subject, 2*n_loci);


  if(*verbose){
    printf("geno matrix:\n");
    for(i=0;i< *n_subject;i++){
      for (j=0; j< (2*n_loci); j++) {
        printf("%i ",geno[i][j]);
      }
      printf("\n");
    }
  }

  prior = (double *) Calloc(*max_haps, double);
 
  if(prior==NULL){
    errmsg("could not alloc mem for prior");
  }

  /* array to keep track of loci used at any point in time */
  loci_used = (long *) Calloc(n_loci, long);
  if(loci_used==NULL){
    errmsg("could not alloc mem for loci_used");
  }

  /* array of pointers to haplo information */
  hap_list = (HAP **) Calloc(*max_haps, HAP* );
  if(hap_list==NULL){
    errmsg("could not alloc mem for hap_list");
  }


  /* put geno data into haplo list */
  
  pair_id = - 1;
  n_hap=0;
  for(i=0;i< *n_subject;i++){

    pair_id++;

    h1 = new_hap(i, pair_id, weight[i], 0.0, 1.0);
    if(!h1){
      errmsg("could not alloc mem for new_hap");
    }

 
    h2 = new_hap(i, pair_id, weight[i], 0.0, 1.0);
    if(!h2){
      errmsg("could not alloc mem for new hap");
    }


    k=0;
    for (j=0; j< n_loci; j++) {
	  h1 ->loci[j] = geno[i][k++];
	  h2 ->loci[j] = geno[i][k++];
    }

    hap_list[n_hap++] = h1;
    hap_list[n_hap++] = h2;

  }

  /* set seed for random unif generator, if needed */

  if(*random_start){
    ranAS183_seed(*iseed1, *iseed2, *iseed3);
  }



  /* begin haplotype algorithm */

  if(*verbose){

    printf("min_posterior = %8.5f\n",*min_posterior);
    printf("min_prior     = %8.5f\n",*min_prior);

    printf("loci_insert_order = ");
    for(i=0;i<n_loci;i++){
      printf("%i ",loci_insert_order[i]);
    }
    printf("\n\n");
  }


  is=0;
  ie= imin(*insert_batch_size, n_loci);
  n_batch=0;

  do {
    n_batch++;

    if(*verbose){
      printf("Inserting batch %i, loci= ",n_batch);
        for(i=is; i < ie ;i++){
          printf("%i ",loci_insert_order[i]);
         }
         printf("\n");
     }

     /* sort according to subj id, hap pair id, before insert & expand */
     qsort(hap_list, n_hap, sizeof(HAP *), cmp_subId_hapPairId);

     /* insert batch of loci */
     for(i=is; i < ie ;i++){
       insert_loc = loci_insert_order[i];
       n_hap = hap_enum(&hap_list, &prior, max_haps, n_alleles, insert_loc, n_hap, &pair_id);
     }


    /* sort according to subj id, hap pair id, after insert & expand */
    qsort(hap_list, n_hap, sizeof(HAP *), cmp_subId_hapPairId);

 
    /* set post for newly expanded haplos */
    set_posterior(n_hap, hap_list, random_start);

  
    if(*verbose){
      printf("\nhap_list after insert batch %i & set_post, before code haplo\n\n",n_batch);
      write_hap_list(hap_list, n_hap);
    }

    /* sort according to haplotype order - needed for computing unique haplos and their codes */
    qsort(hap_list, n_hap, sizeof(HAP *), cmp_hap);

 
    /* compute hap codes when computing n_u_hap */
    n_u_hap= code_haps(n_hap, hap_list);

  
    /* last sort before EM, according to subj id, then hap pair id */
     qsort(hap_list, n_hap, sizeof(HAP *), cmp_subId_hapPairId);

    if(*verbose){
      printf("\nhap_list after code haplo, before EM\n\n",n_batch);
      write_hap_list(hap_list, n_hap);
    }


    /* begin EM iterations */

    n_iter = 0;
    lnlike = 0.0;
    (*converge) = 0;

    for(iter = 0; iter< (*max_iter); iter++){

      n_iter++;

      hap_prior(n_hap, hap_list, prior, n_u_hap, *min_prior); 

      n_trim = hap_posterior(n_hap, hap_list, prior, n_u_hap, *min_posterior, &lnlike);

      if(*verbose){
        printf("\nprior probabilities\n\n");
        write_prior(n_u_hap, prior);
        printf("\nhap_list after compute posterior (n_trim = %ld))\n\n",n_trim);
        write_hap_list(hap_list, n_hap);
        printf("     iter = %3i, max_iter=%3i, lnlike = %f\n",iter, *max_iter, lnlike);
      }

       /* check for convergence */
       if(iter==0) {
          lnlike_old = lnlike;
          continue;
       } 
        else {
          if (fabs(lnlike - lnlike_old) < *tol){
            (*converge) = 1;
             break;
          } 
          else {
	    lnlike_old = lnlike;
	  }
	}

     } /* end of EM loop */

      if( (*converge)==0){
        PROBLEM "failed to converge for batch %i after %i iterations", n_batch, n_iter
		RECOVER(NULL_ENTRY);
      }

       divideKeep(hap_list, n_hap, &len_hap_list);
       n_hap = len_hap_list;

       if(*verbose){
	 printf("\nhap_list after EM and after divideKeep \n\n",n_trim);
	 write_hap_list(hap_list, n_hap);
       }
     
   
      /* update priors, in case haplos were trimmed during posteior calculations */

      hap_prior(n_hap, hap_list, prior, n_u_hap, *min_prior); 

      if(*verbose){
        if( (*converge)==1){
          printf("\n\nConverged after batch insertion, lnlike = %8.5f, n_iter = %i\n\n", lnlike, n_iter);
        }
      }

      is = ie;
      ie = is + imin(*insert_batch_size, (n_loci-is));

  } while(is < n_loci); /* end inserting all loci */


  /* copy values to return to S */
 
  *S_lnlike = lnlike;

  *n_hap_pairs = n_hap/2;

  /* Because some haplotypes may have been trimmed, we need to determine the number
     of unique haplotypes remaining after trimming */
 
  /* sort according to haplotype code - needed for computing unique haplos */
 
  qsort(hap_list, n_hap, sizeof(HAP *), cmp_hap_code);
 
  n_u_hap = count_unique_haps(n_hap, hap_list);
  *S_n_u_hap = n_u_hap;

  
 /* prepare to return info for unique haplotypes */
  u_hap_list = (HAPUNIQUE **) Calloc(n_u_hap, HAPUNIQUE *);
  if (!u_hap_list){
    errmsg("could not alloc mem for unique haplo");
  }

  unique_haps(n_hap, hap_list, u_hap_list, prior);

  if(*verbose){
    printf("\nn_u_hap = %ld\n",n_u_hap);
    printf("\nunique haps\n\n"); 
    write_unique_hap_list(u_hap_list, n_u_hap); 
  }

  ret_u_hap_list = u_hap_list;

  /* prepare to return info for subjects */

  qsort(hap_list, n_hap, sizeof(HAP *), cmp_subId_hapPairId);
  ret_hap_list = hap_list;

  if(*verbose){
    printf("ret_hap_list\n");
    write_hap_list(ret_hap_list, n_hap);
  }


  /* the following ret (return) values are used for array sizes, for 
     copying data into S arrays, and for later freeing memory */

  ret_n_hap=n_hap;
  ret_n_u_hap = n_u_hap;
  ret_max_haps = (*max_haps);

  /* Free memory */

  Free(prior);
  prior = NULL;

  Free(loci_used);
  loci_used = NULL;

  for(i=0; i< (*n_subject); i++){
    Free(geno[i]);
  }
  Free(geno);
  geno = NULL;


}

/***********************************************************************************/

static HAP* new_hap(long id, long pair_id, double wt, double prior, double post){
  HAP *result;
  long *loc;

  result = (HAP *) Calloc(1, HAP);
 
 if (!result){
  errmsg("could not alloc mem for new hap");
  }

  result->id = id;
  result->pair_id = pair_id;
  result->wt = wt;
  result->post  = post;
  result->keep = 1;

  loc = (long *) Calloc(n_loci, long);
  if (!loc) {
    errmsg("could not alloc mem for new hap");
    Free(result);
  }

  result->loci = loc; 

  return result;
}

/***********************************************************************************/

static void write_hap_list(HAP** so, long n_hap){
  long i,j;

  printf("subID     wt hapPairID hapCode keep");
  for(i=0;i<n_loci;i++){
     if(loci_used[i]==0) continue; 
    printf(" L%2ld",i);
  }
  printf("    post\n");

  for(i=0; i< n_hap ;i++){
    printf("%5ld %6.4f %9ld %7ld %4i", so[i]->id, so[i]->wt, so[i]->pair_id,so[i]->code,so[i]->keep);
    for(j=0;j<n_loci;j++){
       if(loci_used[j]==0) continue; 
      printf("%4ld",so[i]->loci[j]);
    }

    printf("    %6.4f", so[i]->post);
    printf("\n");
  }
}

/***********************************************************************************/

static void write_unique_hap_list(HAPUNIQUE** so, long n_hap){
  long i,j;

  printf("hapCode keep");
  for(i=0;i<n_loci;i++){
     if(loci_used[i]==0) continue; 
    printf(" L%2ld",i);
  }
  printf("  prior\n");

  for(i=0; i< n_hap ;i++){
    printf("%6ld %4i",so[i]->code,so[i]->keep);
    for(j=0;j<n_loci;j++){
       if(loci_used[j]==0) continue; 
      printf("%4ld",so[i]->loci[j]);
    }

    printf("    %6.4f", so[i]->prior);
    printf("\n");
  }
}

/***********************************************************************************/

static void write_prior(long n, double *prior){
  long i;

  printf("hapCode  prior\n");
  for(i=0;i<n;i++){
    printf(" %5ld  %6.4f\n", i, prior[i]);
  }

}

/***********************************************************************************/

static int CDECL cmp_hap(const void *to_one, const void *to_two){
  long i;
  long *loc1, *loc2;
  long a1, a2;
  HAP *one, *two;
  one = * (HAP **) to_one;
  two = * (HAP **) to_two;
  loc1 = one->loci;
  loc2 = two->loci;
  for (i=0; i<n_loci; i++) {
    if(loci_used[i]==0) continue;
    a1 = loc1[i];
    a2 = loc2[i];
    if (a1<a2) return -1;
    if (a1>a2) return +1;
  }
  return 0;
}

/***********************************************************************************/

static int CDECL cmp_hap_code(const void *to_one, const void *to_two){
  HAP *one, *two;
  one = * (HAP **) to_one;
  two = * (HAP **) to_two;
  if((one->code)  <  (two->code)) return -1;
  if((one->code)  >  (two->code)) return  1;
  return 0;
}


/***********************************************************************************/

static int CDECL cmp_subId_hapPairId(const void *to_one, const void *to_two){
 
  /* Using this comparision function with qsort results in a sort first on subj id, 
    then on hap pair_id */

  HAP *one, *two;
  one = * (HAP **) to_one;
  two = * (HAP **) to_two;
 
  if((one->id)  < (two->id)) return -1;
  if((one->id)  > (two->id)) return  1;
  
 
  if( (one->pair_id) < (two->pair_id)) return -1;
  if( (one->pair_id) > (two->pair_id)) return 1;

  return 0;
}




/***********************************************************************************/

static void unique_haps(long n_hap, HAP **hap_list, HAPUNIQUE **u_hap_list,
                        double *prior) {

 /* assumes hap_list is sorted by either haplotype (cmp_hap)
     or haplotype code (cmp_trim) */

  HAP **hs, **he, **h;

  hs = hap_list;
  he = hap_list + n_hap;
 
  while (hs < he) {
    h = hs;
    do {
      h++;
    } while ( (h<he) &&  ((*hs)->code == (*h)->code) );
    *u_hap_list++ = copy_hap_unique(*hs, prior);
    hs = h;
  }

}

/***********************************************************************************/

static HAPUNIQUE* copy_hap_unique(HAP *old, double *prior) {
  HAPUNIQUE *result;
  long i;
  result = (HAPUNIQUE *) Calloc(1, HAPUNIQUE);
  if (result) {
    result->code    = old->code;
    result->prior   = prior[old->code];
    result->keep    = old->keep;
    result->loci = (long *) Calloc(n_loci, long);
    if (result->loci==NULL) {
      errmsg("could not alloc mem for copy_hap_unique");
      Free(result);
    }
    for (i=0; i<n_loci; i++) 
	result->loci[i] = old->loci[i];
  } 
  return result;
}


/***********************************************************************************/

static long code_haps(long n_hap, HAP **hap_list) {

  /* assumes hap_list is sorted by either haplotype (cmp_hap)
     or haplotype code (cmp_trim) */

  HAP **hs, **he, **h;
  long res = 0;
  hs = hap_list;
  he = hap_list + n_hap;
  while (hs < he) {
    h = hs;
    do {
      (*h)->code = res;
      h++;
    } while ((h<he) && (cmp_hap(hs, h)==0));
    res++;
    hs = h;
  }
  return res;
}
/***********************************************************************************/

static long count_unique_haps(long n_hap, HAP **hap_list) {

  /* assumes hap_list is sorted by either haplotype (cmp_hap)
     or haplotype code (cmp_hap_code) */

  HAP **hs, **he, **h;
  long res = 0;
  hs = hap_list;
  he = hap_list + n_hap;
  while (hs < he) {
    h = hs;
    do {
      h++;
    } while ( (h<he) && ( (*hs)->code == (*h)->code) );
    res++;
    hs = h;
  }
  return res;
}


/***********************************************************************************/

static long hap_enum(HAP ***hap_list_ptr, double **prior_ptr, long *max_haps, long *n_alleles, long insert_loc, 
                long n_hap, long *pair_id_ptr){

  long i,j, a_poss,a1_poss,a2_poss, a1, a2,a1_new,a2_new;
  long n_al, n_miss;
  HAP *h1, *h2, *h1_new, *h2_new;
 
  j = n_hap - 1;

  loci_used[insert_loc] = 1;
 
  for(i=0;i<(n_hap-1);i+=2){

    h1 = (*hap_list_ptr)[i];
    h2 = (*hap_list_ptr)[i+1];
    a1 = h1->loci[insert_loc];
    a2 = h2->loci[insert_loc];

    /* fill in missing allele values */

    n_al = n_alleles[insert_loc];
    n_miss = (a1==0) + (a2==0);
 
  switch(n_miss){
 
  case 0:
    if((a1!=a2) && (num_het(h1,h2) > 1) ){

      /* note that nhet = number of het loci that are currently in use,
         including new insert locus. Only need to consider reciprocal
         haplotype allele insertion if het at current insert locus, and
         total num of hets across all used loci > 1 */

      insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a2, a1, pair_id_ptr, &j);

 
    }
    break;

  case 1:
       /* over-write haps for first possible alleles, to fill in one
	  possible haplotype pair, and expand if needed */

      if(a1==0){
        a1_new = 1;
        a2_new = a2;
       } else {
        a1_new = a1;
        a2_new = 1;
       }
      h1 ->loci[insert_loc] = a1_new;
      h2 ->loci[insert_loc] = a2_new;
      (*hap_list_ptr)[i]   = h1;
      (*hap_list_ptr)[i+1] = h2;

     if((a1_new!=a2_new) && (num_het(h1,h2) > 1) ){

        insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a2_new, a1_new, pair_id_ptr, &j); 
   
     }
  
     /* now consider all other values of missing allele */

     for(a_poss=2;a_poss<=n_al;a_poss++){
       
        if(a1==0){
          a1_new = a_poss;
          a2_new = a2;
         } else {
          a1_new = a1;
          a2_new = a_poss;
         }

        insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a1_new, a2_new, pair_id_ptr, &j); 
       
        /* pull out the newly inserted hap pair, to be used to
           determine number of heterozous sites, to determine
           whether new hap pair needs to be inserted */
 
         h1_new = (*hap_list_ptr)[j-1];
         h2_new = (*hap_list_ptr)[j];
 
        if((a1_new!=a2_new) && (num_het(h1_new,h2_new) > 1) ){
  
           insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a2_new, a1_new, pair_id_ptr, &j); 
           
        }
     }
     break;
  
  case 2:
      /* over-write haps for first possible alleles, to fill in one
	  possible haplotype pair, and expand if needed */
    
     h1 ->loci[insert_loc] = 1;
     h2 ->loci[insert_loc] = 1;
     (*hap_list_ptr)[i]   = h1;
     (*hap_list_ptr)[i+1] = h2;

     for(a1_poss=1;a1_poss<=n_al;a1_poss++){
       for(a2_poss=a1_poss; a2_poss<=n_al;a2_poss++){

	 if( (a1_poss==1) && (a2_poss==1) ) continue; /* did this case above */
 
         a1_new = a1_poss;
         a2_new = a2_poss;
 
         insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a1_new, a2_new, pair_id_ptr, &j); 
  
         h1_new = (*hap_list_ptr)[j-1];
         h2_new = (*hap_list_ptr)[j];

         if((a1_new!=a2_new) && (num_het(h1_new,h2_new) > 1) ){

             insert_new_hap_pair(hap_list_ptr, prior_ptr, max_haps, insert_loc,
                       h1, h2, a2_new, a1_new, pair_id_ptr, &j); 

	 }
       }
     }
     break;

   default:
    errmsg("error for number missing alleles");
  }
  }

  return (j + 1); /* return value is new number of haplos after all expanded */

}

/***********************************************************************************/

static HAP* copy_hap(HAP *old) {
  HAP *result;
  long i;
  result = (HAP *) Calloc(1, HAP);
  if (result) {
    result->id      = old->id;
    result->pair_id = old->pair_id;
    result->wt      = old->wt;
    result->post    = old->post;
    result->code    = old->code;
    result->keep    = old->keep;
    result->loci = (long *) Calloc(n_loci, long);
    if (result->loci==NULL) {
      errmsg("could not alloc mem for copy_hap");
      Free(result);
    }
    for (i=0; i<n_loci; i++) 
	result->loci[i] = old->loci[i];
  } 
  return result;
}


/***********************************************************************************/

static long num_het(HAP* h1, HAP* h2){
  long i, nhet;
  nhet = 0;
  for(i=0;i<n_loci;i++){
    if( (loci_used[i]==1) && (h1->loci[i]!=h2->loci[i]) )
      nhet++;
  }
  return nhet;
}

/***********************************************************************************/

static void hap_prior(long n_hap, HAP** hap_list, double *prior, long n_u_hap,
                      double min_prior) {

  double total, a;
  long i;

  for(i=0; i<n_u_hap; i++){
    prior[i] = 0.0;
  }


  total = 0.0;
  for(i =0; i<n_hap; i++){
    a = hap_list[i]->post * hap_list[i]->wt * hap_list[i]->keep;
    total += a;
    prior[hap_list[i]->code] += a;
  }

  for(i=0;i<n_u_hap;i++){
    prior[i] = prior[i]/total;

    if(prior[i] < min_prior){
      prior[i] = 0.0;
    }

  }

 
}

/***********************************************************************************/

static long hap_posterior(long n_hap, HAP **hap_list, double *prior, 
                         long n_u_hap, double min_posterior, double *lnlike) {


  HAP **hs, **he, **hn, **h, **h2;
  long id;
  double subtotal, gp, tmp_wt;
  long keep;
  long n_trim, total_trim;

  hs = hap_list;
  he = hap_list + n_hap;
  total_trim = 0;
  (*lnlike) = 0.0;
 
  while (hs < he) {

    h = hs;
    tmp_wt = (*h)->wt;
    subtotal = 0.0;
    n_trim = 0;
 
    /* numerator of post prob */
    do {
      id = (*h)->id;
      h2 = h+1;
      gp = prior[(*h)->code] * prior[(*h2)->code] ;
      if ((*h)->code != (*h2)->code)
	gp *= 2.0;
      subtotal += gp;
      (*h)->post = (*h2)->post = gp;
      h = h2+1;
    } while ((h<he) && (((*h)->id)==id));
 
    hn = h;
    
 
    if(subtotal > 0.0){

      /* check if need to trim by post */
      for (h=hs; h<hn; h+=2) {
       keep = ((*h)->post/subtotal  < min_posterior) ? 0 : 1;
       if(keep==0){ /* trim pair of haps */
          n_trim +=2;
          subtotal -= (*h)->post ;
          (*h)->post =  0.0;
          (*h)->keep = 0;
          (*(h+1))->post = 0.0;
          (*(h+1))->keep = 0;
       }
      }
       /* rescale if new subtotal > 0.0 */
        if(subtotal > 0.0){
          for (h=hs; h<hn; h++) {
	    (*h)->post = (*h)->post/subtotal;
	  }
	}
       /* zero post and trim all if new subtotal <= 0.0 */
       else {
        for (h=hs; h<hn; h++) {
	  (*h)->post =0.0;
          (*h)->keep = 0;
	}
       }
    }
    /* if original subtotal <= 0, zero post and trim all */
    else {
      for (h=hs; h<hn; h++) {
	  (*h)->post = 0.0;
          (*h)->keep = 0;
      }
    }

    (*lnlike) += (subtotal > 0.0) ? tmp_wt * log(subtotal) : 0.0;
    total_trim += n_trim;
    hs = hn;
  }

    return total_trim;
}

/*********************************************************************************/
static void set_posterior(long n_hap, HAP **hap_list, long *random_start){
  HAP **hs, **he, **hn, **h, **h1, **h2;
  double u, subtotal, post;
  long id;

  hs = hap_list;
  he = hap_list + n_hap;


  /* fill numertators of post */
  if(! (*random_start) )
  {
     while(hs < he){
       h1=hs;
       hs++;
       h2=hs;
       (*h1)->post = (*h2)->post = 1.0;
       hs++;
     }
  } 
  else {
     while(hs < he){
       u = ranAS183();
       h1=hs;
       hs++;
       h2=hs;
       (*h1)->post = (*h2)->post = u;
       hs++;
     }
  }


  /* standardize so post sums to 1 per subject */
 
    hs = hap_list;
    he = hap_list + n_hap;

    while(hs < he){

      subtotal = 0.0;
      h = hs;

      do {
        id = (*h)->id;
        subtotal += (*h)->post;
        h += 2;
      } while ( (h<he) && ( (*h)->id == id ) );

      hn = h; /* new end for a subject */

   
      for(h = hs; h < hn; h+=2){
        post = (*h)->post/subtotal;
        (*h)->post = post;
	(*(h+1))->post = post;
      }

      hs = hn; /* new begin for next subject */

    }


    return ;
}

/*********************************************************************************/

static long **long_matrix(long nrow, long ncol){
/* allocate long matrix with subscript range m[0 ..(nrow-1)][0..(ncol-1)] */
        long i;
        long **m;

        /* allocate pointers to rows */
        m=(long **) Calloc(nrow, long *);
        if (!m) errmsg("mem alloc failure 1 in long_matrix");
  
	/* allocate vec of memory for each row */
        for(i=0;i<nrow;i++) {
          m[i]=(long *) Calloc(ncol, long);
          if(!m[i]) errmsg("mem alloc failure 2 in long_matrix");
	}

        /* return pointer to array of pointers to rows */
        return m;
}

/***********************************************************************************/

static long **long_vec_to_mat(long *Yvec, long nrow, long ncol){

   long i,j,k;
   long **Y;

   Y=long_matrix(nrow,ncol);
   k=0;
   for (j=0;j<ncol;j++){
      for (i=0;i<nrow;i++){
         Y[i][j]=Yvec[k];
         k++;
      }
   }
   return Y;
}

/***********************************************************************************/

void haplo_em_ret_info(
   long   *n_u_hap,      /* number of unique hapoltypes                           */
   long   *S_n_loci,     /* number of loci                                        */
   long   *n_pairs,      /* number of pairs of loci over all subjects             */
   double *hap_prob,     /* probabilities for unique haplotypes, length= n_u_hap  */
   long   *u_hap,        /* unique haplotype, length=n_u_hap * n_loci             */
   long   *u_hap_code,   /* code for unique haplotypes, length=n_u_hap            */
   long   *subj_id,      /* subject id = index of subject                         */
   double *post,         /* posterior probability of pair of haplotypes           */
   long   *hap1_code,    /* code for haplotype-1 of a pair, length=n_pairs        */
   long   *hap2_code     /* code for haplotype-2 of a pair, length=n_pairs        */
  )
{


  long i,j,k;
  HAP **h;
  k= -1;
  for(i=0;i<*n_u_hap;i++){
    hap_prob[i] = ret_u_hap_list[i]->prior;
    u_hap_code[i] = ret_u_hap_list[i]->code;
    for(j=0;j<*S_n_loci;j++){
      k++;
      u_hap[k] = ret_u_hap_list[i]->loci[j];
    }
  }

  h = ret_hap_list;
  for(i=0; i<*n_pairs; i++){
    subj_id[i] = (*h)->id;
    post[i] = (*h)->post;
    hap1_code[i] = (*h)->code;
    h++;
    hap2_code[i] = (*h)->code;
    h++;
  }


  return;
}

/***********************************************************************************/

void haplo_free_memory(void){

  /* free memory saved for returned info */

  long i;


  for(i=0;i<ret_max_haps;i++){
    if(ret_hap_list[i] != NULL) {
      if(ret_hap_list[i]->loci != NULL) Free( ret_hap_list[i]->loci );
      Free( ret_hap_list[i]); 
    }
  }


  Free(ret_hap_list);

  ret_hap_list = NULL;

  for(i=0;i<ret_n_u_hap;i++){
    if(ret_u_hap_list[i] != NULL){
      if(ret_u_hap_list[i]->loci != NULL) Free( ret_u_hap_list[i]->loci );
       Free( ret_u_hap_list[i] );
    }
  }

  Free(ret_u_hap_list);

  ret_u_hap_list = NULL;

  return;
}

/***********************************************************************************/

/*
     Algorithm AS 183 Appl. Statist. (1982) vol.31, no.2

     Returns a pseudo-random number rectangularly distributed
     between 0 and 1.   The cycle length is 6.95E+12 (See page 123
     of Applied Statistics (1984) vol.33), not as claimed in the
     original article.

     ix, iy and iz should be set to integer values between 1 and
     30000 before the first entry. To do this, 
     first call ranAS183_seed(iseed1,iseed2,iseed3), where iseed#
     are 3 long int seeds between 1 and 30000. The 3  seeds are
     saved, but ix,iy,iz can change.

    Translated from fortran to C.
*/

static long ix, iy, iz;

static int ranAS183_seed(int iseed1, int iseed2, int iseed3)
{
  int error;

  error=1;
  if( ( (iseed1 >=1) && (iseed1 <=30000)) && ( (iseed2 >=1) && (iseed2 <=30000) ) && 
      ( (iseed3 >=1) && (iseed3 <=30000) )) error=0;
  if(error) return (error);
  ix = iseed1;
  iy = iseed2;
  iz = iseed3;
  return (error);
}

/***********************************************************************************/

static double ranAS183()
{
   double u;

   ix = (171*ix) % 30269;
   iy = (172*iy) % 30307;
   iz = (170*iz) % 30323;
   u  = (double)ix/30269.0 + (double)iy/30307.0 + (double)iz/30323.0;
   return ( u - (int) u );
}

/***********************************************************************************/
static void errmsg(char *string){

  /* Function to emulate "stop" of S+ - see page 134, S Programing, by
     Venables and Ripley */

   PROBLEM "%s", string RECOVER(NULL_ENTRY);
}


/***********************************************************************************/

static void divideKeep(HAP **hap_list, long n, long *nReturn)
{
  long i,j;
  HAP *temp;
  long nValid = 0;


 i = -1;
  for(j = i+1; j<n; j++){
    if(hap_list[j]->keep !=0){
        i++;
        temp = hap_list[i];
        hap_list[i] = hap_list[j];
	hap_list[j] = temp;
  }
  }

 
  for(i = 0; i<n; i++){
    if( hap_list[i]->keep == 0) continue;
     nValid++;
  }


  *nReturn = nValid;


  return;
}

/***********************************************************************************/

static void add_more_memory(HAP ***hap_list, double **prior,long *max_haps){

  *max_haps = 2 * (*max_haps);

  *prior =  (double *) Realloc(*prior, *max_haps, double);
  if(prior==NULL){
    errmsg("could not realloc mem for prior");
  }

  *hap_list = (HAP **) Realloc(*hap_list, *max_haps, HAP* );
  if(hap_list==NULL){
    errmsg("could not realloc mem for hap_list");
  }

  return;
}

/***********************************************************************************/

static void insert_new_hap_pair(HAP ***hap_list_ptr, double **prior_ptr, 
                                long *max_haps, long insert_loc,
                                HAP *h1_old, HAP *h2_old, 
                                long a1_new, long a2_new,
                                long *pair_id_ptr, long *j){  
  HAP *h1_new, *h2_new;

 
  loci_used[insert_loc] = 1;

  h1_new = copy_hap(h1_old);
  h2_new = copy_hap(h2_old);

  h1_new->loci[insert_loc] = a1_new;
  h2_new->loci[insert_loc] = a2_new;

  (*pair_id_ptr) ++;
  h1_new->pair_id = (*pair_id_ptr);
  h2_new->pair_id = (*pair_id_ptr);

  if(  ((*j)+2)  >= (*max_haps) ){
     add_more_memory(hap_list_ptr, prior_ptr, max_haps);
  }

  (*j)++;
  (*hap_list_ptr)[*j] = h1_new;

  (*j)++;
  (*hap_list_ptr)[*j] = h2_new;

}
