/* $Author: schaid $ */
/* $Date: 2005/03/01 23:06:04 $ */
/* $Header: /people/biostat3/sinnwell/Rdir/Make/RCS/haplo_em_pin.h,v 1.6 2005/03/01 23:06:04 schaid Exp $ */
/* $Locker:  $ */

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
*email: schaid@mayo.edu 
*/



/* redefine long to be int to be compatible with R */

#define long int

typedef struct HAP_T {
  long id;
  long code;
  long pair_id;
  long keep;
  long *loci;
  double post, wt;
} HAP;


typedef struct HAPUNIQUE_T {
  long code;
  long keep;
  long *loci;
  double prior;
} HAPUNIQUE;

static int iminarg1, iminarg2;

# define imin(a,b) (iminarg1=(a), iminarg2=(b), \
		    (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2) )

/* Windows compatibility */

#ifdef _WINDOWS
#define CDECL __cdecl
#else
#define CDECL
#endif


/************************** Function prototypes ***********************************/

static HAP* new_hap(long id, long pair_id, double wt, double prior, double post);

static void write_hap_list(HAP** so, long n_hap);


static int CDECL cmp_hap(const void *to_one, const void *to_two);

static int CDECL cmp_subId_hapPairId(const void *to_one, const void *to_two);

static int CDECL cmp_hap_code(const void *to_one, const void *to_two);

static long code_haps(long n_hap, HAP **hap_list);

static long hap_enum(HAP ***hap_list_ptr, double **prior_ptr, long *max_haps, long *n_alleles, long insert_loc, 
		     long n_hap, long *pair_id);

static HAP* copy_hap(HAP *old);

static long num_het(HAP* h1,HAP* h2);

static void hap_prior(long n_hap, HAP** hap_list, double *prior, long n_u_hap,
                      double min_prior);

static long hap_posterior(long n_hap, HAP **hap_list, double *prior, 
			  long n_u_hap, double min_posterior, double *lnlike);

static long **long_vec_to_mat(long *Yvec, long nrow, long ncol);

static long **long_matrix(long nrow, long ncol);

static void set_posterior(long n_hap, HAP **hap_list, long *random_start);

static int ranAS183_seed(int iseed1, int iseed2, int iseed3);

static double ranAS183();

static void errmsg(char *string);

static void compact(HAP **hap_list, long n, long *nReturn);

static HAPUNIQUE* copy_hap_unique(HAP *old, double *prior);

static void unique_haps(long n_hap, HAP **hap_list, HAPUNIQUE **u_hap_list, double *prior);

static long count_unique_haps(long n_hap, HAP **hap_list);

static void write_prior(long n, double *prior);

static void write_unique_hap_list(HAPUNIQUE** so, long n_hap);

static void divideKeep(HAP **hap_list, long n, long *nReturn);

static void add_more_memory(HAP ***hap_list, double **prior,long *max_haps);

static void insert_new_hap_pair(HAP ***hap_list_ptr, double **prior_ptr, 
                                long *max_haps, long insert_loc,
                                HAP *h1_old, HAP *h2_old, 
                                long a1_new, long a2_new,
                                long *pair_id, long *j);

static void overwrite_hap(HAP *new, HAP *old);
