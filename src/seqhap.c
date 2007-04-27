/* $Author: sinnwell $ */
/* $Date: 2007/04/06 19:30:01 $ */
/* $Header: /people/biostat3/sinnwell/Haplo/Make/RCS/seqhap.c,v 1.1 2007/04/06 19:30:01 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: seqhap.c,v $
 * Revision 1.1  2007/04/06 19:30:01  sinnwell
 * Initial revision
 * * 
 */
/*Author: Yu*/
/*Date: 2007/04/02*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

static long ranAS183_seed(long, long, long);
static double ranAS183();
static long creatsubhap();
static double chisq2h(long *);
static double mantel(long *);
static double r2test(long, long);
static double combine(long *, long, long);
static void permute(long *, double *, long *);
static double pchisq(double, double);

static long *subjid, *hap1code, *hap2code, *disease, *newhap1code, *newhap2code,
  *inlist, inlist_length;
static long **haplist;
static double *post, *pos;
static long *newhap1codesingle, *newhap2codesingle;

static double r2_threshold, mh_threshold, haplo_freq_min; 
static long nsnp, nhap, nsub, npost, newnhap, flag, hap_df;


void seqhap(
	    long *nsnp_c,
	    long *nsub_c,
	    long *npost_c,
	    long *nhap_c,
	    long *subjid_c,
	    long *hap_c,
	    long *hap1_c,
	    long *hap2_c,
	    long *disease_c,

	    double *post_c,
	    double *pos_c,

	    long *seed_c,
	    long *nperm_c,
	    double *lamda_c,
	    double *r2_threshold_c,
	    double *haplo_freq_min_c,

	    long *inlist_c,
	    long *hap_df_c,
	    double *hap_chi_c,
	    double *hap_p_point_c,
	    double *hap_p_region_c,

	    long *sum_df_c,
	    double *sum_chi_c,
	    double *sum_p_point_c,
	    double *sum_p_region_c,

	    double *chi_chi_c,
	    double *chi_p_point_c,
	    double *chi_p_region_c
)
{

  /*read data passed from S*/
  nsnp = *nsnp_c;
  r2_threshold = *r2_threshold_c;
  mh_threshold = *lamda_c;
  long N_PERM = *nperm_c;
  nsub = *nsub_c;
  npost = *npost_c;
  nhap = *nhap_c;
  haplo_freq_min = *haplo_freq_min_c;

  FILE *fp;
  long nperm, si, i, j, k;
  double chi_stat0i[nsnp], chi_max=0, chi_max0=0, chi_pi[nsnp], chi_p=0,
    hap_stat0i[nsnp], hap_min=1, hap_min0=1, hap_pi[nsnp], hap_p=0,
    sum_stat0i[nsnp], sum_min=1, sum_min0=1, sum_pi[nsnp], sum_p=0; 
  for(i=0; i<nsnp; i++){
    chi_pi[i]=0; hap_pi[i]=0; sum_pi[i]=0;}

  long seed[3];
  for(i=0; i<3; i++) seed[i]=seed_c[i];
  double ran1();
  ranAS183_seed(seed[0],seed[1],seed[2]);

  haplist = (long **) malloc (nhap *sizeof (long *));
  for (i=0; i<nhap; i++)
      haplist[i] = (long *) malloc ((unsigned) nsnp * sizeof (long));

  /*read haplotype list passed from S*/
  k=0;
  for(i=0; i<nsnp; i++)
      for(j=0; j<nhap; j++)
	{
	  haplist[j][i]=hap_c[k];
	  k++;
	}

  subjid= (long *) malloc(npost *sizeof(long));
  hap1code = (long *) malloc(npost *sizeof(long));
  hap2code = (long *) malloc(npost *sizeof(long));
  disease = (long *)  malloc(npost *sizeof(long));
  newhap1code = (long *) malloc(npost *sizeof(long));
  newhap2code = (long *) malloc(npost *sizeof(long));
  inlist = (long *) malloc(nsnp *sizeof(long));
  post = (double *) malloc(npost *sizeof(double));
  pos = (double *) malloc(nsnp *sizeof(double));
  newhap1codesingle = (long *) malloc(npost *sizeof(long));
  newhap2codesingle = (long *) malloc(npost *sizeof(long));

  /*read data passed from S*/
  for(i=0; i<npost; i++)
    {
      subjid[i]=subjid_c[i];
      hap1code[i]=hap1_c[i]; hap2code[i]=hap2_c[i];
      disease[i]=disease_c[i];
      post[i]=post_c[i];
    }

  for(i=0; i<nsnp; i++)  pos[i]=pos_c[i];
  double doublek; 
  j=0;

  /*calculation using original disease status*/
  printf("\n");
  for(si=0; si<nsnp; si++)
    {
      inlist_length=1; inlist[inlist_length-1]=si;
      flag=1, k=1, hap_df=1; 
      newnhap=creatsubhap();
      chi_stat0i[si]=chisq2h(disease);
      sum_stat0i[si] = chi_stat0i[si];
      while(flag!=0 && (si+k<nsnp || si-k>=0) )
	{
	  sum_stat0i[si]=sum_stat0i[si] + combine(disease, si, k);
	  k++;
	}
      if(inlist_length==1) {hap_stat0i[si]=chi_stat0i[si]; hap_df=1;}
      else hap_stat0i[si] = chisq2h(disease);
      /*
      printf("snp%d\n",si+1);
      printf("single-locus method:  chi-square %lf, d.f. %d\n", chi_stat0i[si],1);
      printf("Sequential method combines snps ");
      for(i=0; i<inlist_length; i++) printf("%d ",inlist[i]+1); printf("\n");
      printf("Seq-haplotype method: chi-square %lf, d.f. %d\n", hap_stat0i[si],hap_df);
      printf("Seq-summary method:   chi-square %lf, d.f. %d\n\n", sum_stat0i[si],inlist_length);
      */

      /*save for Splus*/
      for(i=0; i<inlist_length; i++) {inlist_c[j]=inlist[i]+1;j++;}
      for(i=inlist_length; i<nsnp; i++) {inlist_c[j]=0; j++;}
      hap_df_c[si]=hap_df;
      hap_chi_c[si]=hap_stat0i[si];
      sum_df_c[si]=inlist_length;
      sum_chi_c[si]=sum_stat0i[si];
      chi_chi_c[si]=chi_stat0i[si];

      doublek = inlist_length;  
      sum_stat0i[si] = pchisq(doublek/2, sum_stat0i[si]/2);
      doublek = hap_df; 
      hap_stat0i[si] = pchisq(doublek/2, hap_stat0i[si]/2);
      if(chi_stat0i[si]>chi_max0) chi_max0=chi_stat0i[si];
      if(hap_stat0i[si]<hap_min0) hap_min0=hap_stat0i[si];
      if(sum_stat0i[si]<sum_min0) sum_min0=sum_stat0i[si];
    }
 
 /*calculate statistics using permutated data*/
  double perm_rand[nsub];
  long perm_disease[npost], nondup_disease[npost], subduplicated[npost];
  for(i=0; i<npost; i++)
    { subduplicated[i]=0;
      j = 0;
      while(j<i && subduplicated[i]==0) {
	if(subjid[i]==subjid[j]) subduplicated[i] = 1;
	j++;}
    }

  j = 0;
  for(i=0; i<npost; i++)
    if(subduplicated[i]==0){
	nondup_disease[j] = disease[i];
	j ++;}

  double chi_tmp, hap_tmp, sum_tmp;
  for(nperm=0; nperm<N_PERM; nperm++)
    {
      /*permute disease status*/
      for(i=0; i<nsub; i++)
	perm_rand[i]=ranAS183();
      permute(nondup_disease, perm_rand, perm_disease);
      chi_max=0; hap_min=1; sum_min=1;

      for(si=0; si<nsnp; si++)
	{
	  inlist_length=1; inlist[inlist_length-1]=si;
	  flag=1, k=1, hap_df=1; 
	  newnhap=creatsubhap();
	  chi_tmp=chisq2h(perm_disease);
	  sum_tmp = chi_tmp;
	  while(flag!=0 && (si+k<nsnp || si-k>=0) )
	    {
	      sum_tmp=sum_tmp + combine(perm_disease, si, k);
	      k++;
	    }
	  doublek = inlist_length; 
	  sum_tmp = pchisq(doublek/2, sum_tmp/2);
	  if(inlist_length==1) {hap_tmp=chi_tmp; hap_df=1;}
	  else hap_tmp = chisq2h(perm_disease);
	  doublek = hap_df; 
	  hap_tmp = pchisq(doublek/2, hap_tmp/2);

	  /*global*/
	  if(chi_tmp>chi_max) chi_max=chi_tmp;
	  if(hap_tmp<hap_min) hap_min=hap_tmp;
	  if(sum_tmp<sum_min) sum_min=sum_tmp;	  

	  /*pointwise*/
	  if(chi_tmp>chi_stat0i[si]) chi_pi[si]++;
	  if(hap_tmp<hap_stat0i[si]) hap_pi[si]++;
	  if(sum_tmp<sum_stat0i[si]) sum_pi[si]++;
	}
      /*global*/
      if(chi_max>chi_max0) chi_p++;
      if(hap_min<hap_min0) hap_p++;
      if(sum_min<sum_min0) sum_p++;
    }

  /*global*/
  *chi_p_region_c=chi_p/N_PERM;
  *hap_p_region_c=hap_p/N_PERM;
  *sum_p_region_c=sum_p/N_PERM;

  /*pointwisw*/
  for(si=0; si<nsnp; si++)/*3rd and beyond are pointwise p-values*/
    {
      chi_p_point_c[si]=chi_pi[si]/N_PERM;
      hap_p_point_c[si]=hap_pi[si]/N_PERM;
      sum_p_point_c[si]=sum_pi[si]/N_PERM;
    }


  free(subjid); free(hap1code); free(hap2code); free(disease); 
  free(newhap1code); free(newhap2code);free(inlist); 
  free(haplist); free(post); free(pos);
  free(newhap1codesingle); free(newhap2codesingle);

}

/*the combine function performs the sequential search*/
static double combine(long *d, long si, long k)
{
  long j, flagl, flagr;
  double mh_sum=0, mh_tmp;
	  if(si-k>=0 && si+k<nsnp)
	    if(pos[si]-pos[si-k] > pos[si+k]-pos[si])
	      {
		/*check si+k*/
		flagr=1;
		if(r2test(si, si+k)<r2_threshold) 
		  {
		    for(j=0; j<npost; j++)
		      {
			newhap1codesingle[j] = haplist[hap1code[j]-1][si+k];
			newhap2codesingle[j] = haplist[hap2code[j]-1][si+k];
		      }
		    mh_tmp=mantel(d);
		    if(mh_tmp>mh_threshold){
		      mh_sum=mh_sum + mh_tmp; 
		      inlist_length++;
		      inlist[inlist_length-1]=si+k;
		      newnhap=creatsubhap();}
		    else flagr=0;
		  }

		/*if r2>r2_threshold, do NOT add si+k into inlist, but continue search;
		  if r2<r2_threshold && mh_tmp>mh_threshold, add si+k into inlist;
		  if r2<r2_threshold && mh_tmp<mh_threshold, do NOT add si+k into inlist,
  		  but continue search*/

		/*check si-k */
		flagl=1;
		if(r2test(si,si-k)<r2_threshold) 
		  {
		    for(j=0; j<npost; j++)
		      {
			newhap1codesingle[j] = haplist[hap1code[j]-1][si-k];
			newhap2codesingle[j] = haplist[hap2code[j]-1][si-k];
		      }
		    mh_tmp=mantel(d);
		    if(mh_tmp>mh_threshold){
		      mh_sum=mh_sum + mh_tmp;
		      inlist_length++;
		      inlist[inlist_length-1]=si-k;
		      newnhap=creatsubhap();}
		    else flagl=0;
		  }

		if(flagl!=1 && flagr!=1) flag=0;
	      }
	    else
	      {
		/*check si-k*/
		flagl=1;
		if(r2test(si,si-k)<r2_threshold) 
		  {
		    for(j=0; j<npost; j++)
		      {
			newhap1codesingle[j] = haplist[hap1code[j]-1][si-k];
			newhap2codesingle[j] = haplist[hap2code[j]-1][si-k];
		      }
		    mh_tmp=mantel(d);
		    if(mh_tmp>mh_threshold){
		      mh_sum=mh_sum + mh_tmp;
		      inlist_length++;
		      inlist[inlist_length-1]=si-k;
		      newnhap=creatsubhap();}
		    else flagl=0;
		  }

		/*check si+k */
		flagr=1;
		if(r2test(si, si+k)<r2_threshold) 
		  {
		    for(j=0; j<npost; j++)
		      {
			newhap1codesingle[j] = haplist[hap1code[j]-1][si+k];
			newhap2codesingle[j] = haplist[hap2code[j]-1][si+k];
		      }
		    mh_tmp=mantel(d);
		    if(mh_tmp>mh_threshold){
		      mh_sum=mh_sum + mh_tmp;
		      inlist_length++;
		      inlist[inlist_length-1]=si+k;
		      newnhap=creatsubhap();}
		    else flagr=0;
		  }

		if(flagl!=1 && flagr!=1) flag=0;
	      }

	  if(si-k<0 && si+k<nsnp)
	    {
	      /*check si+k*/
	      flagr=1;
	      if(r2test(si, si+k)<r2_threshold) 
		{
		  for(j=0; j<npost; j++)
		    {
		      newhap1codesingle[j] = haplist[hap1code[j]-1][si+k];
		      newhap2codesingle[j] = haplist[hap2code[j]-1][si+k];
		    }
		  mh_tmp=mantel(d);
		  if(mh_tmp>mh_threshold){
		    mh_sum=mh_sum + mh_tmp;
		    inlist_length++;
		    inlist[inlist_length-1]=si+k;
		    newnhap=creatsubhap();}
		  else flagr=0;
		}
	      if(flagr!=1) flag=0;
	    }
	  if(si+k>=nsnp && si-k>0)
	    {
	      /*check si-k*/
	      flagl=1;
	      if(r2test(si, si-k)<r2_threshold) 
		{
		  for(j=0; j<npost; j++)
		    {
		      newhap1codesingle[j] = haplist[hap1code[j]-1][si-k];
		      newhap2codesingle[j] = haplist[hap2code[j]-1][si-k];
		    }
		  mh_tmp=mantel(d);
		  if(mh_tmp>mh_threshold){
		    mh_sum=mh_sum + mh_tmp;
		    inlist_length++;
		    inlist[inlist_length-1]=si-k;
		    newnhap=creatsubhap();}
		  else flagl=0;
		}
	      if(flagl!=1) flag=0;
	    }
	  return(mh_sum);
}

/*the permute funciton permutes disease status*/
static void permute(long *nondup_disease, double *perm_rand, long *perm_disease)
{
  long i, j, tmpint, perm_order[nsub];
  double tmp;
  for(i=0; i<nsub; i++)  perm_order[i] = i; /*initialize*/

  for(i=0; i<(nsub-1); i++) 
    for(j=0; j<(nsub-1-i); j++)
      if(perm_rand[j+1] < perm_rand[j])
	{
	  tmp = perm_rand[j];
	  perm_rand[j] = perm_rand[j+1];
	  perm_rand[j+1] = tmp;
	  tmpint = perm_order[j];
	  perm_order[j] = perm_order[j+1];
	  perm_order[j+1] = tmpint;
	}
  for(i=0; i<npost; i++) perm_disease[i] = nondup_disease[perm_order[subjid[i]-1]];  
}

/*the r2test calculates the square of the Pearson's correlation*/
static double r2test(long site1, long site2)
{
  long i;
  double n00=0, n01=0, n10=0, n11=0, r2;
  for(i=0; i<npost; i++)
    {
      switch(haplist[hap1code[i]-1][site1])
	{
	case 0: 
	  if(haplist[hap1code[i]-1][site2] ==0) n00=n00+post[i]; else n01=n01+post[i]; break;
	case 1:  
	  if(haplist[hap1code[i]-1][site2] ==0) n10=n10+post[i]; else n11=n11+post[i]; break;
	}
      switch(haplist[hap2code[i]-1][site1])
	{
	case 0: 
	  if(haplist[hap2code[i]-1][site2] ==0) n00=n00+post[i]; else n01=n01+post[i]; break;
	case 1:  
	  if(haplist[hap2code[i]-1][site2] ==0) n10=n10+post[i]; else n11=n11+post[i]; break;
	}
    }
	
  r2 = pow((n00*n11-n01*n10),2)/((n00+n01)*(n00+n10)*(n01+n11)*(n10+n11)); 
  return(r2);
}

/*the mantel function calculates the Mantel-Haenszel statistic*/
static double mantel(long *d)
{
  long i,j;/*a stratum is labelled if a stratum has counts <=2*/
  double mh_stat=0, maxcount=0, num=0, denom=0, n00[newnhap], n0plus[newnhap], 
    n1plus[newnhap], nplus0[newnhap],nplus1[newnhap], ntotal[newnhap], n=0;
  /*initialize*/
  for(i=0; i<newnhap; i++)
    {
      n00[i]=0; n0plus[i]=0; n1plus[i]=0; nplus0[i]=0;
      nplus1[i]=0; ntotal[i]=0;
    }
  /*count*/
  for(i=0; i<npost; i++)
    {
      switch(d[i])
	{
	case 0: {
	  if(newhap1codesingle[i]==0){
	    n00[newhap1code[i]-1] = n00[newhap1code[i]-1] + post[i];
	    n0plus[newhap1code[i]-1] = n0plus[newhap1code[i]-1] + post[i];
	    nplus0[newhap1code[i]-1] = nplus0[newhap1code[i]-1] + post[i];
	    ntotal[newhap1code[i]-1] = ntotal[newhap1code[i]-1] + post[i];
	  }
	  else {
	    n1plus[newhap1code[i]-1] = n1plus[newhap1code[i]-1] + post[i];
	    nplus0[newhap1code[i]-1] = nplus0[newhap1code[i]-1] + post[i];
	    ntotal[newhap1code[i]-1] = ntotal[newhap1code[i]-1] + post[i];
	  }
	  if(newhap2codesingle[i]==0){
	    n00[newhap2code[i]-1] = n00[newhap2code[i]-1] + post[i];
	    n0plus[newhap2code[i]-1] = n0plus[newhap2code[i]-1] + post[i];
	    nplus0[newhap2code[i]-1] = nplus0[newhap2code[i]-1] + post[i];
	    ntotal[newhap2code[i]-1] = ntotal[newhap2code[i]-1] + post[i];
	  }
	  else {
	    n1plus[newhap2code[i]-1] = n1plus[newhap2code[i]-1] + post[i];
	    nplus0[newhap2code[i]-1] = nplus0[newhap2code[i]-1] + post[i];
	    ntotal[newhap2code[i]-1] = ntotal[newhap2code[i]-1] + post[i];
	  }
	  break;}
	case 1: {
	  if(newhap1codesingle[i]==0){
	    n0plus[newhap1code[i]-1] = n0plus[newhap1code[i]-1] + post[i];
	    nplus1[newhap1code[i]-1] = nplus1[newhap1code[i]-1] + post[i];
	    ntotal[newhap1code[i]-1] = ntotal[newhap1code[i]-1] + post[i];
	  }
	  else {
	    n1plus[newhap1code[i]-1] = n1plus[newhap1code[i]-1] + post[i];
	    nplus1[newhap1code[i]-1] = nplus1[newhap1code[i]-1] + post[i];
	    ntotal[newhap1code[i]-1] = ntotal[newhap1code[i]-1] + post[i];
	  }
	  if(newhap2codesingle[i]==0){
	    n0plus[newhap2code[i]-1] = n0plus[newhap2code[i]-1] + post[i];
	    nplus1[newhap2code[i]-1] = nplus1[newhap2code[i]-1] + post[i];
	    ntotal[newhap2code[i]-1] = ntotal[newhap2code[i]-1] + post[i];
	  }
	  else {
	    n1plus[newhap2code[i]-1] = n1plus[newhap2code[i]-1] + post[i];
	    nplus1[newhap2code[i]-1] = nplus1[newhap2code[i]-1] + post[i];
	    ntotal[newhap2code[i]-1] = ntotal[newhap2code[i]-1] + post[i];
	  }
	  break;}
	}
    }
  for(i=0; i<newnhap; i++)
    {
      if(n0plus[i]>0.5 && n1plus[i]>0.5 && nplus0[i]>0.5 && nplus1[i]>0.5) 
	{
	  num = num + n00[i] - n0plus[i] * nplus0[i]/ntotal[i]; 
	  denom = denom + n0plus[i]*n1plus[i]*nplus0[i]*nplus1[i]/
	    (ntotal[i]*ntotal[i]*(ntotal[i]-1));
	}
      else 
	{
	  if(maxcount < ntotal[i])
	    maxcount = ntotal[i];
	}
      n = n + ntotal[i];    
    }

  if(denom==0) return(denom);
  mh_stat = pow(num,2)/denom;
  return(mh_stat);
  /*
    for(i=0; i<newnhap; i++)
      {
	printf("n00=%lf n0plus=%lf n1plus=%lf nplus0=%lf nplus1=%lf\n",
	       n00[i], n0plus[i], n1plus[i], nplus0[i], nplus1[i]);
      }
	printf("num=%lf, denom=%lf, mh_stat=%lf\n", num, denom, mh_stat);
	  printf("mh_stat\n%lf\n", mh_stat);
  */
}

/*the chisq2h function calculates the chisquare statistic for a 2-by-h table */
static double chisq2h(long *d)
{
  double n0[newnhap], n1[newnhap], n0sum=0, n1sum=0, nhapsum[newnhap];
  double en0[newnhap], en1[newnhap], ntotal=0, chi2h_stat=0;
  long i, j;

  /*initialize*/
  for(i=0; i<newnhap; i++)
    {n0[i]=0; n1[i]=0; nhapsum[i]=0; en0[i]=0; en1[i]=0;}
  for(i=0; i<npost; i++)
    switch(d[i])
      {
      case 0: {
	n0[newhap1code[i]-1] = n0[newhap1code[i]-1] + post[i]; 
	n0[newhap2code[i]-1] = n0[newhap2code[i]-1] + post[i]; 
	nhapsum[newhap1code[i]-1] = nhapsum[newhap1code[i]-1] + post[i];
	nhapsum[newhap2code[i]-1] = nhapsum[newhap2code[i]-1] + post[i];
	n0sum = n0sum + 2*post[i];
	break;}
      case 1: {
	n1[newhap1code[i]-1] = n1[newhap1code[i]-1] + post[i]; 
	n1[newhap2code[i]-1] = n1[newhap2code[i]-1] + post[i]; 
	nhapsum[newhap1code[i]-1] = nhapsum[newhap1code[i]-1] + post[i];
	nhapsum[newhap2code[i]-1] = nhapsum[newhap2code[i]-1] + post[i];
	n1sum = n1sum + 2*post[i];
	break;}
      }
  ntotal = n0sum + n1sum;

  long ntrunhap=0; 
  double trunhapsum[newnhap], trunn0[newnhap], trunn1[newnhap]; 
  for(i=0; i<newnhap; i++)
    {
      if(nhapsum[i]>haplo_freq_min*ntotal)
	{
	  trunhapsum[ntrunhap]=nhapsum[i];
	  trunn0[ntrunhap] = n0[i];
	  trunn1[ntrunhap] = n1[i];
	  ntrunhap++;
	}
      else
	{
	  n0sum = n0sum - n0[i];
	  n1sum = n1sum - n1[i];
	}
    }
  ntotal = n0sum + n1sum;

  for(i=0; i<ntrunhap; i++)
    {
      en0[i] = trunhapsum[i] * n0sum /ntotal;
      en1[i] = trunhapsum[i] * n1sum /ntotal;
      chi2h_stat = chi2h_stat + pow(trunn0[i]-en0[i],2)/en0[i] + 
	pow(trunn1[i]-en1[i],2)/en1[i];
    }
  /*
  printf("k=%d\n", k);
  printf("chisq2htest:\n");
  for(i=0; i<ntrunhap; i++)
    printf("n0=%lf n1=%lf en0=%lf en1=%lf %lf\n", trunn0[i], trunn1[i], en0[i], en1[i], trunhapsum[i]);
  printf("n0sum=%lf n1sum=%lf ntotal=%lf\n", n0sum, n1sum, n0sum+n1sum);
  printf("chi2h=%lf\n", chi2h_stat);
  */
  hap_df=ntrunhap-1;
  return(chi2h_stat);
}

/*the creatsubhap function identifies the subhaplotypes for a chosen set of loci*/
static long creatsubhap()
{
  long i, j, a, duplicated[nhap];
  newnhap =0;
  char temp1[inlist_length];
  char temp2[inlist_length];
  long newhapcode[nhap];

  /*initialize duplicated*/
  for(i=0; i<nhap; i++)  duplicated[i] = 0;
 
  for(i=0; i<nhap; i++)
      {
	duplicated[i] =0;
	j = 0;
	while(j<i && duplicated[i]==0)
	  {
	    for(a=0; a<inlist_length; a++){
	      sprintf(&temp1[a], "%d", haplist[i][inlist[a]]);
	      sprintf(&temp2[a], "%d", haplist[j][inlist[a]]);}
	    if(strncmp(temp1, temp2,inlist_length)==0) 
	      {
		duplicated[i] = 1;
		newhapcode[i] = newhapcode[j];
	      }
	    j++;
	  }
	if(duplicated[i]==0)  
	  {
	    newnhap++; 
	    newhapcode[i] = newnhap;
	  }
      }

  /*label the haplotype code with new value*/
  for(i=0; i<npost; i++)
    {
      newhap1code[i] = newhapcode[hap1code[i]-1];
      newhap2code[i] = newhapcode[hap2code[i]-1];
    }

  /*
  printf("\n\nprintout from creatsub\n");
  for(i=0; i<nhap; i++) printf("%d ", duplicated[i]); printf("newnhap=%d\n", newnhap);
  for(i=0; i<nhap; i++) printf("%d,%d ", i+1,newhapcode[i]); printf("\n");
  for(i=0; i<nhap; i++) printf("%d ", newhap1code[i]); printf("\n");
  for(i=0; i<nhap; i++) printf("%d ", newhap2code[i]); printf("\n");
  printf("createsub:\nnewnhap=%d\n",newnhap);*/

  return(newnhap);

}

/*the pchisq function calculates the cumulative distribution function of 
a chi-square distribution*/
double pchisq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;
/*
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
*/
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 3#Q)$. */

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	long j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 3#Q)$. */
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
/*
	void nrerror(char error_text[]);
*/
	long i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
/*
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
*/
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software 3#Q)$. */

#define ITMAX 100
#define EPS 3.0e-7

void gser (double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
/*
	void nrerror(char error_text[]);
*/
	long n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
/*
		if (x < 0.0) nrerror("x less than 0 in routine gser");
*/
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
/*
		nrerror("a too large, ITMAX too small in routine gser");
*/
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 3#Q)$. */

static long ix, iy, iz;

static long ranAS183_seed(long iseed1, long iseed2, long iseed3) {
  long error;

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
