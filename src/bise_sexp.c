# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "bise_sexp.h"
# include <R.h>
# include <Rinternals.h>

SEXP bise(SEXP Rd0, SEXP Rn0, SEXP RSP, SEXP RCUT, SEXP Rll, SEXP Rdb, SEXP Rnb, SEXP Rni, SEXP Rnirm)

 {
 
  int i,k,u,lc,lb,lt,lv1,l20,lg,vraimax,dsa,dsa2,kk;
  int d1[1100],dma[500],dmb[500],dg[366],db1[366],db2[366];
  double mb,sma,sma2,absma;
  double n1[1100],nma[500],nmb[500],b[500],s[250],ng[366];
  double nb1[366],nb2[366],s1[366],s2[366],dif[366],dif2[366];

//  SEXP Rlist = R_NilValue;
  SEXP Rval = R_NilValue;

//  PROTECT(Rlist = NEW_LIST(2));
// reinitialisation of brought in variables
  double *n0, *nb, *ni, *nirm;
  int *d0 = INTEGER(Rd0), *db = INTEGER(Rdb);
  int SP = asInteger(RSP), ll = asInteger(Rll), CUT = asInteger(RCUT);

 // fill copied variables
  n0 = REAL(Rn0);
  nb = REAL(Rnb); ni = REAL(Rni); nirm = REAL(Rnirm);

  PROTECT(Rval = allocVector(REALSXP, 366));
  
  lv1 = 2*SP+1;
  lc = ll;
  lt = ll*3;
 
 // triplement de la serie
 for(i=0;i<ll;i++)
   { d1[i]=d0[i];
     n1[i]=n0[i];}
 for(i=ll;i<(2*ll);i++)
   { d1[i]=d0[i-ll]+ll;
     n1[i]=n0[i-ll];}
 for(i=(2*ll);i<(3*ll);i++)
   { d1[i]=d0[i-2*ll]+2*ll;
     n1[i]=n0[i-2*ll];}

 for(i=0;i<(3*ll);i++) { if(n1[i]<0.09) n1[i]=0.09; }
  
// for lists to be constructed SET_ELEMENT(Rlist, 1, Rval1);
// recherche du 1er max
// cherche le (les) max de n1 a partir de offset SP et sur lc+offest

 mmax(n1,SP,lc,&dma[0],&nma[0]);

 vraimax=0;
// est-ce un vrai max????
//   si sur une periode 2*SP autour de ce point,
//  la moyenne des 20% valeurs les + elevees est inferieure a
//   ce point-0.25, on decide que ce max n'est pas le vrai max:

 while(vraimax==0)
   {
   k=0;
   for(i=dma[0]+ll-SP-1;i<dma[0]+ll+SP;i++)
     {
     dmb[k]=d1[i];
     nmb[k]=n1[i];
     k=k+1;
     }
   nmb[SP]=0;

   sortir(nmb,lv1,b);

   l20=(short)(lv1/5);

   mb=moy(b,l20);

     if((nma[0]-mb)>0.25)  // ce n'est pas un vrai max !!!
     {
     n1[dma[0]-1]=0.0;
     n1[dma[0]-1+ll]=0.0;
     n1[dma[0]-1+2*ll]=0.0;
   // recherche du 2eme max, puis 3eme, 4eme, etc
     mmax(n1,SP,lc,&dma[0],&nma[0]);
     }
   else
     vraimax=1;
   }

 // si pas de vraimax, ce n'est pas la peine de continuer
 if(vraimax==0)
   {
   for(k=1;k<ll;k++) {ni[k]=0.09;nirm[k]=0.09;}
   // *lbi=0;
   db[1]=1;   // idiot puisque pas de pt, mais j'ai peur qu'on ne puisse pas passer de array vides ...
   nb[1]=0.0;
   }

  else   // calcul de la pente sur la SP a partir du dernier point garde
   {
   u=0;
   while(dma[u]<(dma[0]+ll+ll-SP))
     {
     vraimax=0;
     while(vraimax==0)
       {
       if(nma[u]>0.09)
         {
         for(i=0;i<SP;i++)
           s[i]=(n1[dma[u]-1+i+1]-n1[dma[u]-1])/(i+1);
         mmax(s,0,SP,&dsa,&sma);
         for(i=dsa-1;i<SP;i++)
            s[i]=-9999.0;
         for(i=0;i<dsa-1;i++)
            if(n1[dma[u]-1+i+1]<=0.09)
              s[i]=-9999.0;

         mmax(s,0,SP,&dsa2,&sma2);
         absma=sma;
         if(sma<0) absma=-sma;
         if(((sma-sma2)/absma<0.4)||((n1[dma[u]+dsa-1]==0.09)&&(n1[dma[u]+dsa2-1]>0.09)))
           dsa=dsa2;
         dma[u+1]=dma[u]+dsa;
         nma[u+1]=n1[dma[u+1]-1];
         }
       else
         {
         i=0;
         while((n1[dma[u]+1+i]<=0.09)&&(i<SP))
            i=i+1;
         dma[u+1]=dma[u]+i+2;
         nma[u+1]=n1[dma[u+1]-1];
         }
      // co avant, on va regarder si c'est un vrai max ...
       k=0;
       for(i=dma[u+1]+ll-SP-1;i<dma[u+1]+ll+SP;i++)
         {
         dmb[k]=d1[i];
         nmb[k]=n1[i];
         k=k+1;
         }
       nmb[SP]=0;
       sortir(nmb,lv1,b);
       mb=moy(b,l20);

    // Elimination des points hauts
       if((nma[u+1]>0.2)&&(nma[u+1] > (double)CUT*mb)) // ce n'est pas un vrai max !!!
         {
         n1[dma[u+1]-1]=0.0;
         n1[dma[u+1]-1+ll]=0.0;
         n1[dma[u+1]-1+2*ll]=0.0;
         }
       else
         vraimax=1;
       }    // fin de   while(vraimax==0)
     u=u+1; 
     } // fin de while(dma[u]<(dma[0]+ll+(short)(ll/2)-SP)), c.a.d loop sur les pts hauts
 
   // selection des points retenus sur l'annee
   k=0;
   for(i=0;i<u;i++)
     if((dma[i]>ll)&&(dma[i]<(2*ll+1)))
       {
       db1[k]=dma[i]-ll;
       nb1[k]=nma[i];
       k=k+1;
       }
   lb=k;
   // si pas assez de pts retenus (ex.seulement 2), on met tout a zero!
   if(lb<3)
     {
     for(i=1;i<ll;i++)  {ni[i]=0.09;nirm[i]=0.09;}
     db[1]=1;
     nb[1]=0.0;
     // *lbi=0;
     }
   else
     {
     for(i=1;i<k;i++)
       {
       db[i]=db1[i];
       nb[i]=nb1[i];
       }
     // *lbi=lb;

     // annee etendue pour les besoins de l'interp.
     dg[0]=db[lb-1]-ll;
     ng[0]=nb[lb-1];
     for(i=0;i<lb;i++)
       { dg[i+1]=db[i];
         ng[i+1]=nb[i]; }
     dg[lb+1]=ll+db[1];
     ng[lb+1]=nb[1];
     lg=lb+2;

     // interpolation
     finterpol(dg,lg,ng,ll,ni);
     runmean(ni,ll,nirm);
     }
   } // end of pent loop
  
  for(i=0;i<ll;i++) {
    REAL(Rval)[i] = ni[i];
  }
  UNPROTECT(1);
  return(Rval);
}
// *************************************************************
 void mmax(double a[1500], int offset, int l, int *dada, double *mama)
  {
  int i;

  *mama=-99999.0;
  for(i=0+offset;i<(l+offset);i++)  if(a[i]>*mama) { *mama=a[i]; *dada=i+1;}
  }

// *************************************************************
 double dmax(double a[1500], int ll)
  {
  int i;
  double ma;

  ma=-99999.0;
  for(i=0;i<ll;i++)  if(a[i]>ma) ma=a[i];
  return(ma);
  }
// ************************************************************* 
 void sortir(double a[500], int l, double b[])
  {
  int i,d,lv2;
  double v;

  for(i=0;i<l;i++)
    {
    lv2=l-i;
    mmax(a,0,lv2,&d,&v);
    b[i]=v;
    a[d-1]=0;
    }
  }
// *************************************************************
 double moy(double a[500], int l20)
  {
  int i;
  double m;

  m=0;
  for(i=0;i<l20;i++)
    m=m+a[i];
  m=m/l20;
  return(m);
  }
// *************************************************************
 void finterpol(int *pa, int la, double *na, int ll, double *ni)
 {
  int i,ii;
  double nb[1000],da,db;

  for(i=1;i<la;i++)
    {
    da=pa[i]-pa[i-1];
    db=na[i]-na[i-1];
    for(ii=pa[i-1];ii<pa[i];ii++)
      {
      nb[ii-pa[0]]=na[i-1]+(ii-pa[i-1])*db/da;
      if((ii>0)&&(ii<(ll+1)))
         ni[ii-1]=nb[ii-pa[0]];
      }
   }
 }
// *************************************************************
 void runmean(double *ni, int ll, double *nirm)
  {
  int i,prn,demip,u;
  double nc[800],val;

  prn=31;demip=15;

  for(i=0;i<demip;i++)
     nc[i]=ni[i+ll-demip];
  for(i=demip;i<ll+demip;i++)
     nc[i]=ni[i-demip];
  for(i=ll+demip;i<ll+2*demip;i++)
     nc[i]=ni[i-ll-demip];

  for(i=1;i<ll;i++)
    {
    val=0.0;
    for(u=0;u<prn;u++)
       val=val+nc[i+u];
    nirm[i]=val/prn;
    }

  }
