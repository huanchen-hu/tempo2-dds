#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <math.h>
#include "tempo2.h"
#include <stdlib.h>
#include <iostream>
using namespace std;
/* Timing model       */
/* Based on bnrydds.f */
/* Orginally developed by M. Kramer */

double DDSmodel(pulsar *psr,int p,int ipos,int param)
{
    double an;
    double pb,k;
    double rad2deg = 180.0/M_PI;
    double SUNMASS = 4.925490947e-6;
    double delta = 0.0;
    double epsNum = 1.0e-12;   // numerical precision (seconds)
    double m2,tt0,t0,tt,x,ecc,er,eth,dr,dth,xdot,edot,am2,ct;
    double pbdot,xpbdot,phase,u,du,gamma;
    double orbits;
    int norbits;
    double cu,onemecu,cae,sae,ae,ae1,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
    double sqr1me2,cume,dlogbr,ds,da,a0,b0,delta_old,torb,diff;
    double csigma,ce,cx,comega,cgamma,cm2,csi;
    double shapmax,cshapmax,sdds,shaphof,cshaphof,brace,braceho;
    double epsLen,epsVel,dRotDefl;
    double droe,psi,cpsi,spsi,d2bar;
    const char *CVS_verNum = "$Id$";

    if (displayCVSversion == 1) CVSdisplayVersion("DDSmodel.C","DDSmodel()",CVS_verNum);


    if (psr[p].param[param_shapmax].paramSet[0]==1) shapmax = psr[p].param[param_shapmax].val[0];
    else shapmax = 0.0;

    if (psr[p].param[param_shaphof].paramSet[0]==1) shaphof = psr[p].param[param_shaphof].val[0];
    else shaphof = 0;
    
    //if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
    //else si = 0.0;

    if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
    else am2 = 0.0;

    if (psr[p].param[param_gamma].paramSet[0]==1) gamma = psr[p].param[param_gamma].val[0];
    else gamma = 0.0;

    if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0];
    else omz = 0.0;

    if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
    else xdot  = 0.0;

    if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
    else pbdot = 0.0;

    if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
    else edot = 0.0;

    if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
    else xpbdot = 0.0;

    if (psr[p].param[param_dr].paramSet[0] == 1) dr = psr[p].param[param_dr].val[0];
    else dr = 0.0;

    if (psr[p].param[param_dtheta].paramSet[0] == 1) dth = psr[p].param[param_dtheta].val[0];
    else dth = 0.0;

    if (psr[p].param[param_a0].paramSet[0] == 1) a0 = psr[p].param[param_a0].val[0];
    else a0 = 0.0;

    if (psr[p].param[param_b0].paramSet[0] == 1) b0 = psr[p].param[param_b0].val[0];
    else b0 = 0.0;
    
    pb = psr[p].param[param_pb].val[0]*SECDAY;
    an = 2.0*M_PI/pb;
    k = psr[p].param[param_omdot].val[0]/(rad2deg*365.25*86400.0*an);

    m2 = am2*SUNMASS;
    t0 = psr[p].param[param_t0].val[0];
    ct = psr[p].obsn[ipos].bbat;    

    tt0 = (ct-t0)*SECDAY;
    
    x = psr[p].param[param_a1].val[0]+xdot*tt0;
    ecc = psr[p].param[param_ecc].val[0]+edot*tt0;
    er = ecc*(1.0+dr);
    eth = ecc*(1.0+dth);
    
    sqr1me2=sqrt(1-pow(ecc,2));

do {
    delta_old = delta;
    tt = tt0-delta;
    orbits = tt/pb - 0.5*(pbdot+xpbdot)*pow(tt/pb,2);
    //orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
    norbits = (int)orbits;
    if (orbits<0.0) norbits--;
    phase=2.0*M_PI*(orbits-norbits);
    /*  Compute eccentric anomaly u by iterating Kepler's equation. */
    u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
    do {
        du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
        u=u+du;
    } while (fabs(du)>1.0e-14);

    /*  DD equations 17b, 17c, 29, and 46 through 52 */
    su=sin(u);
    cu=cos(u);
    onemecu=1.0-ecc*cu;
    cume=cu-ecc;
    
    cae=cume/onemecu;
    sae=sqr1me2*su/onemecu;
    ae1=atan2(sae,cae);
    if(ae1<0.0) ae1=ae1+2.0*M_PI;
    ae = ae1 + norbits*2.0*M_PI;
    //ae=2.0*M_PI*orbits + ae - phase;
    omega=omz/rad2deg + k*ae;
    sw=sin(omega);
    cw=cos(omega);
    
    psi  = omega + ae1; // angle w.r.t. ascending node
    spsi = sin(psi);
    cpsi = cos(psi);
    
    /* Roemer delay */
    alpha=x*sw;
    beta=x*sqrt(1-pow(eth,2))*cw;
    bg=beta+gamma;
    dre=alpha*(cu-er) + bg*su;
    drep=-alpha*su + bg*cu;
    drepp=-alpha*cu - bg*su;
    anhat=an/onemecu;

    /* Shapiro delay -- New additions for DDS model */
    sdds = 1.0-exp(-1.0*shapmax);    /* sidds in tempo */
    brace=onemecu-sdds*(sw*cume+sqr1me2*cw*su);
    /* **************************************** */
    if (psr[p].param[param_shaphof].paramSet[0]==0){
      dlogbr=log(brace);
      ds=-2.0*m2*dlogbr;

   }else{ /* Higher order corrections related to light propagation (for Double Pulsar only!) */
      double ratiompmc = 1.0714;  /* mass ratio mp/mc (= xc/xp) for Double Pulsar */
      double xc = x*ratiompmc;
      double xR = x+xc;         /* aR*sini/c [s] */
      double aR = xR/sdds;        /* aR [s] */

      /* Account for lensing contribution to propagation time
         simplified version (cf. Zschocke & Klioner 2010, eq. (73)) */
      epsLen = 2.0*m2/aR;

      /* 1.5pN contribution to Shapiro, i.e. leading order velocity dependence
         (Kopeikin & Schäfer 1999, eq. (130)) */
      epsVel = an*x/sdds*ratiompmc*ecc*su
             - an*x*sdds*ratiompmc/sqr1me2
             * (sw*cume + sqr1me2*cw*su)
             * (ecc*cw + (cw*cume - sqr1me2*sw*su)/onemecu);

      /* Shapiro delay incl. higher order corrections */
      braceho =  brace + (epsLen + epsVel) * shaphof;
      dlogbr  =  log(braceho);
      double dShaho  = -2.0*m2*dlogbr;

      /* Bending/lens-rotational delay
         - Doroshenko & Kopeikin 1995 approximation,
           including retardation to leading order (shift in c's position)
         - assumes pulsar spin to be perpendicular to the orbital plane,
           i.e. nearly perpendicular to the line-of-sight */
      double dfdt     = an*sqr1me2/pow(onemecu,2);
      double dpsi     = dfdt*ratiompmc*droe;
      double cpsiRet  = cpsi - spsi*dpsi;
      double braceRet = brace + epsVel * shaphof;
      dRotDefl = 2.0*m2/(2.0*M_PI*psr[p].param[param_f].val[0])/xR * cpsiRet/braceRet;

      /* Sum of all the signal propagation contributions */
      ds = dShaho + dRotDefl * shaphof;
    }
   
    /*** Aberration (DD) */
    da=a0*(spsi + ecc*sw) + b0*(cpsi + ecc*cw);

    // Sum of delays
    delta = dre + ds + da;
    diff = fabs(delta - delta_old);

  } while (diff > epsNum); // Inversion of timing model by iteration: end of loop              
    torb=-delta;
    
    /*  Now compute d2bar, the orbital time correction in DD equation 42. */
    //d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp - 0.5*ecc*su*dre*drep/onemecu)) + ds + da;
    //torb=-d2bar;

    if (param==-1) return torb;

    /*  Now we need the partial derivatives. Use DD equations 62a - 62k. */
    csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
    cx=sw*cume+sqr1me2*cw*su;
    comega=x*(cw*cume-sqr1me2*sw*su);
    cgamma=su;
    cm2=-2*dlogbr;
    double cdth;
    cdth=-ecc*ecc*x*cw*su/sqr1me2;
    
    /* New additions for the DDS model */
    if (psr[p].param[param_shaphof].paramSet[0]==0)
    {
        csi=2.0*m2*(sw*cume+sqr1me2*cw*su)/brace;
        cshaphof = 0.0;
    }
    else if (psr[p].param[param_shaphof].paramSet[0]==1)
    {
        csi=2.0*m2*(sw*cume+sqr1me2*cw*su)/braceho;
        cshaphof = -2.0*m2/braceho*(epsLen + epsVel) + dRotDefl;
    }
    cshapmax = csi * exp(-1.0*shapmax);
    
    cerr << "huhu: shaphof=" << shaphof << "huhu: cshaphof=" << cshaphof << endl;

    if (param==param_pb)
        return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
    else if (param==param_a1)
        return cx;
    else if (param==param_ecc)
        return ce;
    else if (param==param_om)
        return comega;
    else if (param==param_omdot)
        return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
    else if (param==param_t0)
        return -csigma*an*SECDAY;
    else if (param==param_pbdot)
        return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
    else if (param==param_sini)
        return csi;
    else if (param==param_gamma)
        return cgamma;
    else if (param==param_m2)
        return cm2*SUNMASS;
    else if (param==param_shapmax) /* New for DDS */
        return cshapmax;
    else if (param == param_shaphof)
        return cshaphof;
    else if (param==param_dtheta)
        return cdth;
    else if (param==param_a1dot) /* Also known as xdot */
        return cx*tt0;
    else if (param==param_edot)
        return ce*tt0;

    return 0.0;
}


void updateDDS(pulsar *psr,double val,double err,int pos)
{
    if (pos==param_pb)
    {
        psr->param[param_pb].val[0] += val/SECDAY;
        psr->param[param_pb].err[0]  = err/SECDAY;
    }
    else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_edot || pos==param_m2
            || pos==param_gamma || pos==param_shapmax || pos==param_shaphof || pos==param_dtheta)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_om)
    {
        psr->param[pos].val[0] += val*180.0/M_PI;
        psr->param[pos].err[0]  = err*180.0/M_PI;
    }
    else if (pos==param_pbdot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_a1dot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_omdot)
    {
        psr->param[pos].val[0] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
        psr->param[pos].err[0]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
    }
}
