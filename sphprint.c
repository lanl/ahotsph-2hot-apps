#include <stdio.h>
#include "physics_sph.h"
#include "protos.h"
#include "vop.h"

/* Make this return a ptr to static data so we can finesse the */
/* problem of what kind of FILE *! */
char *
PrintSPHCellContents(const SPHcell *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "\tmass: %.4g, nd:%d, bmax:%.2g, rcrit:%.2g\n"
		"\t" Sinfix("%.4f", " "), 
		p->mass, p->daughters, p->bmax, p->rcrit, Vinfix(p->pos, COMMA));
    }
    return contents_string;
}

char *
PrintSPHBodyContents(const SPHbody *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "id:%d mass:%.4g nbrs:%d\n"
		" pos:" Sinfix("%.4f", " ") " vel:" Sinfix("%.4f", " ")
		"\nh:%g hdot:%g rho:%g rho_est:%g drho_dt:%g\n"
		"u:%g udot:%g pr:%g vsound:%g",
		p->ident, p->mass, p->nbrs,
		Vinfix(p->pos, COMMA), Vinfix(p->vel, COMMA), 
		p->h, p->hdot, p->rho, p->rho_est, p->drho_dt,
		p->u, p->udot, p->pr, p->vsound);
    }
    return contents_string;
}

/* For out of bits confirmation */
char *
PrintSPHBodyContentsLong(const SPHbody *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "\t@%#lx %7.4g\n"
		"\t" Sinfix("%14.10f", " ")
		" " Sinfix("%x", " "),
		(unsigned long)p, p->mass, Vinfix(p->pos, COMMA),
		Vinfix( *(int *)&p->pos, COMMA));
    }
    return contents_string;
}

char *PrintSPHBranch(const SPHcofmdata *cmp){
    static char ret[512];
    sprintf(ret, "Br: mass: %.3g, ndaughters=%d, pos=(%.3f %.3f %.3f), bmax:%.2g\n",
	   cmp->mass, cmp->ndaughters, 
	    cmp->pos[0], cmp->pos[1], cmp->pos[2],
	    cmp->bmax);
    return ret;
}
