#include <stdio.h>
#include "physics.h"
#include "protos.h"
#include "vop.h"

/* Make this return a ptr to static data so we can finesse the */
/* problem of what kind of FILE *! */
char *
PrintCellContents(const cell *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "\tmass: %.4g, nd:%d, bmax:%.2g, rcrit:%.2g\n"
		"\t" Sinfix("%.4f", " "), 
		Mass(p), p->daughters, p->bmax, p->rcrit, Vinfix(Pos(p), COMMA));
    }
    return contents_string;
}

char *
PrintBodyContents(const body *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "\tid:%d, mass:%.4g\n"
		"\t" Sinfix("%.4f", " "),
		p->ident, Mass(p), 
		Vinfix(Pos(p), COMMA));
    }
    return contents_string;
}

/* For out of bits confirmation */
char *
PrintBodyContentsLong(const body *p)
{
    static char contents_string[256];

    if( p == NULL){
	sprintf(contents_string, "\t(null)");
    }else{
	sprintf(contents_string, "\t@%#lx %7.4g %u\n"
		"\t" Sinfix("%14.10f", " ")
		" " Sinfix("%x", " "),
		(unsigned long)p, Mass(p), 
		(p->ident & (1<<31)) ? (p->ident & ~(1<<31)) : p->ident, 
		Vinfix(Pos(p), COMMA),
		Vinfix( *(int *)&p->pos, COMMA));
    }
    return contents_string;
}

char *PrintBranch(const cofmdata *cmp){
    static char ret[512];
    sprintf(ret, "Br: mass: %.3g, ndaughters=%d, pos=(%.3f %.3f %.3f), bmax:%.2g\n",
	   cmp->mass, cmp->ndaughters, 
	    cmp->pos[0], cmp->pos[1], cmp->pos[2],
	    cmp->bmax);
    return ret;
}
