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
	sprintf(contents_string, "\tstr: %.4g, nd:%d, nu:%d, nv:%d\n"
		"   " Sinfix("%.3f", " "), 
		p->strength, p->daughters, p->nu, p->nv, 
		Vinfix(p->pos, COMMA));
	if (p->nu) 
	    sprintf(contents_string, "%s   ffsf%d: %.8f %.8f %.8f %.8f",
		    contents_string,
		    p->nu, 
		    p->ffsf[0].r, p->ffsf[0].i, p->ffsf[1].r, p->ffsf[1].i);
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
	sprintf(contents_string, "\tid:%d, str:%.4g\n"
		" b " Sinfix("%.4f", " "),
		p->ident, p->strength, 
		Vinfix(p->pos, COMMA));
    }
    return contents_string;
}

char *PrintBranch(const cofmdata *cmp){
    static char ret[512];
    sprintf(ret, "Br: str: %.3g, ndaughters=%d, pos=(%.3f %.3f %.3f)\n",
	   cmp->strength, cmp->ndaughters, 
	    cmp->pos[0], cmp->pos[1], cmp->pos[2]);
    return ret;
}
