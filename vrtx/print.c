#include <stdio.h>
#include "physics_vrtx.h"
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
	sprintf(contents_string, 
		"Str: %e %e %e, Pos: %e %e %e, rc2: %e",
		Strength(p)[0], Strength(p)[1], Strength(p)[2], 
		Pos(p)[0], Pos(p)[1], Pos(p)[2], p->rcrit2);
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
	sprintf(contents_string, 
		"Str: %e %e %e, Pos: %e %e %e",
		Strength(p)[0], Strength(p)[1], Strength(p)[2], 
		Pos(p)[0], Pos(p)[1], Pos(p)[2]);
    }
    return contents_string;
}

char *PrintBranch(const cofm_data *cmp){
    static char ret[512];
    sprintf(ret, "Br: daughters=%d, pos=(%.3f %.3f %.3f), bmax:%.2g\n",
	    cmp->daughters, 
	    cmp->pos[0], cmp->pos[1], cmp->pos[2],
	    cmp->bmax);
    return ret;
}
