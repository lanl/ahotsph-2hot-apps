#define KEYdotC
/* Most of the definitions are in key.h */
#include <stdio.h>		/* just for sprintf */
#include "key.h"
#include "protos.h"

/* Non-inlined definitions go here. */

char *
PrintKey(Key_t key)
{
    static char str[128];

#ifndef LONG_LONG_KEYS
    if (NK == 1)
      sprintf(str, "%lo", key.k[0]);
    else
#ifndef __alpha__
      /* This only works for NDIM==3 */
      sprintf(str, "%010lo.%01lo%010lo", key.k[1] >> 2, 
	      ((key.k[1] & 01) << 2) | (key.k[0] >> 30),
	      key.k[0] & ~(3U << 30));
#else
      /* This only works for NDIM==3 */
      sprintf(str, "%021lo.%01lo%021lo", key.k[1] >> 2, 
	      ((key.k[1] & 3) << 1) | (key.k[0] >> 63),
	      key.k[0] & ~(1L << 63));
#endif
#else
    sprintf(str, "%010o.%01o%010o", 
	    (unsigned long int)(key.k[0]>>33), /* bits 63-33 */
	    (unsigned long int)((key.k[0]>>30) & 7), /* bits 32,31,30 */
	    (unsigned long int)(key.k[0] & ((1<<30)-1))); /* 29-0 */
#endif
		
    return str;
}

int
TreeLevel(Key_t key, int ndim)
{
    int level;
    int chubits = (KEYBITS-1)/ndim;
    Key_t testkey;

    /* First check whether it's a 'body' (at the deepest level.)  */
    /* This will save considerable time... */
    testkey = KeyLshift(KeyInt(1), chubits*ndim);
    if( KeyEQ( testkey, KeyAnd(testkey, key) ) )
	return chubits;

    /* Now start looking from low levels */
    testkey = KeyInt(1);
    for (level = 0; level<chubits; level++){
	if( KeyEQ(key, testkey) )
	    return level;
	key = KeyRshift(key, ndim);
    }
    
    return -1;
}

/* Find the common level between two 'body'-keys */
int
CommonLev(Key_t bkey1, Key_t bkey2, int ndim)
{
    int level = (KEYBITS-1)/ndim;
    Key_t key0 = KeyInt(0);

    for (bkey1 = KeyXOR(bkey1, bkey2); KeyNEQ(bkey1, key0); level--)
      bkey1 = KeyRshift(bkey1, ndim);
    
    return(level);
}

int
KeyContained(Key_t outer, Key_t key, int ndim)
{
    int ret;
    int difference;
    
    /* What if difference is negative?! */
    difference = TreeLevel(key, ndim) - TreeLevel(outer, ndim);
    key = KeyRshift(key, ndim*difference);
    ret = KeyEQ(key, outer);
    return(ret);
}
