/*
 * hash --
 *
 * 	This file declares hash table functions and data structures.
 *
 *  Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: hash.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 *************************************************************************
 */

#ifndef HASH_H_
#define HASH_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This type of structure is used to put a string key and associated
 * address into a hash table.
 */

struct Hash_Info {
    char *key;
    long val;
};

/*
 * Hash table entry structure.
 */

struct Hash_Entry {
    char *key;				/* String identifier */
    long val;				/* Value associated with string */
    struct Hash_Entry *next;		/* Pointer to next entry in bucket
					 * chain */
};

/*
 * Hash table structure.
 */

struct Hash_Table {
    struct Hash_Entry *entries;		/* Entries array */
    unsigned n_entries;			/* Number of entries */
    struct Hash_Entry **buckets;	/* Bucket array.  Each element is a
					 * linked list of entries */
    unsigned n_buckets;
};

/*
 * Global hash table functions.
 */

void Hash_InitTable(struct Hash_Table *);
void Hash_ClearTable(struct Hash_Table *);
int Hash_FillTable(struct Hash_Table *, unsigned, struct Hash_Info *);
int Hash_GetVal(struct Hash_Table *, char *, long *);
    
#ifdef __cplusplus
}
#endif

#endif
