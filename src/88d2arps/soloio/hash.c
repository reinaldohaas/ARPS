/*
 * hash --
 *
 *         This file defines hash table functions.
 *
 *         Reference:
 *         Kernighan, Brian W. and Rob Pike
 *         The practice of programming
 *         Reading, Massachusetts
 *         1999
 *
 *  Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: hash.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 *************************************************************************
 */

#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "errMsg.h"
#include "hash.h"

#define HASH_X 31

static unsigned hash(char *, unsigned);

/*
 *----------------------------------------------------------------------
 *
 * hash -
 *
 *         Retrieve an index in a hash table given the key.
 *
 * Arguments:
 *         char *k                - string key
 *         unsigned n        - number of buckets in hash table
 *
 *----------------------------------------------------------------------
 */

static unsigned hash(char *k, unsigned n)
{
    unsigned h;

    for (h = 0 ; *k != '\0'; k++) {
        h = HASH_X * h + (unsigned)*k;
    }
    return h % n;
}

/*
 *----------------------------------------------------------------------
 *
 * Hash_InitTable --
 *
 *         This constructor initializes a new hash table.
 *
 * Arguments:
 *         struct Hash_Table *tblP        - pointer to space for a hash table.
 *                                          Contents are assumed to be garbage.
 *
 * Side effects:
 *         Members in the struct are set to bogus values.  Upon return, it is
 *        safe to give the table to Hash_ClearTable.
 *
 *----------------------------------------------------------------------
 */

void Hash_InitTable(struct Hash_Table *tblP)
{
    tblP->entries = NULL;
    tblP->n_entries = 0;
    tblP->buckets = NULL;
    tblP->n_buckets = 0;
}

/*
 *----------------------------------------------------------------------
 *
 * Hash_ClearTable --
 *
 *         This destructor empties a hash table.
 *
 * Arguments:
 *         struct Hash_Table *tblP        - pointer to a hash table.
 *
 * Side effects:
 *         Memory associated with the table is freed and the table is
 *         reinitialized.
 *
 *----------------------------------------------------------------------
 */

void Hash_ClearTable(struct Hash_Table *tblP)
{
    struct Hash_Entry *ep, *entries;
    unsigned n_entries;

    if ( !tblP )  return;

    entries = tblP->entries;
    n_entries = tblP->n_entries;
    if (entries) {
      for (ep = entries; ep < entries + n_entries; ep++) {
          FREE(ep->key);
      }
    }
    FREE(tblP->entries);
    FREE(tblP->buckets);
    Hash_InitTable(tblP);
}

/*
 *----------------------------------------------------------------------
 *
 * Hash_FillTable --
 *
 *         Puts a set of entries into a hash table.
 *
 * Arguments:
 *         struct Hash_Table *tblP        - pointer to initialized hash table.
 *         unsigned n_entries                - number of entries
 *         struct Hash_Info *entries        - entry keys and data
 *
 * Results:
 *         The hash table is cleared and then repopulated with the given entries.
 *         Return value is true if table was successfully updated, otherwise
 *         return value is false.  If something goes wrong, the global error
 *         message is modified.
 *
 *----------------------------------------------------------------------
 */

int Hash_FillTable(struct Hash_Table *tblP, unsigned n_entries,
        struct Hash_Info *info)
{
    struct Hash_Entry *entries;                /* tblP->entries */
    struct Hash_Entry **buckets;        /* tblP->buckets */
    unsigned n_buckets;                        /* tblP->n_buckets */
    struct Hash_Info *ip;                /* Pointer into hash info array */
    struct Hash_Entry *ep;                /* Pointer into entry array */
    struct Hash_Entry **bp;                /* Pointer into bucket array */
    unsigned b;                                /* Index into buckets array */

    if ( !info || n_entries == 0) {
        return 1;
    }
    if (tblP->n_entries > 0) {
        Hash_ClearTable(tblP);
    }

    /*
     * Allocate and initialize memory.
     */

    n_buckets = n_entries;
    if (n_buckets % HASH_X == 0) {
        n_buckets++;
    }
    buckets = (struct Hash_Entry **)MALLOC(n_buckets * sizeof(struct Hash_Entry *));
    for (bp = buckets; bp < buckets + n_buckets; bp++) {
        *bp = NULL;
    }
    entries = (struct Hash_Entry *)MALLOC(n_entries * sizeof(struct Hash_Entry));
    for (ep = entries; ep < entries + n_entries; ep++) {
        ep->key = NULL;
        ep->next = NULL;
    }

    /*
     * Copy keys and values from info array to entries array.
     * Put entries into their buckets.
     */

    for (ip = info, ep = entries; ip < info + n_entries; ip++, ep++) {
        size_t len;
        char *s, *s1, *d;

        len = strlen(ip->key);
        ep->key = (char *)MALLOC(len + 1);

        for (s = ip->key, s1 = s + len, d = ep->key; s < s1; s++, d++) {
            *d = *s;
        }
        *d = '\0';
        ep->val = ip->val;
        b = hash(ep->key, n_buckets);
        ep->next = buckets[b];
        buckets[b] = ep;
    }

    tblP->entries = entries;
    tblP->n_entries = n_entries;
    tblP->buckets = buckets;
    tblP->n_buckets = n_buckets;
    return 1;
}

/*
 *----------------------------------------------------------------------
 *
 * Hash_GetVal --
 *
 *         This function retrieves a value from a hash table.
 *
 *
 * Arguments:
 *        struct Hash_Table *tblP        - pointer to space for a hash table.
 *        char *key                - string key.
 *        long *val                - recipient for the value from the table
 *                                  corresponding to key.
 *
 * Results:
 *        Return value is true if there is a value for key, otherwise return
 *        value is false.
 *
 *----------------------------------------------------------------------
 */

int Hash_GetVal(struct Hash_Table *tblP, char *key, long *val)
{
    unsigned b;                        /* Index into buckets array */
    struct Hash_Entry *ep;        /* Hash entry */

    b = hash(key, tblP->n_buckets);
    for (ep = tblP->buckets[b]; ep; ep = ep->next) {
        if (strcmp(ep->key, key) == 0) {
            *val = ep->val;
            return 1;
        }
    }
    return 0;
}
