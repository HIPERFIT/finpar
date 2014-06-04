/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

/* Type definitions for the C representation of variants */
#ifndef MLFI_VARIANT_H
#define MLFI_VARIANT_H

#include "mlfi_pricing.h"

/* String buffers */

typedef struct mlfi_strbuf {
  int len;
  int size;
  char *data;
} mlfi_strbuf;

mlfi_strbuf *mlfi_strbuf_create();
mlfi_strbuf *mlfi_strbuf_printf(mlfi_strbuf *, const char *, ...);
void mlfi_strbuf_ensure_space(mlfi_strbuf *, int);
void mlfi_strbuf_puts(mlfi_strbuf *, const char*);
void mlfi_strbuf_putc(mlfi_strbuf *, const char);
mlfi_strbuf *mlfi_strbuf_readfile(mlfi_strbuf *, const char *);
char *mlfi_strbuf_getfree(mlfi_strbuf *);

#define bprintf(...) mlfi_strbuf_getfree(mlfi_strbuf_printf(mlfi_strbuf_create(),__VA_ARGS__))

/* Variants */

typedef enum varkind {
  vunit, vbool, vint, vfloat, vstring, vdate,
  vtuple, vlist, varray, voption, vrecord, vconstructor,
  vvariant
} varkind;

struct variant;

typedef struct variant_field {
  char *name;
  struct variant *value;
} variant_field;

typedef struct variant_constructor {
  char *name;
  struct variant **arguments;
} variant_constructor;

typedef struct variant {
  varkind kind;
  int len; /* For record, tuple, list, array, constructor, option */
  union {
    int vbool;
    long vint;
    mlfi_date vdate;
    int vunit;
    double vfloat;
    char *vstring;
    struct variant **vtuple;
    struct variant **vlist;
    struct variant **varray;
    struct variant *voption;
    variant_field **vrecord;
    variant_constructor *vconstructor;
    struct variant *vvariant;
  };
} variant;


void mlfi_dump_variant(const variant *);
variant *mlfi_parse_variant(const char *); /* NULL if parse error */
void mlfi_print_variant(mlfi_strbuf *, const variant *);
char *mlfi_string_of_variant(const variant *);

/* Accessors */

int mlfi_variant_len(const variant *); /* For record, tuple, list, array, constructor, option */
variant *mlfi_variant_nth(const variant *, int, char *, int); /* For tuple, list, array, constructor */
variant *mlfi_variant_field(const variant *, char *, char *, int); /* For record */
variant_constructor *mlfi_variant_constructor(const variant *, char *, int); /* For constructor */
variant *mlfi_variant_kind(variant *, varkind, char *, char *, int);

int mlfi_variant_list(const variant *, char *, int);
int mlfi_variant_bool(const variant *, char *, int);
long mlfi_variant_int(const variant *, char *, int);
variant* mlfi_variant_variant(const variant *, char *, int);
mlfi_date mlfi_variant_date(const variant *, char *, int);
double mlfi_variant_float(const variant *, char *, int);
char *mlfi_variant_string(const variant *, char *, int);

double *mlfi_variant_float_vector(const variant*, int*, char*, int);
mlfi_date *mlfi_variant_date_vector(const variant*, int*, char*, int);
char **mlfi_variant_string_vector(const variant*, int*, char*, int);

#define _field(v,name) mlfi_variant_field(v,#name,__FILE__,__LINE__)
#define _variant(v) mlfi_variant_variant(v,__FILE__,__LINE__)
#define _int(v) mlfi_variant_int(v,__FILE__,__LINE__)
#define _bool(v) mlfi_variant_bool(v,__FILE__,__LINE__)
#define _string(v) mlfi_variant_string(v,__FILE__,__LINE__)
#define _date(v) mlfi_variant_date(v,__FILE__,__LINE__)
#define _float(v) mlfi_variant_float(v,__FILE__,__LINE__)
#define _constructor(v) mlfi_variant_constructor(v,__FILE__,__LINE__)
#define _record(v) mlfi_variant_kind(v,vrecord,"record",__FILE__,__LINE__)
#define _array(v) mlfi_variant_kind(v,varray,"array",__FILE__,__LINE__)
#define _list(v) mlfi_variant_kind(v,vlist,"list",__FILE__,__LINE__)
#define _tuple(v) mlfi_variant_kind(v,vtuple,"tuple",__FILE__,__LINE__)
#define _nth(v,i) mlfi_variant_nth(v,i,__FILE__,__LINE__)
#define _fst(v) _nth(v,0)
#define _snd(v) _nth(v,1)
#define _thd(v) _nth(v,2)
#define _float_array(v,len) mlfi_variant_float_vector(_array(v),len,__FILE__,__LINE__)
#define _string_array_len(v,len) (len = (-1), mlfi_variant_string_vector(_array(v),&len,__FILE__,__LINE__))
#define _date_array(v,len) mlfi_variant_date_vector(_array(v),len,__FILE__,__LINE__)
#define _date_array_len(v,len) (len = (-1), mlfi_variant_date_vector(_array(v),&len,__FILE__,__LINE__))


/* Constructors */

variant *mlfi_variant_create_tuple(varkind k, int);
variant *mlfi_variant_mk_tuple_lit(varkind k, variant *, ...);
variant *mlfi_variant_mk_tuple_arr(varkind k, int, variant **);
variant *mlfi_variant_mk_record_lit(char *, variant *, ...);
variant *mlfi_variant_mk_constructor(char *, variant *, ...);
variant *mlfi_variant_mk_option(variant *);
variant *mlfi_variant_mk_variant(variant *);

variant *mlfi_variant_mk_bool(int);
variant *mlfi_variant_mk_int(long);
variant *mlfi_variant_mk_date(mlfi_date);
variant *mlfi_variant_mk_float(double);
variant *mlfi_variant_mk_string(char*);

variant *mlfi_variant_mk_float_list(int, double*);
variant *mlfi_variant_mk_float_list_list(int, int, double**);

#define create_array(len) mlfi_variant_create_tuple(varray,len)
#define create_list(len) mlfi_variant_create_tuple(vlist,len)
#define create_tuple(len) mlfi_variant_create_tuple(vtuple,len)
#define bool(x) mlfi_variant_mk_bool(x)
#define int(x) mlfi_variant_mk_int(x)
#define date(x) mlfi_variant_mk_date(x)
#define float(x) mlfi_variant_mk_float(x)
#define string(x) mlfi_variant_mk_string(x)
#define stringf(...) mlfi_variant_mk_string(bprintf(__VA_ARGS__))
#define tuple(...) mlfi_variant_mk_tuple_lit(vtuple, __VA_ARGS__, NULL)
#define list(...) mlfi_variant_mk_tuple_lit(vlist, __VA_ARGS__, NULL)
#define array(...) mlfi_variant_mk_tuple_lit(varray, __VA_ARGS__, NULL)
#define record(...)  mlfi_variant_mk_record_lit(__VA_ARGS__, NULL)
#define tuple_arr(len,a) mlfi_variant_mk_tuple_arr(vtuple, len, a)
#define list_arr(len,a) mlfi_variant_mk_tuple_arr(vlist, len, a)
#define array_arr(len,a) mlfi_variant_mk_tuple_arr(varray, len, a)
#define elist() mlfi_variant_mk_tuple_arr(vlist, 0, NULL)
#define earray() mlfi_variant_mk_tuple_arr(vlist, 0, NULL)
#define constructor(...) mlfi_variant_mk_constructor(__VA_ARGS__, NULL)
#define option(x) mlfi_variant_mk_option(x)
#define some(x) mlfi_variant_mk_option(x)
#define none() mlfi_variant_mk_option(NULL)
#define mk_variant(x) mlfi_variant_mk_variant(x)
#define float_list(len, data) mlfi_variant_mk_float_list(len, data)
#define float_list_list(len1, len2, data) mlfi_variant_mk_float_list_list(len1, len2, data)

#define mlfi_err(...) (fprintf(stderr, "File %s, line %u:", file, line), fprintf(stderr, __VA_ARGS__), fprintf(stderr, "\n"), fflush(stderr), exit(2))

#endif
