/***************************************************************************/
/*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               */
/*                                                                         */
/*  No part of this document may be reproduced or transmitted in any       */
/*  form or for any purpose without the express permission of LexiFi SAS.  */
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "mlfi_variant.h"

/* String buffer */

mlfi_strbuf *mlfi_strbuf_create() {
  mlfi_strbuf *buf = (mlfi_strbuf*) malloc(sizeof(mlfi_strbuf));
  assert(buf);
  buf->len = 0;
  buf->size = 10;
  buf->data = (char*) malloc(buf->size);
  assert(buf->data);
  memset(buf->data,0,buf->size);
  return buf;
}

void mlfi_strbuf_ensure_space(mlfi_strbuf *buf, int n) {
  /* Always keep an extra byte for the trailing 0 */
  if (buf->size <= buf->len + n) {
    int nsize = buf->size;
    char *ndata;
    do { nsize *= 2; } while (nsize <= buf->len + n);
    ndata = (char*) malloc(nsize);
    assert(ndata);
    memset(ndata,0,nsize);
    memcpy(ndata,buf->data,buf->len);
    free(buf->data);
    buf->data = ndata;
    buf->size = nsize;
  }
}

mlfi_strbuf *mlfi_strbuf_printf(mlfi_strbuf *buf, const char *format, ...) {
  int maxlen,len;
  va_list p;

 retry:
  maxlen = buf->size - buf->len;
  va_start(p, format);
  len = vsnprintf(buf->data + buf->len, maxlen, format, p);
  va_end(p);
  if (len >= maxlen || len < 0) {
    mlfi_strbuf_ensure_space(buf, (len >= 0 ? len + 10 : buf->size + 50));
    goto retry;
  }

  buf->len += len;
  return buf;
}

void mlfi_strbuf_puts(mlfi_strbuf *buf, const char *s) {
  int l = strlen(s);
  mlfi_strbuf_ensure_space(buf, l);
  memcpy(buf->data + buf->len, s, l);
  buf->len += l;
}

void mlfi_strbuf_putc(mlfi_strbuf *buf, const char c) {
  mlfi_strbuf_ensure_space(buf, 1);
  buf->data[buf->len++] = c;
}

char *mlfi_strbuf_getfree(mlfi_strbuf *buf) {
  char *s = buf->data;
  free(buf);
  return s;
}

mlfi_strbuf *mlfi_strbuf_readfile(mlfi_strbuf *buf, const char *filename) {
  FILE *file = fopen(filename, "r");
  size_t read;
  assert(file && filename);
  do {
    mlfi_strbuf_ensure_space(buf, 1024);
    read = fread(buf->data + buf->len, 1, 1024, file);
    buf->len += read;
  } while (read == 1024);
  fclose(file);
  return buf;
}


/* Parsing (string -> variant) */

/* This code assumes that all currency names have 3 chars. */

/* TODO:
   - more liberal syntax (optional ; at the end of sequences, extra whitespace,
   ...)
   - error messages
   - in case of error, free malloc'ed data
*/

#define is_idchar(c) (is_digit(c) || ('A' <= (c) && (c) <= 'Z') || (c) == '_' || ('a' <= (c) && (c) <= 'z') || ((c) == '\''))
#define is_digit(c) ((c) >= '0' && (c) <= '9')
#define is_hexdigit(c) (((c) >= '0' && (c) <= '9') || ((c) >= 'a' && (c) <= 'f') || ((c) >= 'A' && (c) <= 'F'))
#define TRUE 1
#define FALSE 0
#define try(f,r) if (!f(&c,&r)) return FALSE;
#define _char(d) if (*c != d) return FALSE; else c++;
#define _ischar(d) (*c == d ? (c++, TRUE) : FALSE)
#define _lit(s) parse_lit(&c,s)
#define __int(r) try(parse_int, r)
#define success(res) { *result = res; *cp = c; return TRUE; }
#define _parse_variant(r) try(parse_variant, r)
#define _variant_ws(r) try(parse_variant_ws, r)

static int parse_int(char **cp, long *result) {
  long res = 0;
  char *c = *cp;
  int sign = 1;

  if (_ischar('-')) sign = -1;
  if (!is_digit(*c)) return FALSE;

  while (is_digit(*c)) res = 10 * res + (*c++ - '0');
  success(sign * res);
}

static int hexdigit(char c) {
  if (is_digit(c)) return (c - '0');
  if ('a' <= c && c <= 'z') return (10 + c - 'a');
  if ('A' <= c && c <= 'Z') return (10 + c - 'A');
  assert(0); return 0;
}

static int parse_float(char **cp, double *result) {
  char *c = *cp;
  double f;
  if (*c == '-') c++;
  if (*c != '.') {
    long i;
    __int(i);
    if (*c != '.' && *c != 'e' && *c != 'E') return FALSE;
  }
  f = strtod(*cp,&c);
  if (c == *cp) return FALSE;
  success(f);
}

static int parse_lit(char **cp, char *s) {
  char *c = *cp;
  while (*s != 0) if (*c == 0 || *c++ != *s++) return FALSE;
  *cp = c;
  return TRUE;
}

static int parse_string(char **cp, char **result) {
  int len = 0;
  int bufsize = 128;
  char d;
  char *buf, *res;
  char *c = *cp;

  _char('\"');
  buf = (char*) malloc(bufsize);
  while ((d = *c++) != '\"') {
    if (d == 0) return FALSE;
    if (d == '\\') {
      d = *c++;
      if (d == 'n') d = '\n';
      else if (d == 'r') d = '\r';
      else if (d == 'b') d = '\b';
      else if (d == 't') d = '\t';
      else if (d == '\\' || d == '\'' || d == '"' || d == ' ') { }
      else if (d == '\n' || d == '\r') {
        if (d == '\r' && *c == '\n') c++;
        while (*c == ' ' || *c == '\t') c++;
        continue;
      } else if (is_digit(d) && is_digit(*c) && is_digit(*(c+1))) {
        d = 100 * (d - '0') + 10 * (*c - '0') + (*(c+1) - '0');
        c+=2;
      } else if (d == 'x' && is_hexdigit(*c) && is_hexdigit(*(c+1))) {
        d = (char) (16 * hexdigit(*c) + hexdigit(*(c+1)));
        c+=2;
      } else return FALSE; /* Illegal backslash escape */
    }
    if (len == bufsize) {
      char *newbuf;
      int newsize = 2 * bufsize;
      newbuf = (char*) malloc(newsize);
      memcpy(newbuf,buf,bufsize);
      free(buf);
      buf = newbuf;
      bufsize = newsize;
    }
    buf[len++] = d;
  }
  res = (char*) malloc(len+1);
  memcpy(res,buf,len);
  free(buf);
  res[len] = 0;
  success(res);
}

static int parse_date(char **cp, mlfi_date *result) {
  char *c = *cp;
  long y,m,d;
  long hr = 12, mn = 0, sec = 0;

  __int(y); _char('-'); __int(m); _char('-'); __int(d);
  if (*c == 'T') {
    c++;
    __int(hr); _char(':'); __int(mn);
    if (*c == ':') { c++; __int(sec); }
  }

  /* date_of_gregorian */
  success
    ((((m == 1 || m == 2) ?
       ( 1461 * ( y + 4800 - 1 ) ) / 4 +
       ( 367 * ( m + 10 ) ) / 12 -
       ( 3 * ( ( y + 4900 - 1 ) / 100 ) ) / 4
       :
       ( 1461 * ( y + 4800 ) ) / 4 +
       ( 367 * ( m - 2 ) ) / 12 -
       ( 3 * ( ( y + 4900 ) / 100 ) ) / 4)
      + d - 32075 - 2444238) * 1440
     + hr * 60 + mn);
}


#define mkvar_len(k,v,l)                                                     \
  { variant *res=(variant*) malloc(sizeof(variant)); res->kind = k; res->len = l; res->k = v; success(res); }
#define mkvar(k,v) mkvar_len(k,v,0)

typedef struct listcell {
  void *data;
  struct listcell *next;
} listcell;

static listcell *cons(void *data, listcell *next) {
  listcell *c = (listcell*) malloc(sizeof(listcell));
  c->data = data;
  c->next = next;
  return c;
}

static void **array_of_list(listcell *l, int n) {
  void **r = (void**) malloc(sizeof(void*) * n);
  r += n;
  while (n-- > 0) { listcell *next = l->next; r--; *r = l->data; free(l); l = next; }
  return r;
}

#define skip_ws() while (*c == ' ' || *c == '\t' || *c == '\r' || *c == '\n') c++

static int parse_variant(char **cp, variant **result);

static int parse_variant_ws(char **cp, variant **result) {
  if (parse_variant(cp, result)) {
    char *c = *cp; skip_ws(); *cp = c; return TRUE;
  }
  else return FALSE;
}

static int parse_variant(char **cp, variant **result) {
  char *c = *cp;
  skip_ws();

  { mlfi_date i; if (parse_date(&c, &i)) mkvar(vdate,i); }
  { double f; if (parse_float(&c, &f)) mkvar(vfloat,f); }
  { long i; if (parse_int(&c, &i)) mkvar(vint,i); }
  if (_lit("()")) mkvar(vunit,0);
  if (_lit("true")) mkvar(vbool,TRUE);
  if (_lit("false")) mkvar(vbool,FALSE);
  if (_lit("infinity")) mkvar(vfloat,1e100);/*1e100*/
  if (_lit("variant")) { variant *v; _parse_variant(v); mkvar(vvariant, v); }
  if (*c == '\"') { char *s; if (parse_string(&c, &s)) mkvar(vstring,s); }
  if (_ischar('(')) {
    variant *v;
    int nb;
    listcell *l;
    _variant_ws(v);
    if (_ischar(')')) success(v);
    l = cons((void*) v, NULL); nb = 1;
    while (_ischar(',')) { _variant_ws(v);  l = cons((void*) v, l); nb++; }
    _char(')');
    mkvar_len(vtuple, (variant**) array_of_list(l, nb), nb);
  }
  if (_ischar('[')) {
    variant *v;
    int nb = 0;
    listcell *l = NULL;
    int array = 0;
    if (_ischar('|')) {
      skip_ws();
      if (_lit("|]")) mkvar_len(varray, NULL, 0);
      array = 1;
   } else {
      skip_ws();
      if (_ischar(']')) mkvar_len(vlist, NULL, 0);
    }

    do { _variant_ws(v);  l = cons((void*) v, l); nb++; }  while (_ischar(';'));
    if (array) { _char('|'); }
    _char(']');
    if (array) mkvar_len(varray, (variant**) array_of_list(l, nb), nb);
    mkvar_len(vlist, (variant**) array_of_list(l, nb), nb);
  }
  if (_ischar('{')) {
    int nb = 0;
    listcell *l = NULL;

    skip_ws();
    if (_ischar('}')) mkvar_len(vrecord, NULL, 0);
    do {
      variant_field *f = (variant_field*) malloc(sizeof(variant_field));
      char *name_start;

      skip_ws();
      name_start = c;
      while (is_idchar(*c)) c++;
      f->name = (char*) malloc(c - name_start + 1);
      memcpy(f->name,name_start,c-name_start);
      f->name[c-name_start] = 0;

      skip_ws();
      _char('=');
      _variant_ws(f->value);
      l = cons((void*) f, l); nb++;
    } while (_ischar(';'));
    _char('}');
    mkvar_len(vrecord, (variant_field**) array_of_list(l,nb), nb);
  }
  if (_lit("Some")) {
    variant *v;
    _parse_variant(v);
    mkvar_len(voption, v, 1);
  }
  if (_lit("None")) mkvar_len(voption, NULL, 0);
  if ('A' <= *c && *c <= 'Z') {
    char *name_start = c;
    int len;
    int nb;
    variant_constructor *constr;
    variant **args;

    while (is_idchar(*c)) c++;
    len = c - name_start;
    /* TODO: all the possible characters in constructors */

    skip_ws();
    if (_ischar('(')) {
      variant *v;
      listcell *l;
      _variant_ws(v);
      l = cons((void*) v, NULL); nb = 1;
      while (_ischar(',')) { _variant_ws(v);  l = cons((void*) v, l); nb++; }
      _char(')');
      args = (variant**) array_of_list(l, nb);
    } else {
      args = NULL;
      nb = 0;
    }
    constr = (variant_constructor*) malloc (sizeof(*constr));
    constr->name = (char*) malloc (len + 1);
    memcpy(constr->name, name_start, len);
    constr->name[len] = 0;
    constr->arguments = args;
    mkvar_len(vconstructor, constr, nb);
  }
  return FALSE;
}

variant *mlfi_parse_variant(const char *s) {
  variant *r;
  char *c = (char*) s;
  if (!parse_variant_ws(&c, &r)) return NULL;
  if (*c != 0) return NULL;
  return r;
}

/* Dump the representation (for debugging) */

static void dump_variants(int len, variant **l) {
  int i;
  if (len == 0) { printf("()"); return; }
  printf("(");
  for (i = 0; i < len; i++) { if (0 < i) printf(";"); mlfi_dump_variant(l[i]); }
  printf(")");
}

static void dump_fields(int len, variant_field **l) {
  int i;
  if (len == 0) { printf("{}"); return; }
  printf("{");
  for (i = 0; i < len; i++) {
    if (0 < i) printf(";");
    printf(" %s:", l[i]->name);
    mlfi_dump_variant(l[i]->value);
  }
  printf(" }");
}

void mlfi_dump_variant(const variant *r) {
  if (r->kind == vint)
    printf("int %ld", r->vint);
  else if (r->kind == vfloat)
    printf("float %f", r->vfloat);
  else if (r->kind == vdate)
    printf("date %ld", r->vdate);
  else if (r->kind == vunit)
    printf("unit");
  else if (r->kind == vbool)
    printf("bool %i", r->vbool);
  else if (r->kind == vtuple)
    { printf("tuple"); dump_variants(r->len,r->vtuple); }
  else if (r->kind == vlist)
    { printf("list"); dump_variants(r->len,r->vlist); }
  else if (r->kind == varray)
    { printf("array"); dump_variants(r->len,r->varray); }
  else if (r->kind == voption && r->voption)
    { printf("some "); mlfi_dump_variant(r->voption); }
  else if (r->kind == voption && !r->voption)
    printf("none");
  else if (r->kind == vconstructor)
    { printf("constr %s", r->vconstructor->name); dump_variants(r->len,r->vconstructor->arguments); }
  else if (r->kind == vstring)
    printf("string %s", r->vstring);
  else if (r->kind == vrecord)
    { printf("record"); dump_fields(r->len,r->vrecord); }
  else if (r->kind == vvariant)
    { printf("variant "); mlfi_dump_variant(r->vvariant); }
  else
    printf("other (%i)", r->kind);
}

/* Pretty-print */

#define bputs(s) mlfi_strbuf_puts(buf, s)
#define bputc(c) mlfi_strbuf_putc(buf, c)
#define bprint(...) mlfi_strbuf_printf(buf, __VA_ARGS__)

static void print_variants(mlfi_strbuf *buf, char *open, char *sep, char *close, int len, variant **l) {
  int i;
  bputs(open);
  if (len)
    for (i = 0; i < len; i++) {
      if (i) bputs(sep);
      mlfi_print_variant(buf, l[i]);
    }
  bputs(close);
}

void mlfi_print_variant(mlfi_strbuf *buf, const variant *r) {
  if (r->kind == vint)
    bprint("%ld", r->vint);
  else if (r->kind == vfloat) {
    if (r->vfloat == r->vfloat)
      bprint("%#.16g", r->vfloat);
    else
      bputs("nan");
  }
  else if (r->kind == vdate) {
    int daytime = r->vdate % 1440;
    int jul = r->vdate / 1440;
    int hr = daytime / 60;
    int mn = daytime % 60;
    int l, n, i, j, d, m, y;
    l = jul + 68569 + 2444238;
    n = ( 4 * l ) / 146097;
    l = l - ( 146097 * n + 3 ) / 4;
    i = ( 4000 * ( l + 1 ) ) / 1461001;
    l = l - ( 1461 * i ) / 4 + 31;
    j = ( 80 * l ) / 2447;
    d = l - ( 2447 * j ) / 80;
    l = j / 11;
    m = j + 2 - ( 12 * l );
    y = 100 * ( n - 49 ) + i + l;
    bprint("%04i-%02i-%02iT%02i:%02i:00", y, m, d, hr, mn);
  }
  else if (r->kind == vunit)
    bputs("()");
  else if (r->kind == vbool && r->vbool)
    bputs("true");
  else if (r->kind == vbool && !r->vbool)
    bputs("false");
  else if (r->kind == vtuple)
    print_variants(buf,"(",", ",")",r->len,r->vtuple);
  else if (r->kind == vlist)
    print_variants(buf,"[","; ","]",r->len,r->vtuple);
  else if (r->kind == varray)
    print_variants(buf,"[|","; ","|]",r->len,r->vtuple);
  else if (r->kind == voption && r->voption)
    { bputs("(Some "); mlfi_print_variant(buf, r->voption); bputc(')'); }
  else if (r->kind == voption && !r->voption)
    bputs("None");
  else if (r->kind == vconstructor && !r->len)
    bputs(r->vconstructor->name);
  else if (r->kind == vconstructor && r->len)
    { bprint("(%s", r->vconstructor->name);
      print_variants(buf,"(",", ",")",r->len,r->vconstructor->arguments);
      bputc(')');
    }
  else if (r->kind == vstring)
    {
      char *c = r->vstring;
      bputc('"');
      while (*c)  {
        if (*c == '\"' || *c == '\\') { bputc('\\'); bputc(*c++); }
        else if (*c == '\n') { bputc('\\'); bputc('n'); c++; }
        else if (*c == '\r') { bputc('\\'); bputc('r'); c++; }
        else if (*c == '\t') { bputc('\\'); bputc('t'); c++; }
        else if (isprint(*c)) bputc(*c++);
        else bprint("%03d",*c++);
      }
      bputc('"');
    }
  else if (r->kind == vrecord) {
    int i = 0;
    bputs("{");
    if (r->len)
      for (i = 0; i < r->len; i++) {
        if (i) bputc(';');
        bprint("%s=", r->vrecord[i]->name);
        mlfi_print_variant(buf,r->vrecord[i]->value);
      }
    bputs("}");
  } else if (r->kind == vvariant)
    { bputs("(variant "); mlfi_print_variant(buf, r->vvariant); bputs(")"); }
  else {
    bprintf("ERROR %i", r->kind);
    printf("ERROR %i", r->kind);
  }
}

#undef bputs
#undef bputc
#undef bprint

/* Accessors */

void mlfi_print_variant_stderr(const variant *v) {
  mlfi_strbuf *buf = mlfi_strbuf_create();
  mlfi_print_variant(buf,v);
  fprintf(stderr, "%s\n", buf->data);
  fflush(stderr);
}

variant *mlfi_variant_nth(const variant *v, int n, char *file, int line) {
  if (!(0 <= n && n < v->len)) mlfi_err("variant does not have enough elements");
  if (v->kind == vconstructor)
    return v->vconstructor->arguments[n];
  if (!(v->kind == vtuple || v->kind == varray || v->kind == vlist))
    mlfi_err("variant does not have sub elements");
  return v->vtuple[n];
}

variant *mlfi_variant_kind(variant *v, varkind k, char *s, char *file, int line) {
  if (v->kind != k) {
    mlfi_print_variant_stderr(v);
    mlfi_err("variant is not a %s", s);
  }

  return v;
}

variant *mlfi_variant_field(const variant *v, char *name, char *file, int line) {
  /* TODO: sort fields + dichotomic lookup? */
  int i;
  if (v->kind != vrecord) {
    mlfi_print_variant_stderr(v);
    mlfi_err("variant is not a record");
  }
  for (i = 0; i < v->len; i++) {
    if (strcmp(v->vrecord[i]->name,name) == 0) return v->vrecord[i]->value;
  }
  mlfi_err("field %s not present in record", name);
  exit(2);
}

variant_constructor *mlfi_variant_constructor(const variant *v, char *file, int line) {
  if (v->kind != vconstructor) mlfi_err("variant is not a constructor");
  return v->vconstructor;
}

int mlfi_variant_len(const variant *v) {
  return v->len;
}

#define extr(ty,na,k,e) \
  ty mlfi_variant_##na(const variant *v, char *file, int line) { \
    if (v->kind !=k) { mlfi_print_variant_stderr(v); mlfi_err("variant is not a " #na); } \
    return v->e; }
#define mk(ty,na,k,l) \
variant *mlfi_variant_mk_##na(ty v) \
{ variant *res=(variant*) malloc(sizeof(variant)); res->kind = k; res->len = l; res->k = v; return res; }
#define basic(ty,na,k) extr(ty,na,k,k) mk(ty,na,k,0)

basic(int, bool, vbool);
basic(long, int, vint);
basic(mlfi_date, date, vdate);
basic(double, float, vfloat);
basic(char*, string, vstring);
extr(variant*, variant, vvariant, vvariant);
extr(char*, contructor, vconstructor, vconstructor->name);

double *mlfi_variant_float_vector(const variant *v, int *len, char *file, int line) {
  double *r;
  int i;
  if (len) {
    if (*len < 0)
      *len = v->len;
    else
      if (*len != v->len)
        mlfi_err("Wrong array size %d (should be %d)", v->len, *len);
  }
  r = (double*) malloc(sizeof(double) * v->len);
  for (i=0; i < v->len; i++) r[i] = _float(v->varray[i]);
  return r;
}

char **mlfi_variant_string_vector(const variant *v, int *len, char *file, int line) {
  char **r;
  int i;
  if (len) {
    if (*len < 0)
      *len = v->len;
    else
      if (*len != v->len)
        mlfi_err("Wrong array size %d (should be %d)", v->len, *len);
  }
  r = (char**) malloc(sizeof(char*) * v->len);
  for (i=0; i < v->len; i++) r[i] = _string(v->varray[i]);
  return r;
}

mlfi_date *mlfi_variant_date_vector(const variant *v, int *len, char *file, int line) {
  mlfi_date *r;
  int i;
  if (len) {
    if (*len < 0)
      *len = v->len;
    else
      if (*len != v->len)
        mlfi_err("Wrong array size %d (should be %d)", v->len, *len);
  }
  if (!v->len) return NULL;
  r = (mlfi_date*) malloc(sizeof(mlfi_date) * v->len);
  for (i=0; i < v->len; i++) r[i] = _date(v->varray[i]);
  return r;
}

/* Constructors */

mk(variant*, option, voption, (v ? 1 : 0));
mk(variant*, variant, vvariant, 1);

static variant **mkarray(variant *v, va_list p, va_list q, int *len) {
  variant *x;
  variant **res;
  int n;

  x = v;
  n = 0;
  while (x != NULL) { n++; x = va_arg(p, variant*); }
  va_end(p);
  *len = n;

  if (!n) return NULL;

  res = (variant**) malloc(sizeof(void*) * n);
  x = v;
  n = 0;
  while (x != NULL) { res[n] = x; n++; x = va_arg(q, variant*); }
  va_end(q);

  return res;
}

variant *mlfi_variant_mk_tuple_arr(varkind k, int len, variant **a) {
  variant *res=(variant*) malloc(sizeof(variant));
  res->kind = k;
  res->vtuple = a;
  res->len = len;
  return res;
}

variant *mlfi_variant_create_tuple(varkind k, int len) {
  variant **a = (variant**) malloc(sizeof(variant) * len);
  return mlfi_variant_mk_tuple_arr(k,len,a);
}

variant *mlfi_variant_mk_tuple_lit(varkind k, variant *v, ...) {
  int len;
  variant **a;
  va_list p,q;
  va_start(p, v); va_start(q, v);
  a = mkarray(v, p, q, &len);
  return mlfi_variant_mk_tuple_arr(k, len, a);
}

variant *mlfi_variant_mk_constructor(char *name, variant *v, ...) {
  va_list p, q;
  variant *res=(variant*) malloc(sizeof(variant));
  variant_constructor *constr=(variant_constructor*) malloc(sizeof(variant_constructor));
  res->kind = vconstructor;
  res->vconstructor = constr;
  constr->name = name;
  va_start(p, v); va_start(q, v); constr->arguments = mkarray(v, p, q, &res->len);
  return res;
}

variant *mlfi_variant_mk_record_lit(char *name0, variant *v0, ...) {
  variant *v;
  char *name;
  int n = 0;
  va_list p;
  variant *res=(variant*) malloc(sizeof(variant));

  res->kind = vrecord;
  name = name0;
  v = v0;
  va_start(p, v0);
  while (1) { n++; name = va_arg(p, char*); if (!name) break; v = va_arg(p, variant*); }
  va_end(p);
  res->len = n;

  if (!n) {
    res->vrecord = NULL;
  } else {
    res->vrecord = (variant_field**) malloc(sizeof(variant_field*) * n);
    name = name0;
    v = v0;
    va_start(p, v0);
    for (n = 0; n < res->len; n++) {
      variant_field *f = (variant_field*) malloc(sizeof(variant_field));
      res->vrecord[n] = f;
      if (n > 0) { name = va_arg(p, char*); v = va_arg(p, variant*); }
      f->name = name;
      f->value = v;
    }
  }

  return res;
}


char *mlfi_string_of_variant(const variant *v) {
  mlfi_strbuf *buf = mlfi_strbuf_create();
  mlfi_print_variant(buf, v);
  return buf->data;
}

variant *mlfi_variant_mk_float_list(int len, double *data) {
  variant *v = create_list(len);
  int i;
  for (i = 0; i < len; i++)
    v->vlist[i] = float(data[i]);
  return v;
}

variant *mlfi_variant_mk_float_list_list(int len1, int len2, double **data) {
  variant *v = create_list(len1);
  int i;
  for (i = 0; i < len1; i++)
    v->vlist[i] = mlfi_variant_mk_float_list(len2, data[i]);
  return v;
}

#undef is_idchar
#undef is_digit
#undef is_hexdigit
#undef TRUE
#undef FALSE
#undef try
#undef _char
#undef _ischar
#undef _lit
#undef __int
#undef success
#undef _parse_variant
#undef mkvar_len
#undef mkvar
#undef skip_ws
#undef mlfi_err
