#ifndef DATASET_PARSER
#define DATASET_PARSER

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>


struct array_reader {
    char* elems;
    int64_t n_elems_space;
    int64_t elem_size;
    int64_t n_elems_used;
    int64_t *shape;
    int (*elem_reader)(void*);
};

int peekc() {
    int c = getchar();
    ungetc(c,stdin);
    return c;
}

void skipspaces() {
    int c = getchar();
    if (isspace(c)) {
        skipspaces();
    } else if (c == '/') {
        // Skip to end of line.
        for (; c != '\n' && c != EOF; c = getchar())
           ;        
        skipspaces();  // Next line may have more spaces.
    } else if (c != EOF) {
        ungetc(c, stdin);
    }
}

int read_elem(struct array_reader *reader) {
    int ret;
    if (reader->n_elems_used == reader->n_elems_space) {
        reader->n_elems_space *= 2;
        reader->elems = static_cast<char*> (
          realloc(reader->elems,
                  reader->n_elems_space * reader->elem_size) );
    }

    ret = reader->elem_reader(reader->elems + reader->n_elems_used * reader->elem_size);

    if (ret == 0) {
        reader->n_elems_used++;
    }

    return ret;
}

int read_array_elems(struct array_reader *reader, int dims) {
    int c;
    int ret;
    int first = 1;
    char *knows_dimsize = static_cast<char*> ( calloc(dims,sizeof(char)) );
    int cur_dim = dims-1;
    int64_t *elems_read_in_dim = static_cast<int64_t*>( calloc(dims,sizeof(int64_t)) );
    while (1) {
        skipspaces();

        c = getchar();
        if (c == ']') {
          if (knows_dimsize[cur_dim]) {
            if (reader->shape[cur_dim] != elems_read_in_dim[cur_dim]) {
              ret = 1;
              break;
            }
          } else {
            knows_dimsize[cur_dim] = 1;
            reader->shape[cur_dim] = elems_read_in_dim[cur_dim];
          }
          if (cur_dim == 0) {
            ret = 0;
            break;
          } else {
            cur_dim--;
            elems_read_in_dim[cur_dim]++;
          }
        } else if (c == ',') {
          skipspaces();
          c = getchar();
          if (c == '[') {
            if (cur_dim == dims - 1) {
              ret = 1;
              break;
            }
            first = 1;
            cur_dim++;
            elems_read_in_dim[cur_dim] = 0;
          } else if (cur_dim == dims - 1) {
            ungetc(c, stdin);
            ret = read_elem(reader);
            if (ret != 0) {
              break;
            }
            elems_read_in_dim[cur_dim]++;
          } else {
            ret = 1;
            break;
          }
        } else if (c == EOF) {
          ret = 1;
          break;
        } else if (first) {
          if (c == '[') {
            if (cur_dim == dims - 1) {
              ret = 1;
              break;
            }
            cur_dim++;
            elems_read_in_dim[cur_dim] = 0;
          } else {
            ungetc(c, stdin);
            ret = read_elem(reader);
            if (ret != 0) {
              break;
            }
            elems_read_in_dim[cur_dim]++;
            first = 0;
          }
        } else {
          ret = 1;
          break;
        }
    }

    free(knows_dimsize);
    free(elems_read_in_dim);
    return ret;
}
////////////////////
/// Entry Points ///
////////////////////
int read_array(int64_t elem_size, int (*elem_reader)(void*),
               void **data, int64_t *shape, int64_t dims) {
    int ret;
    struct array_reader reader;
     int64_t read_dims = 0;
    while (1) {
        int c;
        skipspaces();
        c = getchar();
        if (c=='[') {
          read_dims++;
        } else {
          if (c != EOF) {
            ungetc(c, stdin);
          }
          break;
        }
    }

    if (read_dims != dims) {
//        printf("ERROR DIMS: %ld %ld \n\n", read_dims, dims);
        return 1;
    }

    reader.shape = shape;
    reader.n_elems_used = 0;
    reader.elem_size = elem_size;
    reader.n_elems_space = 16;
    reader.elems = static_cast<char*>( calloc(elem_size, reader.n_elems_space) );
    reader.elem_reader = elem_reader;

    ret = read_array_elems(&reader, dims);

    *data = reader.elems;

    return ret;
}

int read_int(void* dest) {
    skipspaces();
    if (scanf("%d", (int*)dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}

int read_char(void* dest) {
    skipspaces();
    if (scanf("%c", (char*)dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}

int read_float(void* dest) {
    skipspaces();
    if (scanf("%f", (float*)dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}

int read_double(void* dest) {
    skipspaces();
    if (scanf("%lf", (double*)dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}


///////////////////////
/// Writing Dataset ///
///////////////////////
void write_scal( const int   * i, const char* msg ) {
    fprintf(stdout, "%d ", *i);
    if( msg ) fprintf(stdout, "\t// %s\n", msg);
}
void write_scal( const double* r, const char* msg ) {
    fprintf(stdout, "%lf ", *r);
    if( msg ) fprintf(stdout, "\t// %s\n", msg);
}
void write_scal( const float * r, const char* msg ) {
    fprintf(stdout, "%f ", *r);
    if( msg ) fprintf(stdout, "\t// %s\n", msg);
}

template<class T>
void write_1Darr( const T* ptr, const int& N, const char* msg ) {
    fprintf(stdout, "\n [ ");
    for( int i = 0; i < N-1; i ++ ) {
        write_scal(&ptr[i], NULL);
        fprintf(stdout, ", ");
    }
    write_scal(&ptr[N-1], NULL);
    if (msg) fprintf(stdout, " ]\t//%s\n\n", msg);
    else     fprintf(stdout, " ]");
}

template<class T>
void write_2Darr( const T* ptr, const int& Nouter, const int& Ninner, const char* msg ) {
    fprintf(stdout, "\n[ ");
    for( int i = 0; i < Nouter-1; i ++ ) {
        write_1Darr( ptr + i*Ninner, Ninner, NULL );
        fprintf(stdout, ",");
    }
    write_1Darr( ptr + (Nouter-1)*Ninner, Ninner, NULL );
    if (msg) fprintf(stdout, "\n]\t//%s\n\n", msg);
    else     fprintf(stdout, "\n]\n");
}

#endif // DATASET_PARSER
