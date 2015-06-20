#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
struct tuple_int_mem_int {
    int elem_0;
    unsigned char* elem_1;
    int elem_2;
};
static struct tuple_int_mem_int futhark_main(int dir_vs_mem_size_27,
                                             int all_md_c_mem_size_71,
                                             int all_md_vols_mem_size_255,
                                             int all_md_drifts_mem_size_345,
                                             int all_md_st_mem_size_354,
                                             int all_md_detvals_mem_size_360,
                                             int all_md_disc_mem_size_363,
                                             int bb_inds_mem_size_369,
                                             int bb_data_mem_size_376,
                                             unsigned char* dir_vs_mem_40,
                                             unsigned char* all_md_c_mem_84,
                                             unsigned char* all_md_vols_mem_268,
                                             unsigned char* all_md_drifts_mem_350,
                                             unsigned char* all_md_st_mem_359,
                                             unsigned char* all_md_detvals_mem_361,
                                             unsigned char* all_md_disc_mem_368,
                                             unsigned char* bb_inds_mem_371,
                                             unsigned char* bb_data_mem_381,
                                             int size_14, int size_18,
                                             int size_19, int size_20,
                                             int size_21, int size_22,
                                             int size_86, int size_89,
                                             int size_90, int size_91,
                                             int size_92, int size_93,
                                             int size_94, int size_95,
                                             int size_111, int size_112,
                                             int size_113, int size_114,
                                             int size_115, int size_116,
                                             int size_117,
                                             int contract_number_118,
                                             int num_mc_it_119,
                                             int num_dates_120, int num_und_121,
                                             int num_models_122,
                                             int num_bits_123);
float futhark_toFloat32(int x)
{
    return x;
}
int futhark_trunc32(float x)
{
    return x;
}
float futhark_log32(float x)
{
    return log(x);
}
float futhark_sqrt32(float x)
{
    return sqrt(x);
}
float futhark_exp32(float x)
{
    return exp(x);
}
double futhark_toFloat64(int x)
{
    return x;
}
int futhark_trunc64(double x)
{
    return x;
}
double futhark_log64(double x)
{
    return log(x);
}
double futhark_sqrt64(double x)
{
    return sqrt(x);
}
double futhark_exp64(double x)
{
    return exp(x);
}
int timeval_subtract(struct timeval* result, struct timeval* t2,
                     struct timeval* t1)
{
    unsigned int resolution = 1000000;
    long diff = t2->tv_usec + resolution * t2->tv_sec - (t1->tv_usec +
                                                         resolution *
                                                         t1->tv_sec);
    
    result->tv_sec = diff / resolution;
    result->tv_usec = diff % resolution;
    return diff < 0;
}
struct array_reader {
    char* elems;
    int64_t n_elems_space;
    int64_t elem_size;
    int64_t n_elems_used;
    int64_t* shape;
    int(* elem_reader)(void*);
};
int peekc()

{
    int c = getchar();
    
    ungetc(c, stdin);
    return c;
}
void skipspaces()

{
    int c = getchar();
    
    if (isspace(c)) {
        skipspaces();
    } else if (c == '/') {
        // Skip to end of line.
        for (; c != '\n' && c != EOF; c = getchar())
            ;
        // Next line may have more spaces.
        skipspaces();
    } else if (c != EOF) {
        ungetc(c, stdin);
    }
}
int read_elem(struct array_reader* reader)
{
    int ret;
    
    if (reader->n_elems_used == reader->n_elems_space) {
        reader->n_elems_space *= 2;
        reader->elems = realloc(reader->elems, reader->n_elems_space *
                                reader->elem_size);
    }
    ret = reader->elem_reader(reader->elems + reader->n_elems_used *
        reader->elem_size);
    if (ret == 0) {
        reader->n_elems_used++;
    }
    return ret;
}
int read_array_elems(struct array_reader* reader, int dims)
{
    int c;
    int ret;
    int first = 1;
    char* knows_dimsize = calloc(dims, sizeof(char));
    int cur_dim = dims - 1;
    int64_t* elems_read_in_dim = calloc(dims, sizeof(int64_t));
    
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
int read_array(int64_t elem_size, int(* elem_reader)(void*), void** data,
               int64_t* shape, int64_t dims)
{
    int ret;
    struct array_reader reader;
    int64_t read_dims = 0;
    
    while (1) {
        int c;
        
        skipspaces();
        c = getchar();
        if (c == '[') {
            read_dims++;
        } else {
            if (c != EOF) {
                ungetc(c, stdin);
            }
            break;
        }
    }
    if (read_dims != dims) {
        return 1;
    }
    reader.shape = shape;
    reader.n_elems_used = 0;
    reader.elem_size = elem_size;
    reader.n_elems_space = 16;
    reader.elems = calloc(elem_size, reader.n_elems_space);
    reader.elem_reader = elem_reader;
    ret = read_array_elems(&reader, dims);
    *data = reader.elems;
    return ret;
}
int read_int(void* dest)
{
    skipspaces();
    if (scanf("%d", (int*) dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}
int read_char(void* dest)
{
    skipspaces();
    if (scanf("%c", (char*) dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}
int read_double(void* dest)
{
    skipspaces();
    if (scanf("%lf", (double*) dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}
int read_float(void* dest)
{
    skipspaces();
    if (scanf("%f", (float*) dest) == 1) {
        return 0;
    } else {
        return 1;
    }
}
int read_bool(void* dest)
{
    char b[4];
    
    skipspaces();
    if (scanf("%4c", b) == 1) {
        if (strncmp(b, "True", 4) == 0) {
            *(int*) dest = 1;
            return 0;
        } else if (strncmp(b, "Fals", 4) == 0 && getchar() == 'e') {
            *(int*) dest = 0;
            return 0;
        } else {
            return 1;
        }
    } else {
        return 1;
    }
}
static struct tuple_int_mem_int futhark_main(int dir_vs_mem_size_27,
                                             int all_md_c_mem_size_71,
                                             int all_md_vols_mem_size_255,
                                             int all_md_drifts_mem_size_345,
                                             int all_md_st_mem_size_354,
                                             int all_md_detvals_mem_size_360,
                                             int all_md_disc_mem_size_363,
                                             int bb_inds_mem_size_369,
                                             int bb_data_mem_size_376,
                                             unsigned char* dir_vs_mem_40,
                                             unsigned char* all_md_c_mem_84,
                                             unsigned char* all_md_vols_mem_268,
                                             unsigned char* all_md_drifts_mem_350,
                                             unsigned char* all_md_st_mem_359,
                                             unsigned char* all_md_detvals_mem_361,
                                             unsigned char* all_md_disc_mem_368,
                                             unsigned char* bb_inds_mem_371,
                                             unsigned char* bb_data_mem_381,
                                             int size_14, int size_18,
                                             int size_19, int size_20,
                                             int size_21, int size_22,
                                             int size_86, int size_89,
                                             int size_90, int size_91,
                                             int size_92, int size_93,
                                             int size_94, int size_95,
                                             int size_111, int size_112,
                                             int size_113, int size_114,
                                             int size_115, int size_116,
                                             int size_117,
                                             int contract_number_118,
                                             int num_mc_it_119,
                                             int num_dates_120, int num_und_121,
                                             int num_models_122,
                                             int num_bits_123)
{
    int out_memsize_24;
    unsigned char* out_mem_23;
    int out_arrsize_26;
    char y_133;
    
    y_133 = 0 < size_19;
    
    int bounds_check_134;
    
    if (!y_133) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_133",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:253:16-253:24");
        abort();
    }
    
    char y_135;
    
    y_135 = 0 < size_22;
    
    int bounds_check_136;
    
    if (!y_135) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_135",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:254:19-254:30");
        abort();
    }
    
    char y_138;
    
    y_138 = 0 < size_90;
    
    int bounds_check_139;
    
    if (!y_138) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_138",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:255:21-255:34");
        abort();
    }
    
    char y_141;
    
    y_141 = 0 < size_93;
    
    int bounds_check_142;
    
    if (!y_141) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_141",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:256:17-256:26");
        abort();
    }
    
    char y_144;
    
    y_144 = 0 < size_95;
    
    int bounds_check_145;
    
    if (!y_144) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_144",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:257:22-257:36");
        abort();
    }
    
    char y_146;
    
    y_146 = 0 < size_112;
    
    int bounds_check_147;
    
    if (!y_146) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_146",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:258:19-258:30");
        abort();
    }
    
    char cond_149;
    
    cond_149 = num_mc_it_119 == 0;
    
    double y_150;
    
    y_150 = futhark_toFloat64(num_bits_123);
    
    double res_151;
    
    res_151 = powl(2.0, y_150);
    
    int bytes_385;
    
    bytes_385 = 4 * num_bits_123;
    
    unsigned char* mem_387;
    
    mem_387 = malloc(bytes_385);
    for (int i_39 = 0; i_39 < num_bits_123; i_39++) {
        *(int*) &mem_387[i_39 * 1 * 4 + 0] = i_39;
    }
    
    char cond_153;
    
    cond_153 = num_und_121 == 0;
    
    int size_154;
    
    if (cond_153) {
        size_154 = 0;
    } else {
        size_154 = num_dates_120;
    }
    
    int size_155;
    
    if (cond_149) {
        size_155 = 0;
    } else {
        size_155 = size_154;
    }
    
    int y_156;
    
    y_156 = num_dates_120 * num_und_121;
    
    char assert_arg_157;
    
    assert_arg_157 = size_14 == y_156;
    
    int shape_ok_158;
    
    if (!assert_arg_157) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_157",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:190:20-190:27");
        abort();
    }
    
    char y_159;
    
    y_159 = 0 < size_114;
    
    int bounds_check_160;
    
    if (!y_159) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_159",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:147:14-147:21");
        abort();
    }
    
    char y_161;
    
    y_161 = 1 < size_114;
    
    int bounds_check_162;
    
    if (!y_161) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_161",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:148:14-148:21");
        abort();
    }
    
    char y_163;
    
    y_163 = 2 < size_114;
    
    int bounds_check_164;
    
    if (!y_163) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_163",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:149:14-149:21");
        abort();
    }
    
    char y_165;
    
    y_165 = 0 < size_116;
    
    int bounds_check_166;
    
    if (!y_165) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_165",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:150:14-150:21");
        abort();
    }
    
    char y_167;
    
    y_167 = 1 < size_116;
    
    int bounds_check_168;
    
    if (!y_167) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_167",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:151:14-151:21");
        abort();
    }
    
    char y_169;
    
    y_169 = 2 < size_116;
    
    int bounds_check_170;
    
    if (!y_169) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_169",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:152:14-152:21");
        abort();
    }
    
    int bytes_399;
    
    bytes_399 = 8 * num_dates_120;
    
    unsigned char* mem_404;
    
    mem_404 = malloc(bytes_399);
    for (int i_56 = 0; i_56 < num_dates_120; i_56++) {
        *(double*) &mem_404[i_56 * 1 * 8 + 0] = 0.0;
    }
    
    char y_172;
    
    y_172 = 0 < size_115;
    
    int bounds_check_173;
    
    if (!y_172) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_172",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:156:16-156:18");
        abort();
    }
    
    int x_174;
    
    x_174 = *(int*) &bb_inds_mem_371[(0 * size_115 + 0 * 1) * 4 + 0];
    
    int idx_175;
    
    idx_175 = x_174 - 1;
    
    char y_176;
    
    y_176 = 0 < size_117;
    
    int bounds_check_177;
    
    if (!y_176) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_176",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:156:28-156:30");
        abort();
    }
    
    double x_178;
    
    x_178 = *(double*) &bb_data_mem_381[(0 * size_117 + 0 * 1) * 8 + 0];
    
    char y_179;
    
    y_179 = 0 < num_dates_120;
    
    int bounds_check_180;
    
    if (!y_179) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "y_179",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:156:36-156:41");
        abort();
    }
    
    char x_181;
    
    x_181 = 0 <= idx_175;
    
    char y_182;
    
    y_182 = idx_175 < num_dates_120;
    
    char assert_arg_183;
    
    assert_arg_183 = x_181 && y_182;
    
    int bounds_check_184;
    
    if (!assert_arg_183) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_183",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:156:5-156:8");
        abort();
    }
    
    int bound_185;
    
    bound_185 = num_dates_120 - 1;
    
    char assert_arg_186;
    
    assert_arg_186 = size_155 == num_dates_120;
    
    int shape_cert_187;
    
    if (!assert_arg_186) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_186",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:261:27-261:41");
        abort();
    }
    
    char cond_188;
    
    cond_188 = size_155 == 0;
    
    int size_189;
    
    if (cond_188) {
        size_189 = 0;
    } else {
        size_189 = size_92;
    }
    
    int size_190;
    
    if (cond_149) {
        size_190 = 0;
    } else {
        size_190 = size_189;
    }
    
    char zip_cmp_192;
    
    zip_cmp_192 = size_155 == size_86;
    
    int zip_assert_193;
    
    if (!zip_cmp_192) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_192",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:218:38-218:41");
        abort();
    }
    
    char zip_cmp_195;
    
    zip_cmp_195 = size_155 == size_91;
    
    int zip_assert_196;
    
    if (!zip_cmp_195) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_195",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:218:38-218:41");
        abort();
    }
    
    char zip_cmp_198;
    
    zip_cmp_198 = num_und_121 == size_89;
    
    int zip_assert_199;
    
    if (!zip_cmp_198) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_198",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:214:33-214:36");
        abort();
    }
    
    char zip_cmp_200;
    
    zip_cmp_200 = size_92 == num_und_121;
    
    int zip_assert_201;
    
    if (!zip_cmp_200) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_200",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:214:13-214:16");
        abort();
    }
    
    char assert_arg_202;
    
    assert_arg_202 = size_92 == size_94;
    
    int shape_cert_203;
    
    if (!assert_arg_202) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_202",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:220:9-220:13");
        abort();
    }
    
    char assert_arg_205;
    
    assert_arg_205 = size_190 == size_92;
    
    int shape_cert_206;
    
    if (!assert_arg_205) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_205",
                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:262:27-262:39");
        abort();
    }
    
    char cond_207;
    
    cond_207 = contract_number_118 == 1;
    
    char cond_208;
    
    cond_208 = contract_number_118 == 2;
    
    char cond_209;
    
    cond_209 = contract_number_118 == 3;
    
    char y_210;
    
    y_210 = 0 < size_111;
    
    char y_211;
    
    y_211 = 0 < size_155;
    
    char y_212;
    
    y_212 = 0 < size_190;
    
    char y_213;
    
    y_213 = 0 < size_113;
    
    char zip_cmp_214;
    
    zip_cmp_214 = size_190 == 3;
    
    char y_215;
    
    y_215 = 1 < size_155;
    
    char y_216;
    
    y_216 = 2 < size_155;
    
    char y_217;
    
    y_217 = 3 < size_155;
    
    char y_218;
    
    y_218 = 4 < size_155;
    
    char y_219;
    
    y_219 = 1 < size_113;
    
    char y_220;
    
    y_220 = 2 < size_113;
    
    char y_221;
    
    y_221 = 3 < size_113;
    
    char y_222;
    
    y_222 = 4 < size_113;
    
    char y_223;
    
    y_223 = 1 < size_190;
    
    char y_224;
    
    y_224 = 2 < size_190;
    
    char y_281;
    
    y_281 = 366 < size_155;
    
    int x_282;
    
    x_282 = 366 * size_190;
    
    int i_283;
    
    i_283 = x_282 / size_92 - (x_282 % size_92 != 0 && x_282 % size_92 < 0 !=
                               size_92 < 0 ? 1 : 0);
    
    int y_284;
    
    y_284 = i_283 * size_92;
    
    int i_remnant_288;
    
    i_remnant_288 = x_282 - y_284;
    
    int flat_index_289;
    
    flat_index_289 = x_282 + 2;
    
    int i_291;
    
    i_291 = flat_index_289 / size_92 - (flat_index_289 % size_92 != 0 &&
                                        flat_index_289 % size_92 < 0 !=
                                        size_92 < 0 ? 1 : 0);
    
    int y_295;
    
    y_295 = i_291 * size_92;
    
    int i_remnant_296;
    
    i_remnant_296 = flat_index_289 - y_295;
    
    int flat_index_297;
    
    flat_index_297 = x_282 + 1;
    
    int i_298;
    
    i_298 = flat_index_297 / size_92 - (flat_index_297 % size_92 != 0 &&
                                        flat_index_297 % size_92 < 0 !=
                                        size_92 < 0 ? 1 : 0);
    
    int y_299;
    
    y_299 = i_298 * size_92;
    
    int i_remnant_300;
    
    i_remnant_300 = flat_index_297 - y_299;
    
    int new_slice_size_286;
    
    new_slice_size_286 = num_und_121 * size_18;
    
    unsigned char* mem_406;
    
    mem_406 = malloc(bytes_385);
    
    unsigned char* mem_418;
    
    mem_418 = malloc(bytes_385);
    
    int x_424;
    
    x_424 = 8 * num_und_121;
    
    int bytes_423;
    
    bytes_423 = x_424 * num_dates_120;
    
    unsigned char* mem_425;
    
    mem_425 = malloc(bytes_423);
    
    unsigned char* mem_433;
    
    mem_433 = malloc(bytes_399);
    
    unsigned char* mem_483;
    
    mem_483 = malloc(bytes_399);
    
    int x_524;
    
    x_524 = 8 * size_155;
    
    int bytes_523;
    
    bytes_523 = x_524 * size_92;
    
    unsigned char* mem_525;
    
    mem_525 = malloc(bytes_523);
    
    int x_535;
    
    x_535 = 8 * size_92;
    
    int bytes_534;
    
    bytes_534 = x_535 * size_155;
    
    unsigned char* mem_536;
    
    mem_536 = malloc(bytes_534);
    
    unsigned char* mem_540;
    
    mem_540 = malloc(x_524);
    
    unsigned char* mem_552;
    
    mem_552 = malloc(24);
    
    unsigned char* mem_554;
    
    mem_554 = malloc(bytes_385);
    
    unsigned char* mem_560;
    
    mem_560 = malloc(bytes_385);
    
    unsigned char* mem_563;
    
    mem_563 = malloc(bytes_423);
    
    unsigned char* mem_567;
    
    mem_567 = malloc(bytes_399);
    
    unsigned char* mem_573;
    
    mem_573 = malloc(bytes_399);
    
    unsigned char* mem_596;
    
    mem_596 = malloc(bytes_523);
    
    unsigned char* mem_607;
    
    mem_607 = malloc(bytes_534);
    
    unsigned char* mem_611;
    
    mem_611 = malloc(x_524);
    
    unsigned char* mem_623;
    
    mem_623 = malloc(x_524);
    
    int bytes_626;
    
    bytes_626 = 8 * size_190;
    
    unsigned char* mem_627;
    
    mem_627 = malloc(bytes_626);
    
    unsigned char* mem_640;
    
    mem_640 = malloc(bytes_385);
    
    unsigned char* mem_649;
    
    mem_649 = malloc(bytes_385);
    
    unsigned char* mem_652;
    
    mem_652 = malloc(bytes_423);
    
    unsigned char* mem_656;
    
    mem_656 = malloc(bytes_399);
    
    unsigned char* mem_662;
    
    mem_662 = malloc(bytes_399);
    
    unsigned char* mem_722;
    
    mem_722 = malloc(bytes_523);
    
    unsigned char* mem_771;
    
    mem_771 = malloc(bytes_534);
    
    unsigned char* mem_777;
    
    mem_777 = malloc(x_524);
    
    double res_2307;
    
    if (cond_207) {
        int bounds_check_301;
        
        if (!y_210) {
            fprintf(stderr, "Assertion %s at %s failed.\n", "y_210",
                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:284:16-284:26");
            abort();
        }
        
        double res_326;
        
        res_326 = *(double*) &all_md_detvals_mem_361[(0 * size_111 + 0 * 1) *
                                                     8 + 0];
        
        int bounds_check_329;
        
        if (!y_211) {
            fprintf(stderr, "Assertion %s at %s failed.\n", "y_211",
                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:285:17-285:20");
            abort();
        }
        
        int bounds_check_330;
        
        if (!y_212) {
            fprintf(stderr, "Assertion %s at %s failed.\n", "y_212",
                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:285:17-285:20");
            abort();
        }
        
        int bounds_check_331;
        
        if (!y_213) {
            fprintf(stderr, "Assertion %s at %s failed.\n", "y_213",
                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
            abort();
        }
        
        double y_332;
        
        y_332 = *(double*) &all_md_disc_mem_368[(0 * size_113 + 0 * 1) * 8 + 0];
        
        double res_989;
        double acc_48;
        
        acc_48 = 0.0;
        for (int i_47 = 0; i_47 < num_mc_it_119; i_47++) {
            int res_388;
            
            res_388 = i_47 + 1;
            
            int x_389;
            
            x_389 = res_388 >> 1;
            
            int res_390;
            
            res_390 = x_389 ^ res_388;
            for (int i_0 = 0; i_0 < num_bits_123; i_0++) {
                int res_392;
                
                res_392 = 1 << i_0;
                
                int x_393;
                
                x_393 = res_390 & res_392;
                
                char res_394;
                
                res_394 = x_393 == res_392;
                
                int partition_equivalence_class_395;
                
                if (res_394) {
                    partition_equivalence_class_395 = 0;
                } else {
                    partition_equivalence_class_395 = 1;
                }
                *(int*) &mem_406[i_0 * 1 * 4 + 0] =
                    partition_equivalence_class_395;
            }
            
            int filter_size_397;
            
            filter_size_397 = 0;
            
            int eqclass_104;
            
            for (int i_103 = 0; i_103 < num_bits_123; i_103++) {
                eqclass_104 = *(int*) &mem_406[i_103 * 1 * 4 + 0];
                if (eqclass_104 == 0) {
                    filter_size_397 = filter_size_397 + 1;
                } else { }
            }
            
            int partition_offset_109;
            
            partition_offset_109 = 0;
            
            int partition_cur_offset_110;
            
            for (int i_103 = 0; i_103 < num_bits_123; i_103++) {
                eqclass_104 = *(int*) &mem_406[i_103 * 1 * 4 + 0];
                if (eqclass_104 == 0) {
                    partition_cur_offset_110 = partition_offset_109;
                    *(int*) &mem_418[partition_cur_offset_110 * 1 * 4 + 0] =
                        *(int*) &mem_387[i_103 * 1 * 4 + 0];
                    partition_offset_109 = partition_offset_109 + 1;
                } else { }
            }
            for (int i_11 = 0; i_11 < num_und_121; i_11++) {
                int y_290;
                
                y_290 = i_11 * size_18;
                for (int i_7 = 0; i_7 < num_dates_120; i_7++) {
                    int x_287;
                    
                    x_287 = i_7 * new_slice_size_286;
                    
                    int x_292;
                    
                    x_292 = x_287 + y_290;
                    
                    int res_443;
                    int acc_5;
                    
                    acc_5 = 0;
                    for (int i_4 = 0; i_4 < filter_size_397; i_4++) {
                        int not_curried_438;
                        
                        not_curried_438 = *(int*) &mem_418[i_4 * 1 * 4 + 0];
                        
                        char y_439;
                        
                        y_439 = not_curried_438 < size_18;
                        
                        int bounds_check_440;
                        
                        if (!y_439) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_439",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:6:35-6:38");
                            abort();
                        }
                        
                        int flat_index_294;
                        
                        flat_index_294 = x_292 + not_curried_438;
                        
                        int i_302;
                        
                        i_302 = flat_index_294 / size_18 - (flat_index_294 %
                                                            size_18 != 0 &&
                                                            flat_index_294 %
                                                            size_18 < 0 !=
                                                            size_18 <
                                                            0 ? 1 : 0);
                        
                        int y_303;
                        
                        y_303 = i_302 * size_18;
                        
                        int i_remnant_304;
                        
                        i_remnant_304 = flat_index_294 - y_303;
                        
                        int res_441;
                        
                        res_441 = *(int*) &dir_vs_mem_40[(i_302 * size_18 +
                                                          i_remnant_304 * 1) *
                                                         4 + 0];
                        
                        int res_442;
                        
                        res_442 = acc_5 ^ res_441;
                        
                        int acc_tmp_226;
                        
                        acc_tmp_226 = res_442;
                        acc_5 = acc_tmp_226;
                    }
                    res_443 = acc_5;
                    
                    double x_444;
                    
                    x_444 = futhark_toFloat64(res_443);
                    
                    double res_445;
                    
                    res_445 = x_444 / res_151;
                    
                    double res_446;
                    
                    res_446 = res_445 - 0.5;
                    
                    char x_447;
                    
                    x_447 = res_446 < 0.0;
                    
                    double x_448;
                    
                    x_448 = 0.0 - res_446;
                    
                    char y_449;
                    
                    y_449 = x_448 <= 0.425;
                    
                    char x_450;
                    
                    x_450 = x_447 && y_449;
                    
                    char x_451;
                    
                    x_451 = 0.0 <= res_446;
                    
                    char y_452;
                    
                    y_452 = res_446 <= 0.425;
                    
                    char y_453;
                    
                    y_453 = x_451 && y_452;
                    
                    char cond_454;
                    
                    cond_454 = x_450 || y_453;
                    
                    double res_455;
                    
                    res_455 = res_446 + 0.5;
                    
                    double res_456;
                    
                    res_456 = 0.5 - res_446;
                    
                    double res_457;
                    
                    if (x_447) {
                        res_457 = res_455;
                    } else {
                        res_457 = res_456;
                    }
                    
                    double y_458;
                    
                    y_458 = res_446 * res_446;
                    
                    double arg_459;
                    
                    arg_459 = 0.180625 - y_458;
                    
                    double x_460;
                    
                    x_460 = arg_459 * 2509.0809287301227;
                    
                    double y_461;
                    
                    y_461 = x_460 + 33430.57558358813;
                    
                    double x_462;
                    
                    x_462 = arg_459 * y_461;
                    
                    double y_463;
                    
                    y_463 = x_462 + 67265.7709270087;
                    
                    double x_464;
                    
                    x_464 = arg_459 * y_463;
                    
                    double y_465;
                    
                    y_465 = x_464 + 45921.95393154987;
                    
                    double x_466;
                    
                    x_466 = arg_459 * y_465;
                    
                    double y_467;
                    
                    y_467 = x_466 + 13731.69376550946;
                    
                    double x_468;
                    
                    x_468 = arg_459 * y_467;
                    
                    double y_469;
                    
                    y_469 = x_468 + 1971.5909503065513;
                    
                    double x_470;
                    
                    x_470 = arg_459 * y_469;
                    
                    double y_471;
                    
                    y_471 = x_470 + 133.14166789178438;
                    
                    double x_472;
                    
                    x_472 = arg_459 * y_471;
                    
                    double x_473;
                    
                    x_473 = x_472 + 3.3871328727963665;
                    
                    double x_474;
                    
                    x_474 = arg_459 * 5226.495278852854;
                    
                    double y_475;
                    
                    y_475 = x_474 + 28729.085735721943;
                    
                    double x_476;
                    
                    x_476 = arg_459 * y_475;
                    
                    double y_496;
                    
                    y_496 = x_476 + 39307.89580009271;
                    
                    double x_497;
                    
                    x_497 = arg_459 * y_496;
                    
                    double y_498;
                    
                    y_498 = x_497 + 21213.794301586597;
                    
                    double x_499;
                    
                    x_499 = arg_459 * y_498;
                    
                    double y_500;
                    
                    y_500 = x_499 + 5394.196021424751;
                    
                    double x_501;
                    
                    x_501 = arg_459 * y_500;
                    
                    double y_502;
                    
                    y_502 = x_501 + 687.1870074920579;
                    
                    double x_503;
                    
                    x_503 = arg_459 * y_502;
                    
                    double y_504;
                    
                    y_504 = x_503 + 42.31333070160091;
                    
                    double x_505;
                    
                    x_505 = arg_459 * y_504;
                    
                    double y_506;
                    
                    y_506 = x_505 + 1.0;
                    
                    double res_715;
                    
                    if (cond_454) {
                        double res_507;
                        
                        res_507 = x_473 / y_506;
                        
                        double res_508;
                        
                        res_508 = res_446 * res_507;
                        res_715 = res_508;
                    } else {
                        double negate_arg_509;
                        
                        negate_arg_509 = futhark_log64(res_457);
                        
                        double arg_510;
                        
                        arg_510 = -negate_arg_509;
                        
                        double res_511;
                        
                        res_511 = futhark_sqrt64(arg_510);
                        
                        char cond_512;
                        
                        cond_512 = res_511 <= 5.0;
                        
                        double arg_513;
                        
                        arg_513 = res_511 - 1.6;
                        
                        double x_514;
                        
                        x_514 = arg_513 * 7.745450142783414e-4;
                        
                        double y_515;
                        
                        y_515 = x_514 + 2.2723844989269184e-2;
                        
                        double x_516;
                        
                        x_516 = arg_513 * y_515;
                        
                        double y_517;
                        
                        y_517 = x_516 + 0.2417807251774506;
                        
                        double x_518;
                        
                        x_518 = arg_513 * y_517;
                        
                        double y_575;
                        
                        y_575 = x_518 + 1.2704582524523684;
                        
                        double x_576;
                        
                        x_576 = arg_513 * y_575;
                        
                        double y_577;
                        
                        y_577 = x_576 + 3.6478483247632045;
                        
                        double x_578;
                        
                        x_578 = arg_513 * y_577;
                        
                        double y_579;
                        
                        y_579 = x_578 + 5.769497221460691;
                        
                        double x_580;
                        
                        x_580 = arg_513 * y_579;
                        
                        double y_581;
                        
                        y_581 = x_580 + 4.630337846156546;
                        
                        double x_582;
                        
                        x_582 = arg_513 * y_581;
                        
                        double x_583;
                        
                        x_583 = x_582 + 1.4234371107496835;
                        
                        double x_633;
                        
                        x_633 = arg_513 * 1.0507500716444169e-9;
                        
                        double y_636;
                        
                        y_636 = x_633 + 5.475938084995345e-4;
                        
                        double x_637;
                        
                        x_637 = arg_513 * y_636;
                        
                        double y_638;
                        
                        y_638 = x_637 + 1.5198666563616457e-2;
                        
                        double x_639;
                        
                        x_639 = arg_513 * y_638;
                        
                        double y_642;
                        
                        y_642 = x_639 + 0.14810397642748008;
                        
                        double x_645;
                        
                        x_645 = arg_513 * y_642;
                        
                        double y_647;
                        
                        y_647 = x_645 + 0.6897673349851;
                        
                        double x_675;
                        
                        x_675 = arg_513 * y_647;
                        
                        double y_676;
                        
                        y_676 = x_675 + 1.6763848301838038;
                        
                        double x_677;
                        
                        x_677 = arg_513 * y_676;
                        
                        double y_678;
                        
                        y_678 = x_677 + 2.053191626637759;
                        
                        double x_679;
                        
                        x_679 = arg_513 * y_678;
                        
                        double y_680;
                        
                        y_680 = x_679 + 1.0;
                        
                        double arg_681;
                        
                        arg_681 = res_511 - 5.0;
                        
                        double x_682;
                        
                        x_682 = arg_681 * 2.0103343992922881e-7;
                        
                        double y_683;
                        
                        y_683 = x_682 + 2.7115555687434876e-5;
                        
                        double x_684;
                        
                        x_684 = arg_681 * y_683;
                        
                        double y_685;
                        
                        y_685 = x_684 + 1.2426609473880784e-3;
                        
                        double x_686;
                        
                        x_686 = arg_681 * y_685;
                        
                        double y_687;
                        
                        y_687 = x_686 + 2.6532189526576124e-2;
                        
                        double x_688;
                        
                        x_688 = arg_681 * y_687;
                        
                        double y_689;
                        
                        y_689 = x_688 + 0.29656057182850487;
                        
                        double x_690;
                        
                        x_690 = arg_681 * y_689;
                        
                        double y_691;
                        
                        y_691 = x_690 + 1.7848265399172913;
                        
                        double x_692;
                        
                        x_692 = arg_681 * y_691;
                        
                        double y_693;
                        
                        y_693 = x_692 + 5.463784911164114;
                        
                        double x_694;
                        
                        x_694 = arg_681 * y_693;
                        
                        double x_695;
                        
                        x_695 = x_694 + 6.657904643501103;
                        
                        double x_696;
                        
                        x_696 = arg_681 * 2.0442631033899398e-5;
                        
                        double y_697;
                        
                        y_697 = x_696 + 1.421511758316446e-7;
                        
                        double x_698;
                        
                        x_698 = arg_681 * y_697;
                        
                        double y_699;
                        
                        y_699 = x_698 + 1.8463183175100548e-5;
                        
                        double x_700;
                        
                        x_700 = arg_681 * y_699;
                        
                        double y_701;
                        
                        y_701 = x_700 + 7.868691311456133e-4;
                        
                        double x_702;
                        
                        x_702 = arg_681 * y_701;
                        
                        double y_703;
                        
                        y_703 = x_702 + 1.4875361290850615e-2;
                        
                        double x_704;
                        
                        x_704 = arg_681 * y_703;
                        
                        double y_705;
                        
                        y_705 = x_704 + 0.1369298809227358;
                        
                        double x_706;
                        
                        x_706 = arg_681 * y_705;
                        
                        double y_707;
                        
                        y_707 = x_706 + 0.599832206555888;
                        
                        double x_708;
                        
                        x_708 = arg_681 * y_707;
                        
                        double y_709;
                        
                        y_709 = x_708 + 1.0;
                        
                        double res_712;
                        
                        if (cond_512) {
                            double res_710;
                            
                            res_710 = x_583 / y_680;
                            res_712 = res_710;
                        } else {
                            double res_711;
                            
                            res_711 = x_695 / y_709;
                            res_712 = res_711;
                        }
                        
                        double res_713;
                        
                        res_713 = 0.0 - res_712;
                        
                        double res_714;
                        
                        if (x_447) {
                            res_714 = res_713;
                        } else {
                            res_714 = res_712;
                        }
                        res_715 = res_714;
                    }
                    *(double*) &mem_433[i_7 * 1 * 8 + 0] = res_715;
                }
                memmove(mem_483 + 0, mem_404 + 0, num_dates_120 *
                        sizeof(double));
                
                double y_719;
                
                y_719 = *(double*) &mem_433[0 * 1 * 8 + 0];
                
                double lw_val_720;
                
                lw_val_720 = x_178 * y_719;
                *(double*) &mem_483[idx_175 * 1 * 8 + 0] = lw_val_720;
                for (int i_251 = 0; i_251 < num_dates_120; i_251++) {
                    *(double*) &mem_425[8 * (i_11 * num_dates_120) + 8 *
                                        i_251] = *(double*) &mem_483[8 * i_251];
                }
                for (int i_723 = 0; i_723 < bound_185; i_723++) {
                    int i_724;
                    
                    i_724 = i_723 + 1;
                    
                    char y_725;
                    
                    y_725 = i_724 < size_115;
                    
                    int bounds_check_727;
                    
                    if (!y_725) {
                        fprintf(stderr, "Assertion %s at %s failed.\n", "y_725",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:160:22-160:24");
                        abort();
                    }
                    
                    int x_728;
                    
                    x_728 = *(int*) &bb_inds_mem_371[(1 * size_115 + i_724 *
                                                      1) * 4 + 0];
                    
                    int res_729;
                    
                    res_729 = x_728 - 1;
                    
                    int x_730;
                    
                    x_730 = *(int*) &bb_inds_mem_371[(2 * size_115 + i_724 *
                                                      1) * 4 + 0];
                    
                    int res_731;
                    
                    res_731 = x_730 - 1;
                    
                    int x_736;
                    
                    x_736 = *(int*) &bb_inds_mem_371[(0 * size_115 + i_724 *
                                                      1) * 4 + 0];
                    
                    int res_737;
                    
                    res_737 = x_736 - 1;
                    
                    char x_738;
                    
                    x_738 = 0 <= res_731;
                    
                    char y_739;
                    
                    y_739 = res_731 < num_dates_120;
                    
                    char assert_arg_740;
                    
                    assert_arg_740 = x_738 && y_739;
                    
                    int bounds_check_741;
                    
                    if (!assert_arg_740) {
                        fprintf(stderr, "Assertion %s at %s failed.\n",
                                "assert_arg_740",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:164:22-164:27");
                        abort();
                    }
                    
                    double res_742;
                    
                    res_742 = *(double*) &mem_425[(i_11 * num_dates_120 +
                                                   res_731 * 1) * 8 + 0];
                    
                    double res_743;
                    
                    res_743 = *(double*) &mem_433[i_724 * 1 * 8 + 0];
                    
                    char y_744;
                    
                    y_744 = i_724 < size_117;
                    
                    int bounds_check_745;
                    
                    if (!y_744) {
                        fprintf(stderr, "Assertion %s at %s failed.\n", "y_744",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:166:22-166:24");
                        abort();
                    }
                    
                    double x_746;
                    
                    x_746 = *(double*) &bb_data_mem_381[(2 * size_117 + i_724 *
                                                         1) * 8 + 0];
                    
                    double x_747;
                    
                    x_747 = x_746 * res_742;
                    
                    double x_748;
                    
                    x_748 = *(double*) &bb_data_mem_381[(0 * size_117 + i_724 *
                                                         1) * 8 + 0];
                    
                    double y_749;
                    
                    y_749 = x_748 * res_743;
                    
                    double res_750;
                    
                    res_750 = x_747 + y_749;
                    
                    int x_751;
                    
                    x_751 = res_729 + 1;
                    
                    char cond_752;
                    
                    cond_752 = x_751 == 0;
                    
                    char x_753;
                    
                    x_753 = 0 <= res_729;
                    
                    char y_754;
                    
                    y_754 = res_729 < num_dates_120;
                    
                    char assert_arg_755;
                    
                    assert_arg_755 = x_753 && y_754;
                    
                    double lw_val_761;
                    
                    if (cond_752) {
                        lw_val_761 = res_750;
                    } else {
                        double x_756;
                        
                        x_756 = *(double*) &bb_data_mem_381[(1 * size_117 +
                                                             i_724 * 1) * 8 +
                                                            0];
                        
                        int bounds_check_757;
                        
                        if (!assert_arg_755) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "assert_arg_755",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:170:51-170:56");
                            abort();
                        }
                        
                        double y_758;
                        
                        y_758 = *(double*) &mem_425[(i_11 * num_dates_120 +
                                                     res_729 * 1) * 8 + 0];
                        
                        double y_759;
                        
                        y_759 = x_756 * y_758;
                        
                        double res_760;
                        
                        res_760 = res_750 + y_759;
                        lw_val_761 = res_760;
                    }
                    
                    char x_762;
                    
                    x_762 = 0 <= res_737;
                    
                    char y_763;
                    
                    y_763 = res_737 < num_dates_120;
                    
                    char assert_arg_764;
                    
                    assert_arg_764 = x_762 && y_763;
                    
                    int bounds_check_768;
                    
                    if (!assert_arg_764) {
                        fprintf(stderr, "Assertion %s at %s failed.\n",
                                "assert_arg_764",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:168:13-168:16");
                        abort();
                    }
                    *(double*) &mem_425[(i_11 * num_dates_120 + res_737 * 1) *
                                        8 + 0] = lw_val_761;
                }
                for (int ii_778 = 0; ii_778 < bound_185; ii_778++) {
                    int y_779;
                    
                    y_779 = ii_778 + 1;
                    
                    int res_780;
                    
                    res_780 = num_dates_120 - y_779;
                    
                    double x_781;
                    
                    x_781 = *(double*) &mem_425[(i_11 * num_dates_120 +
                                                 res_780 * 1) * 8 + 0];
                    
                    int i_782;
                    
                    i_782 = res_780 - 1;
                    
                    double y_783;
                    
                    y_783 = *(double*) &mem_425[(i_11 * num_dates_120 + i_782 *
                                                 1) * 8 + 0];
                    
                    double lw_val_784;
                    
                    lw_val_784 = x_781 - y_783;
                    *(double*) &mem_425[(i_11 * num_dates_120 + res_780 * 1) *
                                        8 + 0] = lw_val_784;
                }
            }
            for (int i_29 = 0; i_29 < size_155; i_29++) {
                int x_308;
                
                x_308 = i_29 * size_92;
                
                int x_317;
                
                x_317 = i_29 * size_89;
                
                int x_143;
                
                x_143 = i_29 * num_und_121;
                for (int i_25 = 0; i_25 < size_92; i_25++) {
                    int flat_index_310;
                    
                    flat_index_310 = x_308 + i_25;
                    
                    int i_311;
                    
                    i_311 = flat_index_310 / size_92 - (flat_index_310 %
                                                        size_92 != 0 &&
                                                        flat_index_310 %
                                                        size_92 < 0 != size_92 <
                                                        0 ? 1 : 0);
                    
                    int y_312;
                    
                    y_312 = i_311 * size_92;
                    
                    int i_remnant_313;
                    
                    i_remnant_313 = flat_index_310 - y_312;
                    
                    double x_942;
                    
                    x_942 = *(double*) &all_md_drifts_mem_350[(0 * (size_91 *
                                                                    size_92) +
                                                               i_311 * size_92 +
                                                               i_remnant_313 *
                                                               1) * 8 + 0];
                    
                    int flat_index_319;
                    
                    flat_index_319 = x_317 + i_25;
                    
                    int i_320;
                    
                    i_320 = flat_index_319 / size_89 - (flat_index_319 %
                                                        size_89 != 0 &&
                                                        flat_index_319 %
                                                        size_89 < 0 != size_89 <
                                                        0 ? 1 : 0);
                    
                    int y_321;
                    
                    y_321 = i_320 * size_89;
                    
                    int i_remnant_322;
                    
                    i_remnant_322 = flat_index_319 - y_321;
                    
                    double y_943;
                    
                    y_943 = *(double*) &all_md_vols_mem_268[(0 * (size_86 *
                                                                  size_89) +
                                                             i_320 * size_89 +
                                                             i_remnant_322 *
                                                             1) * 8 + 0];
                    
                    int arg_944;
                    
                    arg_944 = i_25 + 1;
                    
                    char y_948;
                    
                    y_948 = i_25 < size_20;
                    
                    int bounds_check_949;
                    
                    if (!y_948) {
                        fprintf(stderr, "Assertion %s at %s failed.\n", "y_948",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:205:58-205:62");
                        abort();
                    }
                    
                    char split_index_cnd_951;
                    
                    split_index_cnd_951 = arg_944 <= size_21;
                    
                    int split_index_assert_952;
                    
                    if (!split_index_cnd_951) {
                        fprintf(stderr, "Assertion %s at %s failed.\n",
                                "split_index_cnd_951",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:201:56-201:61");
                        abort();
                    }
                    
                    double res_965;
                    double acc_17;
                    
                    acc_17 = 0.0;
                    for (int i_16 = 0; i_16 < arg_944; i_16++) {
                        int flat_index_194;
                        
                        flat_index_194 = x_143 + i_16;
                        
                        int i_197;
                        
                        i_197 = flat_index_194 / num_und_121 - (flat_index_194 %
                                                                num_und_121 !=
                                                                0 &&
                                                                flat_index_194 %
                                                                num_und_121 <
                                                                0 !=
                                                                num_und_121 <
                                                                0 ? 1 : 0);
                        
                        int y_204;
                        
                        y_204 = i_197 * num_und_121;
                        
                        int i_remnant_225;
                        
                        i_remnant_225 = flat_index_194 - y_204;
                        
                        double x_961;
                        
                        x_961 = *(double*) &mem_425[(i_remnant_225 *
                                                     num_dates_120 + i_197 *
                                                     1) * 8 + 0];
                        
                        double x_962;
                        
                        x_962 = *(double*) &all_md_c_mem_84[(0 * (size_20 *
                                                                  size_21) +
                                                             i_25 * size_21 +
                                                             i_16 * 1) * 8 + 0];
                        
                        double res_963;
                        
                        res_963 = x_961 * x_962;
                        
                        double res_964;
                        
                        res_964 = acc_17 + res_963;
                        
                        double acc_tmp_293;
                        
                        acc_tmp_293 = res_964;
                        acc_17 = acc_tmp_293;
                    }
                    res_965 = acc_17;
                    
                    double res_966;
                    
                    res_966 = res_965 * y_943;
                    
                    double res_967;
                    
                    res_967 = x_942 + res_966;
                    *(double*) &mem_525[(i_29 * size_92 + i_25 * 1) * 8 + 0] =
                        res_967;
                }
            }
            for (int i_43 = 0; i_43 < size_92; i_43++) {
                double x_973;
                
                x_973 = *(double*) &all_md_st_mem_359[(0 * size_94 + i_43 * 1) *
                                                      8 + 0];
                for (int i_33 = 0; i_33 < size_155; i_33++) {
                    double not_curried_974;
                    
                    not_curried_974 = *(double*) &mem_525[(i_33 * size_92 +
                                                           i_43 * 1) * 8 + 0];
                    
                    double curried_fun_result_975;
                    
                    curried_fun_result_975 = futhark_exp64(not_curried_974);
                    *(double*) &mem_540[i_33 * 1 * 8 + 0] =
                        curried_fun_result_975;
                }
                
                double acc_38;
                
                acc_38 = x_973;
                for (int i_37 = 0; i_37 < size_155; i_37++) {
                    double y_978;
                    
                    y_978 = *(double*) &mem_540[i_37 * 1 * 8 + 0];
                    
                    double res_979;
                    
                    res_979 = acc_38 * y_978;
                    *(double*) &mem_536[(i_43 * size_155 + i_37 * 1) * 8 + 0] =
                        res_979;
                    
                    double row_42;
                    
                    row_42 = *(double*) &mem_536[(i_43 * size_155 + i_37 * 1) *
                                                 8 + 0];
                    
                    double acc_tmp_314;
                    
                    acc_tmp_314 = row_42;
                    acc_38 = acc_tmp_314;
                }
            }
            
            double x_982;
            
            x_982 = *(double*) &mem_536[(0 * size_155 + 0 * 1) * 8 + 0];
            
            double x_983;
            
            x_983 = x_982 - 4000.0;
            
            double res_984;
            
            res_984 = x_983 * res_326;
            
            char cond_985;
            
            cond_985 = res_984 < 0.0;
            
            double res_986;
            
            if (cond_985) {
                res_986 = 0.0;
            } else {
                res_986 = res_984;
            }
            
            double res_987;
            
            res_987 = res_986 * y_332;
            
            double res_988;
            
            res_988 = acc_48 + res_987;
            
            double acc_tmp_70;
            
            acc_tmp_70 = res_988;
            acc_48 = acc_tmp_70;
        }
        res_989 = acc_48;
        res_2307 = res_989;
    } else {
        double res_2306;
        
        if (cond_208) {
            *(double*) &mem_552[0 * 1 * 8 + 0] = 2.6609544843735446e-4;
            *(double*) &mem_552[1 * 1 * 8 + 0] = 8.445945945945946e-5;
            *(double*) &mem_552[2 * 1 * 8 + 0] = 8.333333333333334e-4;
            
            int zip_assert_991;
            
            if (!zip_cmp_214) {
                fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_214",
                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:292:56-292:59");
                abort();
            }
            
            int bounds_check_993;
            
            if (!y_212) {
                fprintf(stderr, "Assertion %s at %s failed.\n", "y_212",
                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:329:62-329:65");
                abort();
            }
            
            int bounds_check_994;
            
            if (!y_211) {
                fprintf(stderr, "Assertion %s at %s failed.\n", "y_211",
                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:294:20-294:24");
                abort();
            }
            
            double res_1845;
            double acc_232;
            
            acc_232 = 0.0;
            for (int i_231 = 0; i_231 < num_mc_it_119; i_231++) {
                int res_1000;
                
                res_1000 = i_231 + 1;
                
                int x_1001;
                
                x_1001 = res_1000 >> 1;
                
                int res_1002;
                
                res_1002 = x_1001 ^ res_1000;
                for (int i_50 = 0; i_50 < num_bits_123; i_50++) {
                    int res_1004;
                    
                    res_1004 = 1 << i_50;
                    
                    int x_1005;
                    
                    x_1005 = res_1002 & res_1004;
                    
                    char res_1006;
                    
                    res_1006 = x_1005 == res_1004;
                    
                    int partition_equivalence_class_1007;
                    
                    if (res_1006) {
                        partition_equivalence_class_1007 = 0;
                    } else {
                        partition_equivalence_class_1007 = 1;
                    }
                    *(int*) &mem_554[i_50 * 1 * 4 + 0] =
                        partition_equivalence_class_1007;
                }
                
                int filter_size_1009;
                
                filter_size_1009 = 0;
                
                int eqclass_341;
                
                for (int i_336 = 0; i_336 < num_bits_123; i_336++) {
                    eqclass_341 = *(int*) &mem_554[i_336 * 1 * 4 + 0];
                    if (eqclass_341 == 0) {
                        filter_size_1009 = filter_size_1009 + 1;
                    } else { }
                }
                
                int partition_offset_370;
                
                partition_offset_370 = 0;
                
                int partition_cur_offset_372;
                
                for (int i_336 = 0; i_336 < num_bits_123; i_336++) {
                    eqclass_341 = *(int*) &mem_554[i_336 * 1 * 4 + 0];
                    if (eqclass_341 == 0) {
                        partition_cur_offset_372 = partition_offset_370;
                        *(int*) &mem_560[partition_cur_offset_372 * 1 * 4 + 0] =
                            *(int*) &mem_387[i_336 * 1 * 4 + 0];
                        partition_offset_370 = partition_offset_370 + 1;
                    } else { }
                }
                for (int i_61 = 0; i_61 < num_und_121; i_61++) {
                    int y_334;
                    
                    y_334 = i_61 * size_18;
                    for (int i_57 = 0; i_57 < num_dates_120; i_57++) {
                        int x_333;
                        
                        x_333 = i_57 * new_slice_size_286;
                        
                        int x_335;
                        
                        x_335 = x_333 + y_334;
                        
                        int res_1055;
                        int acc_55;
                        
                        acc_55 = 0;
                        for (int i_54 = 0; i_54 < filter_size_1009; i_54++) {
                            int not_curried_1050;
                            
                            not_curried_1050 = *(int*) &mem_560[i_54 * 1 * 4 +
                                                                0];
                            
                            char y_1051;
                            
                            y_1051 = not_curried_1050 < size_18;
                            
                            int bounds_check_1052;
                            
                            if (!y_1051) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_1051",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:6:35-6:38");
                                abort();
                            }
                            
                            int flat_index_337;
                            
                            flat_index_337 = x_335 + not_curried_1050;
                            
                            int i_338;
                            
                            i_338 = flat_index_337 / size_18 - (flat_index_337 %
                                                                size_18 != 0 &&
                                                                flat_index_337 %
                                                                size_18 < 0 !=
                                                                size_18 <
                                                                0 ? 1 : 0);
                            
                            int y_339;
                            
                            y_339 = i_338 * size_18;
                            
                            int i_remnant_340;
                            
                            i_remnant_340 = flat_index_337 - y_339;
                            
                            int res_1053;
                            
                            res_1053 = *(int*) &dir_vs_mem_40[(i_338 * size_18 +
                                                               i_remnant_340 *
                                                               1) * 4 + 0];
                            
                            int res_1054;
                            
                            res_1054 = acc_55 ^ res_1053;
                            
                            int acc_tmp_413;
                            
                            acc_tmp_413 = res_1054;
                            acc_55 = acc_tmp_413;
                        }
                        res_1055 = acc_55;
                        
                        double x_1056;
                        
                        x_1056 = futhark_toFloat64(res_1055);
                        
                        double res_1057;
                        
                        res_1057 = x_1056 / res_151;
                        
                        double res_1058;
                        
                        res_1058 = res_1057 - 0.5;
                        
                        char x_1059;
                        
                        x_1059 = res_1058 < 0.0;
                        
                        double x_1060;
                        
                        x_1060 = 0.0 - res_1058;
                        
                        char y_1061;
                        
                        y_1061 = x_1060 <= 0.425;
                        
                        char x_1062;
                        
                        x_1062 = x_1059 && y_1061;
                        
                        char x_1063;
                        
                        x_1063 = 0.0 <= res_1058;
                        
                        char y_1064;
                        
                        y_1064 = res_1058 <= 0.425;
                        
                        char y_1065;
                        
                        y_1065 = x_1063 && y_1064;
                        
                        char cond_1066;
                        
                        cond_1066 = x_1062 || y_1065;
                        
                        double res_1067;
                        
                        res_1067 = res_1058 + 0.5;
                        
                        double res_1068;
                        
                        res_1068 = 0.5 - res_1058;
                        
                        double res_1069;
                        
                        if (x_1059) {
                            res_1069 = res_1067;
                        } else {
                            res_1069 = res_1068;
                        }
                        
                        double y_1070;
                        
                        y_1070 = res_1058 * res_1058;
                        
                        double arg_1071;
                        
                        arg_1071 = 0.180625 - y_1070;
                        
                        double x_1072;
                        
                        x_1072 = arg_1071 * 2509.0809287301227;
                        
                        double y_1073;
                        
                        y_1073 = x_1072 + 33430.57558358813;
                        
                        double x_1074;
                        
                        x_1074 = arg_1071 * y_1073;
                        
                        double y_1075;
                        
                        y_1075 = x_1074 + 67265.7709270087;
                        
                        double x_1076;
                        
                        x_1076 = arg_1071 * y_1075;
                        
                        double y_1077;
                        
                        y_1077 = x_1076 + 45921.95393154987;
                        
                        double x_1079;
                        
                        x_1079 = arg_1071 * y_1077;
                        
                        double y_1081;
                        
                        y_1081 = x_1079 + 13731.69376550946;
                        
                        double x_1082;
                        
                        x_1082 = arg_1071 * y_1081;
                        
                        double y_1083;
                        
                        y_1083 = x_1082 + 1971.5909503065513;
                        
                        double x_1084;
                        
                        x_1084 = arg_1071 * y_1083;
                        
                        double y_1085;
                        
                        y_1085 = x_1084 + 133.14166789178438;
                        
                        double x_1086;
                        
                        x_1086 = arg_1071 * y_1085;
                        
                        double x_1087;
                        
                        x_1087 = x_1086 + 3.3871328727963665;
                        
                        double x_1088;
                        
                        x_1088 = arg_1071 * 5226.495278852854;
                        
                        double y_1089;
                        
                        y_1089 = x_1088 + 28729.085735721943;
                        
                        double x_1090;
                        
                        x_1090 = arg_1071 * y_1089;
                        
                        double y_1096;
                        
                        y_1096 = x_1090 + 39307.89580009271;
                        
                        double x_1097;
                        
                        x_1097 = arg_1071 * y_1096;
                        
                        double y_1098;
                        
                        y_1098 = x_1097 + 21213.794301586597;
                        
                        double x_1099;
                        
                        x_1099 = arg_1071 * y_1098;
                        
                        double y_1100;
                        
                        y_1100 = x_1099 + 5394.196021424751;
                        
                        double x_1101;
                        
                        x_1101 = arg_1071 * y_1100;
                        
                        double y_1102;
                        
                        y_1102 = x_1101 + 687.1870074920579;
                        
                        double x_1103;
                        
                        x_1103 = arg_1071 * y_1102;
                        
                        double y_1104;
                        
                        y_1104 = x_1103 + 42.31333070160091;
                        
                        double x_1105;
                        
                        x_1105 = arg_1071 * y_1104;
                        
                        double y_1106;
                        
                        y_1106 = x_1105 + 1.0;
                        
                        double res_1596;
                        
                        if (cond_1066) {
                            double res_1107;
                            
                            res_1107 = x_1087 / y_1106;
                            
                            double res_1108;
                            
                            res_1108 = res_1058 * res_1107;
                            res_1596 = res_1108;
                        } else {
                            double negate_arg_1185;
                            
                            negate_arg_1185 = futhark_log64(res_1069);
                            
                            double arg_1186;
                            
                            arg_1186 = -negate_arg_1185;
                            
                            double res_1187;
                            
                            res_1187 = futhark_sqrt64(arg_1186);
                            
                            char cond_1188;
                            
                            cond_1188 = res_1187 <= 5.0;
                            
                            double arg_1189;
                            
                            arg_1189 = res_1187 - 1.6;
                            
                            double x_1216;
                            
                            x_1216 = arg_1189 * 7.745450142783414e-4;
                            
                            double y_1225;
                            
                            y_1225 = x_1216 + 2.2723844989269184e-2;
                            
                            double x_1229;
                            
                            x_1229 = arg_1189 * y_1225;
                            
                            double y_1230;
                            
                            y_1230 = x_1229 + 0.2417807251774506;
                            
                            double x_1231;
                            
                            x_1231 = arg_1189 * y_1230;
                            
                            double y_1245;
                            
                            y_1245 = x_1231 + 1.2704582524523684;
                            
                            double x_1371;
                            
                            x_1371 = arg_1189 * y_1245;
                            
                            double y_1376;
                            
                            y_1376 = x_1371 + 3.6478483247632045;
                            
                            double x_1380;
                            
                            x_1380 = arg_1189 * y_1376;
                            
                            double y_1381;
                            
                            y_1381 = x_1380 + 5.769497221460691;
                            
                            double x_1382;
                            
                            x_1382 = arg_1189 * y_1381;
                            
                            double y_1413;
                            
                            y_1413 = x_1382 + 4.630337846156546;
                            
                            double x_1540;
                            
                            x_1540 = arg_1189 * y_1413;
                            
                            double x_1545;
                            
                            x_1545 = x_1540 + 1.4234371107496835;
                            
                            double x_1546;
                            
                            x_1546 = arg_1189 * 1.0507500716444169e-9;
                            
                            double y_1547;
                            
                            y_1547 = x_1546 + 5.475938084995345e-4;
                            
                            double x_1548;
                            
                            x_1548 = arg_1189 * y_1547;
                            
                            double y_1549;
                            
                            y_1549 = x_1548 + 1.5198666563616457e-2;
                            
                            double x_1550;
                            
                            x_1550 = arg_1189 * y_1549;
                            
                            double y_1551;
                            
                            y_1551 = x_1550 + 0.14810397642748008;
                            
                            double x_1552;
                            
                            x_1552 = arg_1189 * y_1551;
                            
                            double y_1553;
                            
                            y_1553 = x_1552 + 0.6897673349851;
                            
                            double x_1554;
                            
                            x_1554 = arg_1189 * y_1553;
                            
                            double y_1555;
                            
                            y_1555 = x_1554 + 1.6763848301838038;
                            
                            double x_1556;
                            
                            x_1556 = arg_1189 * y_1555;
                            
                            double y_1557;
                            
                            y_1557 = x_1556 + 2.053191626637759;
                            
                            double x_1558;
                            
                            x_1558 = arg_1189 * y_1557;
                            
                            double y_1559;
                            
                            y_1559 = x_1558 + 1.0;
                            
                            double arg_1560;
                            
                            arg_1560 = res_1187 - 5.0;
                            
                            double x_1561;
                            
                            x_1561 = arg_1560 * 2.0103343992922881e-7;
                            
                            double y_1562;
                            
                            y_1562 = x_1561 + 2.7115555687434876e-5;
                            
                            double x_1563;
                            
                            x_1563 = arg_1560 * y_1562;
                            
                            double y_1564;
                            
                            y_1564 = x_1563 + 1.2426609473880784e-3;
                            
                            double x_1565;
                            
                            x_1565 = arg_1560 * y_1564;
                            
                            double y_1568;
                            
                            y_1568 = x_1565 + 2.6532189526576124e-2;
                            
                            double x_1569;
                            
                            x_1569 = arg_1560 * y_1568;
                            
                            double y_1570;
                            
                            y_1570 = x_1569 + 0.29656057182850487;
                            
                            double x_1571;
                            
                            x_1571 = arg_1560 * y_1570;
                            
                            double y_1572;
                            
                            y_1572 = x_1571 + 1.7848265399172913;
                            
                            double x_1573;
                            
                            x_1573 = arg_1560 * y_1572;
                            
                            double y_1574;
                            
                            y_1574 = x_1573 + 5.463784911164114;
                            
                            double x_1575;
                            
                            x_1575 = arg_1560 * y_1574;
                            
                            double x_1576;
                            
                            x_1576 = x_1575 + 6.657904643501103;
                            
                            double x_1577;
                            
                            x_1577 = arg_1560 * 2.0442631033899398e-5;
                            
                            double y_1578;
                            
                            y_1578 = x_1577 + 1.421511758316446e-7;
                            
                            double x_1579;
                            
                            x_1579 = arg_1560 * y_1578;
                            
                            double y_1580;
                            
                            y_1580 = x_1579 + 1.8463183175100548e-5;
                            
                            double x_1581;
                            
                            x_1581 = arg_1560 * y_1580;
                            
                            double y_1582;
                            
                            y_1582 = x_1581 + 7.868691311456133e-4;
                            
                            double x_1583;
                            
                            x_1583 = arg_1560 * y_1582;
                            
                            double y_1584;
                            
                            y_1584 = x_1583 + 1.4875361290850615e-2;
                            
                            double x_1585;
                            
                            x_1585 = arg_1560 * y_1584;
                            
                            double y_1586;
                            
                            y_1586 = x_1585 + 0.1369298809227358;
                            
                            double x_1587;
                            
                            x_1587 = arg_1560 * y_1586;
                            
                            double y_1588;
                            
                            y_1588 = x_1587 + 0.599832206555888;
                            
                            double x_1589;
                            
                            x_1589 = arg_1560 * y_1588;
                            
                            double y_1590;
                            
                            y_1590 = x_1589 + 1.0;
                            
                            double res_1593;
                            
                            if (cond_1188) {
                                double res_1591;
                                
                                res_1591 = x_1545 / y_1559;
                                res_1593 = res_1591;
                            } else {
                                double res_1592;
                                
                                res_1592 = x_1576 / y_1590;
                                res_1593 = res_1592;
                            }
                            
                            double res_1594;
                            
                            res_1594 = 0.0 - res_1593;
                            
                            double res_1595;
                            
                            if (x_1059) {
                                res_1595 = res_1594;
                            } else {
                                res_1595 = res_1593;
                            }
                            res_1596 = res_1595;
                        }
                        *(double*) &mem_567[i_57 * 1 * 8 + 0] = res_1596;
                    }
                    memmove(mem_573 + 0, mem_404 + 0, num_dates_120 *
                            sizeof(double));
                    
                    double y_1600;
                    
                    y_1600 = *(double*) &mem_567[0 * 1 * 8 + 0];
                    
                    double lw_val_1601;
                    
                    lw_val_1601 = x_178 * y_1600;
                    *(double*) &mem_573[idx_175 * 1 * 8 + 0] = lw_val_1601;
                    for (int i_414 = 0; i_414 < num_dates_120; i_414++) {
                        *(double*) &mem_563[8 * (i_61 * num_dates_120) + 8 *
                                            i_414] = *(double*) &mem_573[8 *
                                                                         i_414];
                    }
                    for (int i_1604 = 0; i_1604 < bound_185; i_1604++) {
                        int i_1605;
                        
                        i_1605 = i_1604 + 1;
                        
                        char y_1606;
                        
                        y_1606 = i_1605 < size_115;
                        
                        int bounds_check_1607;
                        
                        if (!y_1606) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_1606",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:160:22-160:24");
                            abort();
                        }
                        
                        int x_1608;
                        
                        x_1608 = *(int*) &bb_inds_mem_371[(1 * size_115 +
                                                           i_1605 * 1) * 4 + 0];
                        
                        int res_1609;
                        
                        res_1609 = x_1608 - 1;
                        
                        int x_1610;
                        
                        x_1610 = *(int*) &bb_inds_mem_371[(2 * size_115 +
                                                           i_1605 * 1) * 4 + 0];
                        
                        int res_1611;
                        
                        res_1611 = x_1610 - 1;
                        
                        int x_1612;
                        
                        x_1612 = *(int*) &bb_inds_mem_371[(0 * size_115 +
                                                           i_1605 * 1) * 4 + 0];
                        
                        int res_1613;
                        
                        res_1613 = x_1612 - 1;
                        
                        char x_1614;
                        
                        x_1614 = 0 <= res_1611;
                        
                        char y_1615;
                        
                        y_1615 = res_1611 < num_dates_120;
                        
                        char assert_arg_1616;
                        
                        assert_arg_1616 = x_1614 && y_1615;
                        
                        int bounds_check_1617;
                        
                        if (!assert_arg_1616) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "assert_arg_1616",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:164:22-164:27");
                            abort();
                        }
                        
                        double res_1618;
                        
                        res_1618 = *(double*) &mem_563[(i_61 * num_dates_120 +
                                                        res_1611 * 1) * 8 + 0];
                        
                        double res_1619;
                        
                        res_1619 = *(double*) &mem_567[i_1605 * 1 * 8 + 0];
                        
                        char y_1620;
                        
                        y_1620 = i_1605 < size_117;
                        
                        int bounds_check_1621;
                        
                        if (!y_1620) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_1620",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:166:22-166:24");
                            abort();
                        }
                        
                        double x_1622;
                        
                        x_1622 = *(double*) &bb_data_mem_381[(2 * size_117 +
                                                              i_1605 * 1) * 8 +
                                                             0];
                        
                        double x_1623;
                        
                        x_1623 = x_1622 * res_1618;
                        
                        double x_1624;
                        
                        x_1624 = *(double*) &bb_data_mem_381[(0 * size_117 +
                                                              i_1605 * 1) * 8 +
                                                             0];
                        
                        double y_1625;
                        
                        y_1625 = x_1624 * res_1619;
                        
                        double res_1626;
                        
                        res_1626 = x_1623 + y_1625;
                        
                        int x_1627;
                        
                        x_1627 = res_1609 + 1;
                        
                        char cond_1628;
                        
                        cond_1628 = x_1627 == 0;
                        
                        char x_1629;
                        
                        x_1629 = 0 <= res_1609;
                        
                        char y_1630;
                        
                        y_1630 = res_1609 < num_dates_120;
                        
                        char assert_arg_1631;
                        
                        assert_arg_1631 = x_1629 && y_1630;
                        
                        double lw_val_1637;
                        
                        if (cond_1628) {
                            lw_val_1637 = res_1626;
                        } else {
                            double x_1632;
                            
                            x_1632 = *(double*) &bb_data_mem_381[(1 * size_117 +
                                                                  i_1605 * 1) *
                                                                 8 + 0];
                            
                            int bounds_check_1633;
                            
                            if (!assert_arg_1631) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "assert_arg_1631",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:170:51-170:56");
                                abort();
                            }
                            
                            double y_1634;
                            
                            y_1634 = *(double*) &mem_563[(i_61 * num_dates_120 +
                                                          res_1609 * 1) * 8 +
                                                         0];
                            
                            double y_1635;
                            
                            y_1635 = x_1632 * y_1634;
                            
                            double res_1636;
                            
                            res_1636 = res_1626 + y_1635;
                            lw_val_1637 = res_1636;
                        }
                        
                        char x_1638;
                        
                        x_1638 = 0 <= res_1613;
                        
                        char y_1639;
                        
                        y_1639 = res_1613 < num_dates_120;
                        
                        char assert_arg_1640;
                        
                        assert_arg_1640 = x_1638 && y_1639;
                        
                        int bounds_check_1641;
                        
                        if (!assert_arg_1640) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "assert_arg_1640",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:168:13-168:16");
                            abort();
                        }
                        *(double*) &mem_563[(i_61 * num_dates_120 + res_1613 *
                                             1) * 8 + 0] = lw_val_1637;
                    }
                    for (int ii_1645 = 0; ii_1645 < bound_185; ii_1645++) {
                        int y_1646;
                        
                        y_1646 = ii_1645 + 1;
                        
                        int res_1647;
                        
                        res_1647 = num_dates_120 - y_1646;
                        
                        double x_1648;
                        
                        x_1648 = *(double*) &mem_563[(i_61 * num_dates_120 +
                                                      res_1647 * 1) * 8 + 0];
                        
                        int i_1649;
                        
                        i_1649 = res_1647 - 1;
                        
                        double y_1650;
                        
                        y_1650 = *(double*) &mem_563[(i_61 * num_dates_120 +
                                                      i_1649 * 1) * 8 + 0];
                        
                        double lw_val_1651;
                        
                        lw_val_1651 = x_1648 - y_1650;
                        *(double*) &mem_563[(i_61 * num_dates_120 + res_1647 *
                                             1) * 8 + 0] = lw_val_1651;
                    }
                }
                for (int i_73 = 0; i_73 < size_155; i_73++) {
                    int x_344;
                    
                    x_344 = i_73 * size_92;
                    
                    int x_353;
                    
                    x_353 = i_73 * size_89;
                    
                    int x_325;
                    
                    x_325 = i_73 * num_und_121;
                    for (int i_69 = 0; i_69 < size_92; i_69++) {
                        int flat_index_346;
                        
                        flat_index_346 = x_344 + i_69;
                        
                        int i_347;
                        
                        i_347 = flat_index_346 / size_92 - (flat_index_346 %
                                                            size_92 != 0 &&
                                                            flat_index_346 %
                                                            size_92 < 0 !=
                                                            size_92 <
                                                            0 ? 1 : 0);
                        
                        int y_348;
                        
                        y_348 = i_347 * size_92;
                        
                        int i_remnant_349;
                        
                        i_remnant_349 = flat_index_346 - y_348;
                        
                        double x_1688;
                        
                        x_1688 = *(double*) &all_md_drifts_mem_350[(0 *
                                                                    (size_91 *
                                                                     size_92) +
                                                                    i_347 *
                                                                    size_92 +
                                                                    i_remnant_349 *
                                                                    1) * 8 + 0];
                        
                        int flat_index_355;
                        
                        flat_index_355 = x_353 + i_69;
                        
                        int i_356;
                        
                        i_356 = flat_index_355 / size_89 - (flat_index_355 %
                                                            size_89 != 0 &&
                                                            flat_index_355 %
                                                            size_89 < 0 !=
                                                            size_89 <
                                                            0 ? 1 : 0);
                        
                        int y_357;
                        
                        y_357 = i_356 * size_89;
                        
                        int i_remnant_358;
                        
                        i_remnant_358 = flat_index_355 - y_357;
                        
                        double y_1690;
                        
                        y_1690 = *(double*) &all_md_vols_mem_268[(0 * (size_86 *
                                                                       size_89) +
                                                                  i_356 *
                                                                  size_89 +
                                                                  i_remnant_358 *
                                                                  1) * 8 + 0];
                        
                        int arg_1691;
                        
                        arg_1691 = i_69 + 1;
                        
                        char y_1695;
                        
                        y_1695 = i_69 < size_20;
                        
                        int bounds_check_1696;
                        
                        if (!y_1695) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_1695",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:205:58-205:62");
                            abort();
                        }
                        
                        char split_index_cnd_1698;
                        
                        split_index_cnd_1698 = arg_1691 <= size_21;
                        
                        int split_index_assert_1699;
                        
                        if (!split_index_cnd_1698) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "split_index_cnd_1698",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:201:56-201:61");
                            abort();
                        }
                        
                        double res_1711;
                        double acc_66;
                        
                        acc_66 = 0.0;
                        for (int i_65 = 0; i_65 < arg_1691; i_65++) {
                            double x_1707;
                            
                            x_1707 = *(double*) &all_md_c_mem_84[(0 * (size_20 *
                                                                       size_21) +
                                                                  i_69 *
                                                                  size_21 +
                                                                  i_65 * 1) *
                                                                 8 + 0];
                            
                            int flat_index_342;
                            
                            flat_index_342 = x_325 + i_65;
                            
                            int i_343;
                            
                            i_343 = flat_index_342 / num_und_121 -
                                (flat_index_342 % num_und_121 != 0 &&
                                 flat_index_342 % num_und_121 < 0 !=
                                 num_und_121 < 0 ? 1 : 0);
                            
                            int y_351;
                            
                            y_351 = i_343 * num_und_121;
                            
                            int i_remnant_352;
                            
                            i_remnant_352 = flat_index_342 - y_351;
                            
                            double x_1708;
                            
                            x_1708 = *(double*) &mem_563[(i_remnant_352 *
                                                          num_dates_120 +
                                                          i_343 * 1) * 8 + 0];
                            
                            double res_1709;
                            
                            res_1709 = x_1708 * x_1707;
                            
                            double res_1710;
                            
                            res_1710 = acc_66 + res_1709;
                            
                            double acc_tmp_429;
                            
                            acc_tmp_429 = res_1710;
                            acc_66 = acc_tmp_429;
                        }
                        res_1711 = acc_66;
                        
                        double res_1712;
                        
                        res_1712 = res_1711 * y_1690;
                        
                        double res_1713;
                        
                        res_1713 = x_1688 + res_1712;
                        *(double*) &mem_596[(i_73 * size_92 + i_69 * 1) * 8 +
                                            0] = res_1713;
                    }
                }
                for (int i_88 = 0; i_88 < size_92; i_88++) {
                    double x_1719;
                    
                    x_1719 = *(double*) &all_md_st_mem_359[(0 * size_94 + i_88 *
                                                            1) * 8 + 0];
                    for (int i_77 = 0; i_77 < size_155; i_77++) {
                        double not_curried_1720;
                        
                        not_curried_1720 = *(double*) &mem_596[(i_77 * size_92 +
                                                                i_88 * 1) * 8 +
                                                               0];
                        
                        double curried_fun_result_1721;
                        
                        curried_fun_result_1721 =
                            futhark_exp64(not_curried_1720);
                        *(double*) &mem_611[i_77 * 1 * 8 + 0] =
                            curried_fun_result_1721;
                    }
                    
                    double acc_82;
                    
                    acc_82 = x_1719;
                    for (int i_81 = 0; i_81 < size_155; i_81++) {
                        double y_1724;
                        
                        y_1724 = *(double*) &mem_611[i_81 * 1 * 8 + 0];
                        
                        double res_1725;
                        
                        res_1725 = acc_82 * y_1724;
                        *(double*) &mem_607[(i_88 * size_155 + i_81 * 1) * 8 +
                                            0] = res_1725;
                        
                        double row_87;
                        
                        row_87 = *(double*) &mem_607[(i_88 * size_155 + i_81 *
                                                      1) * 8 + 0];
                        
                        double acc_tmp_436;
                        
                        acc_tmp_436 = row_87;
                        acc_82 = acc_tmp_436;
                    }
                }
                for (int i_227 = 0; i_227 < size_155; i_227++) {
                    int x_362;
                    
                    x_362 = i_227 * size_190;
                    for (int i_99 = 0; i_99 < size_190; i_99++) {
                        int flat_index_364;
                        
                        flat_index_364 = x_362 + i_99;
                        
                        int i_365;
                        
                        i_365 = flat_index_364 / size_92 - (flat_index_364 %
                                                            size_92 != 0 &&
                                                            flat_index_364 %
                                                            size_92 < 0 !=
                                                            size_92 <
                                                            0 ? 1 : 0);
                        
                        int y_366;
                        
                        y_366 = i_365 * size_92;
                        
                        int i_remnant_367;
                        
                        i_remnant_367 = flat_index_364 - y_366;
                        
                        double x_1731;
                        
                        x_1731 = *(double*) &mem_607[(i_remnant_367 * size_155 +
                                                      i_365 * 1) * 8 + 0];
                        
                        double y_1732;
                        
                        y_1732 = *(double*) &mem_552[i_99 * 1 * 8 + 0];
                        
                        double res_1733;
                        
                        res_1733 = x_1731 * y_1732;
                        *(double*) &mem_627[i_99 * 1 * 8 + 0] = res_1733;
                    }
                    
                    double reduce_ne_1735;
                    
                    reduce_ne_1735 = *(double*) &mem_627[0 * 1 * 8 + 0];
                    
                    double res_1751;
                    double acc_108;
                    
                    acc_108 = reduce_ne_1735;
                    for (int i_107 = 0; i_107 < size_190; i_107++) {
                        double y_1745;
                        
                        y_1745 = *(double*) &mem_552[i_107 * 1 * 8 + 0];
                        
                        int flat_index_382;
                        
                        flat_index_382 = x_362 + i_107;
                        
                        int i_383;
                        
                        i_383 = flat_index_382 / size_92 - (flat_index_382 %
                                                            size_92 != 0 &&
                                                            flat_index_382 %
                                                            size_92 < 0 !=
                                                            size_92 <
                                                            0 ? 1 : 0);
                        
                        int y_386;
                        
                        y_386 = i_383 * size_92;
                        
                        int i_remnant_391;
                        
                        i_remnant_391 = flat_index_382 - y_386;
                        
                        double x_1746;
                        
                        x_1746 = *(double*) &mem_607[(i_remnant_391 * size_155 +
                                                      i_383 * 1) * 8 + 0];
                        
                        double res_1747;
                        
                        res_1747 = x_1746 * y_1745;
                        
                        char cond_1749;
                        
                        cond_1749 = acc_108 < res_1747;
                        
                        double res_1750;
                        
                        if (cond_1749) {
                            res_1750 = acc_108;
                        } else {
                            res_1750 = res_1747;
                        }
                        
                        double acc_tmp_484;
                        
                        acc_tmp_484 = res_1750;
                        acc_108 = acc_tmp_484;
                    }
                    res_1751 = acc_108;
                    *(double*) &mem_623[i_227 * 1 * 8 + 0] = res_1751;
                }
                
                double y_1758;
                
                y_1758 = *(double*) &mem_623[0 * 1 * 8 + 0];
                
                char cond_1759;
                
                cond_1759 = 1.0 <= y_1758;
                
                double res_1843;
                
                if (cond_1759) {
                    int bounds_check_1760;
                    
                    if (!y_213) {
                        fprintf(stderr, "Assertion %s at %s failed.\n", "y_213",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                        abort();
                    }
                    
                    double y_1761;
                    
                    y_1761 = *(double*) &all_md_disc_mem_368[(0 * size_113 + 0 *
                                                              1) * 8 + 0];
                    
                    double res_1762;
                    
                    res_1762 = 1150.0 * y_1761;
                    res_1843 = res_1762;
                } else {
                    int bounds_check_1763;
                    
                    if (!y_215) {
                        fprintf(stderr, "Assertion %s at %s failed.\n", "y_215",
                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:295:25-295:29");
                        abort();
                    }
                    
                    double y_1764;
                    
                    y_1764 = *(double*) &mem_623[1 * 1 * 8 + 0];
                    
                    char cond_1765;
                    
                    cond_1765 = 1.0 <= y_1764;
                    
                    double res_1842;
                    
                    if (cond_1765) {
                        int bounds_check_1766;
                        
                        if (!y_219) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_219",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                            abort();
                        }
                        
                        double y_1767;
                        
                        y_1767 = *(double*) &all_md_disc_mem_368[(0 * size_113 +
                                                                  1 * 1) * 8 +
                                                                 0];
                        
                        double res_1768;
                        
                        res_1768 = 1300.0 * y_1767;
                        res_1842 = res_1768;
                    } else {
                        int bounds_check_1769;
                        
                        if (!y_216) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_216",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:296:30-296:34");
                            abort();
                        }
                        
                        double y_1770;
                        
                        y_1770 = *(double*) &mem_623[2 * 1 * 8 + 0];
                        
                        char cond_1771;
                        
                        cond_1771 = 1.0 <= y_1770;
                        
                        double res_1840;
                        
                        if (cond_1771) {
                            int bounds_check_1772;
                            
                            if (!y_220) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_220",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                                abort();
                            }
                            
                            double y_1773;
                            
                            y_1773 = *(double*) &all_md_disc_mem_368[(0 *
                                                                      size_113 +
                                                                      2 * 1) *
                                                                     8 + 0];
                            
                            double res_1790;
                            
                            res_1790 = 1450.0 * y_1773;
                            res_1840 = res_1790;
                        } else {
                            int bounds_check_1791;
                            
                            if (!y_217) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_217",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:297:35-297:39");
                                abort();
                            }
                            
                            double y_1792;
                            
                            y_1792 = *(double*) &mem_623[3 * 1 * 8 + 0];
                            
                            char cond_1793;
                            
                            cond_1793 = 1.0 <= y_1792;
                            
                            double res_1827;
                            
                            if (cond_1793) {
                                int bounds_check_1803;
                                
                                if (!y_221) {
                                    fprintf(stderr,
                                            "Assertion %s at %s failed.\n",
                                            "y_221",
                                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                                    abort();
                                }
                                
                                double y_1805;
                                
                                y_1805 = *(double*) &all_md_disc_mem_368[(0 *
                                                                          size_113 +
                                                                          3 *
                                                                          1) *
                                                                         8 + 0];
                                
                                double res_1806;
                                
                                res_1806 = 1600.0 * y_1805;
                                res_1827 = res_1806;
                            } else {
                                int bounds_check_1807;
                                
                                if (!y_218) {
                                    fprintf(stderr,
                                            "Assertion %s at %s failed.\n",
                                            "y_218",
                                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:298:40-298:44");
                                    abort();
                                }
                                
                                double y_1808;
                                
                                y_1808 = *(double*) &mem_623[4 * 1 * 8 + 0];
                                
                                char cond_1809;
                                
                                cond_1809 = 1.0 <= y_1808;
                                
                                char cond_1810;
                                
                                cond_1810 = 0.75 < y_1808;
                                
                                double arg_1812;
                                
                                arg_1812 = 1000.0 * y_1808;
                                
                                double res_1826;
                                
                                if (cond_1809) {
                                    int bounds_check_1814;
                                    
                                    if (!y_222) {
                                        fprintf(stderr,
                                                "Assertion %s at %s failed.\n",
                                                "y_222",
                                                "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                                        abort();
                                    }
                                    
                                    double y_1815;
                                    
                                    y_1815 =
                                        *(double*) &all_md_disc_mem_368[(0 *
                                                                         size_113 +
                                                                         4 *
                                                                         1) *
                                                                        8 + 0];
                                    
                                    double res_1816;
                                    
                                    res_1816 = 1750.0 * y_1815;
                                    res_1826 = res_1816;
                                } else {
                                    double res_1825;
                                    
                                    if (cond_1810) {
                                        int bounds_check_1817;
                                        
                                        if (!y_222) {
                                            fprintf(stderr,
                                                    "Assertion %s at %s failed.\n",
                                                    "y_222",
                                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                                            abort();
                                        }
                                        
                                        double y_1818;
                                        
                                        y_1818 =
                                            *(double*) &all_md_disc_mem_368[(0 *
                                                                             size_113 +
                                                                             4 *
                                                                             1) *
                                                                            8 +
                                                                            0];
                                        
                                        double res_1820;
                                        
                                        res_1820 = 1000.0 * y_1818;
                                        res_1825 = res_1820;
                                    } else {
                                        int bounds_check_1822;
                                        
                                        if (!y_222) {
                                            fprintf(stderr,
                                                    "Assertion %s at %s failed.\n",
                                                    "y_222",
                                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                                            abort();
                                        }
                                        
                                        double y_1823;
                                        
                                        y_1823 =
                                            *(double*) &all_md_disc_mem_368[(0 *
                                                                             size_113 +
                                                                             4 *
                                                                             1) *
                                                                            8 +
                                                                            0];
                                        
                                        double res_1824;
                                        
                                        res_1824 = arg_1812 * y_1823;
                                        res_1825 = res_1824;
                                    }
                                    res_1826 = res_1825;
                                }
                                res_1827 = res_1826;
                            }
                            res_1840 = res_1827;
                        }
                        res_1842 = res_1840;
                    }
                    res_1843 = res_1842;
                }
                
                double res_1844;
                
                res_1844 = acc_232 + res_1843;
                
                double acc_tmp_323;
                
                acc_tmp_323 = res_1844;
                acc_232 = acc_tmp_323;
            }
            res_1845 = acc_232;
            res_2306 = res_1845;
        } else {
            double res_2305;
            
            if (cond_209) {
                int bounds_check_1846;
                
                if (!y_212) {
                    fprintf(stderr, "Assertion %s at %s failed.\n", "y_212",
                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:306:29-306:31");
                    abort();
                }
                
                int bounds_check_1847;
                
                if (!y_223) {
                    fprintf(stderr, "Assertion %s at %s failed.\n", "y_223",
                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:307:29-307:31");
                    abort();
                }
                
                int bounds_check_1848;
                
                if (!y_224) {
                    fprintf(stderr, "Assertion %s at %s failed.\n", "y_224",
                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:308:29-308:31");
                    abort();
                }
                
                int bounds_check_1859;
                
                if (!y_281) {
                    fprintf(stderr, "Assertion %s at %s failed.\n", "y_281",
                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:311:19-311:22");
                    abort();
                }
                
                int bounds_check_1860;
                
                if (!y_213) {
                    fprintf(stderr, "Assertion %s at %s failed.\n", "y_213",
                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                    abort();
                }
                
                double y_1861;
                
                y_1861 = *(double*) &all_md_disc_mem_368[(0 * size_113 + 0 *
                                                          1) * 8 + 0];
                
                double res_1862;
                
                res_1862 = 100.0 * y_1861;
                
                double res_2304;
                double acc_279;
                
                acc_279 = 0.0;
                for (int i_278 = 0; i_278 < num_mc_it_119; i_278++) {
                    int res_1873;
                    
                    res_1873 = i_278 + 1;
                    
                    int x_1917;
                    
                    x_1917 = res_1873 >> 1;
                    
                    int res_1925;
                    
                    res_1925 = x_1917 ^ res_1873;
                    for (int i_234 = 0; i_234 < num_bits_123; i_234++) {
                        int res_1927;
                        
                        res_1927 = 1 << i_234;
                        
                        int x_1928;
                        
                        x_1928 = res_1925 & res_1927;
                        
                        char res_1929;
                        
                        res_1929 = x_1928 == res_1927;
                        
                        int partition_equivalence_class_1930;
                        
                        if (res_1929) {
                            partition_equivalence_class_1930 = 0;
                        } else {
                            partition_equivalence_class_1930 = 1;
                        }
                        *(int*) &mem_640[i_234 * 1 * 4 + 0] =
                            partition_equivalence_class_1930;
                    }
                    
                    int filter_size_1932;
                    
                    filter_size_1932 = 0;
                    
                    int eqclass_492;
                    
                    for (int i_491 = 0; i_491 < num_bits_123; i_491++) {
                        eqclass_492 = *(int*) &mem_640[i_491 * 1 * 4 + 0];
                        if (eqclass_492 == 0) {
                            filter_size_1932 = filter_size_1932 + 1;
                        } else { }
                    }
                    
                    int partition_offset_493;
                    
                    partition_offset_493 = 0;
                    
                    int partition_cur_offset_494;
                    
                    for (int i_491 = 0; i_491 < num_bits_123; i_491++) {
                        eqclass_492 = *(int*) &mem_640[i_491 * 1 * 4 + 0];
                        if (eqclass_492 == 0) {
                            partition_cur_offset_494 = partition_offset_493;
                            *(int*) &mem_649[partition_cur_offset_494 * 1 * 4 +
                                             0] = *(int*) &mem_387[i_491 * 1 *
                                                                   4 + 0];
                            partition_offset_493 = partition_offset_493 + 1;
                        } else { }
                    }
                    for (int i_245 = 0; i_245 < num_und_121; i_245++) {
                        int y_374;
                        
                        y_374 = i_245 * size_18;
                        for (int i_241 = 0; i_241 < num_dates_120; i_241++) {
                            int x_373;
                            
                            x_373 = i_241 * new_slice_size_286;
                            
                            int x_375;
                            
                            x_375 = x_373 + y_374;
                            
                            int res_1978;
                            int acc_239;
                            
                            acc_239 = 0;
                            for (int i_238 = 0; i_238 < filter_size_1932;
                                 i_238++) {
                                int not_curried_1973;
                                
                                not_curried_1973 = *(int*) &mem_649[i_238 * 1 *
                                                                    4 + 0];
                                
                                char y_1974;
                                
                                y_1974 = not_curried_1973 < size_18;
                                
                                int bounds_check_1975;
                                
                                if (!y_1974) {
                                    fprintf(stderr,
                                            "Assertion %s at %s failed.\n",
                                            "y_1974",
                                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:6:35-6:38");
                                    abort();
                                }
                                
                                int flat_index_377;
                                
                                flat_index_377 = x_375 + not_curried_1973;
                                
                                int i_378;
                                
                                i_378 = flat_index_377 / size_18 -
                                    (flat_index_377 % size_18 != 0 &&
                                     flat_index_377 % size_18 < 0 != size_18 <
                                     0 ? 1 : 0);
                                
                                int y_379;
                                
                                y_379 = i_378 * size_18;
                                
                                int i_remnant_380;
                                
                                i_remnant_380 = flat_index_377 - y_379;
                                
                                int res_1976;
                                
                                res_1976 = *(int*) &dir_vs_mem_40[(i_378 *
                                                                   size_18 +
                                                                   i_remnant_380 *
                                                                   1) * 4 + 0];
                                
                                int res_1977;
                                
                                res_1977 = acc_239 ^ res_1976;
                                
                                int acc_tmp_520;
                                
                                acc_tmp_520 = res_1977;
                                acc_239 = acc_tmp_520;
                            }
                            res_1978 = acc_239;
                            
                            double x_1979;
                            
                            x_1979 = futhark_toFloat64(res_1978);
                            
                            double res_1980;
                            
                            res_1980 = x_1979 / res_151;
                            
                            double res_1981;
                            
                            res_1981 = res_1980 - 0.5;
                            
                            char x_1982;
                            
                            x_1982 = res_1981 < 0.0;
                            
                            double x_2000;
                            
                            x_2000 = 0.0 - res_1981;
                            
                            char y_2001;
                            
                            y_2001 = x_2000 <= 0.425;
                            
                            char x_2002;
                            
                            x_2002 = x_1982 && y_2001;
                            
                            char x_2003;
                            
                            x_2003 = 0.0 <= res_1981;
                            
                            char y_2004;
                            
                            y_2004 = res_1981 <= 0.425;
                            
                            char y_2005;
                            
                            y_2005 = x_2003 && y_2004;
                            
                            char cond_2006;
                            
                            cond_2006 = x_2002 || y_2005;
                            
                            double res_2007;
                            
                            res_2007 = res_1981 + 0.5;
                            
                            double res_2008;
                            
                            res_2008 = 0.5 - res_1981;
                            
                            double res_2009;
                            
                            if (x_1982) {
                                res_2009 = res_2007;
                            } else {
                                res_2009 = res_2008;
                            }
                            
                            double y_2010;
                            
                            y_2010 = res_1981 * res_1981;
                            
                            double arg_2011;
                            
                            arg_2011 = 0.180625 - y_2010;
                            
                            double x_2012;
                            
                            x_2012 = arg_2011 * 2509.0809287301227;
                            
                            double y_2013;
                            
                            y_2013 = x_2012 + 33430.57558358813;
                            
                            double x_2014;
                            
                            x_2014 = arg_2011 * y_2013;
                            
                            double y_2015;
                            
                            y_2015 = x_2014 + 67265.7709270087;
                            
                            double x_2016;
                            
                            x_2016 = arg_2011 * y_2015;
                            
                            double y_2017;
                            
                            y_2017 = x_2016 + 45921.95393154987;
                            
                            double x_2019;
                            
                            x_2019 = arg_2011 * y_2017;
                            
                            double y_2020;
                            
                            y_2020 = x_2019 + 13731.69376550946;
                            
                            double x_2021;
                            
                            x_2021 = arg_2011 * y_2020;
                            
                            double y_2022;
                            
                            y_2022 = x_2021 + 1971.5909503065513;
                            
                            double x_2023;
                            
                            x_2023 = arg_2011 * y_2022;
                            
                            double y_2034;
                            
                            y_2034 = x_2023 + 133.14166789178438;
                            
                            double x_2038;
                            
                            x_2038 = arg_2011 * y_2034;
                            
                            double x_2039;
                            
                            x_2039 = x_2038 + 3.3871328727963665;
                            
                            double x_2040;
                            
                            x_2040 = arg_2011 * 5226.495278852854;
                            
                            double y_2041;
                            
                            y_2041 = x_2040 + 28729.085735721943;
                            
                            double x_2042;
                            
                            x_2042 = arg_2011 * y_2041;
                            
                            double y_2043;
                            
                            y_2043 = x_2042 + 39307.89580009271;
                            
                            double x_2044;
                            
                            x_2044 = arg_2011 * y_2043;
                            
                            double y_2045;
                            
                            y_2045 = x_2044 + 21213.794301586597;
                            
                            double x_2046;
                            
                            x_2046 = arg_2011 * y_2045;
                            
                            double y_2047;
                            
                            y_2047 = x_2046 + 5394.196021424751;
                            
                            double x_2048;
                            
                            x_2048 = arg_2011 * y_2047;
                            
                            double y_2049;
                            
                            y_2049 = x_2048 + 687.1870074920579;
                            
                            double x_2050;
                            
                            x_2050 = arg_2011 * y_2049;
                            
                            double y_2051;
                            
                            y_2051 = x_2050 + 42.31333070160091;
                            
                            double x_2052;
                            
                            x_2052 = arg_2011 * y_2051;
                            
                            double y_2053;
                            
                            y_2053 = x_2052 + 1.0;
                            
                            double res_2123;
                            
                            if (cond_2006) {
                                double res_2054;
                                
                                res_2054 = x_2039 / y_2053;
                                
                                double res_2055;
                                
                                res_2055 = res_1981 * res_2054;
                                res_2123 = res_2055;
                            } else {
                                double negate_arg_2056;
                                
                                negate_arg_2056 = futhark_log64(res_2009);
                                
                                double arg_2057;
                                
                                arg_2057 = -negate_arg_2056;
                                
                                double res_2058;
                                
                                res_2058 = futhark_sqrt64(arg_2057);
                                
                                char cond_2059;
                                
                                cond_2059 = res_2058 <= 5.0;
                                
                                double arg_2060;
                                
                                arg_2060 = res_2058 - 1.6;
                                
                                double x_2061;
                                
                                x_2061 = arg_2060 * 7.745450142783414e-4;
                                
                                double y_2062;
                                
                                y_2062 = x_2061 + 2.2723844989269184e-2;
                                
                                double x_2063;
                                
                                x_2063 = arg_2060 * y_2062;
                                
                                double y_2064;
                                
                                y_2064 = x_2063 + 0.2417807251774506;
                                
                                double x_2065;
                                
                                x_2065 = arg_2060 * y_2064;
                                
                                double y_2066;
                                
                                y_2066 = x_2065 + 1.2704582524523684;
                                
                                double x_2067;
                                
                                x_2067 = arg_2060 * y_2066;
                                
                                double y_2068;
                                
                                y_2068 = x_2067 + 3.6478483247632045;
                                
                                double x_2069;
                                
                                x_2069 = arg_2060 * y_2068;
                                
                                double y_2070;
                                
                                y_2070 = x_2069 + 5.769497221460691;
                                
                                double x_2071;
                                
                                x_2071 = arg_2060 * y_2070;
                                
                                double y_2072;
                                
                                y_2072 = x_2071 + 4.630337846156546;
                                
                                double x_2073;
                                
                                x_2073 = arg_2060 * y_2072;
                                
                                double x_2074;
                                
                                x_2074 = x_2073 + 1.4234371107496835;
                                
                                double x_2075;
                                
                                x_2075 = arg_2060 * 1.0507500716444169e-9;
                                
                                double y_2076;
                                
                                y_2076 = x_2075 + 5.475938084995345e-4;
                                
                                double x_2077;
                                
                                x_2077 = arg_2060 * y_2076;
                                
                                double y_2078;
                                
                                y_2078 = x_2077 + 1.5198666563616457e-2;
                                
                                double x_2079;
                                
                                x_2079 = arg_2060 * y_2078;
                                
                                double y_2080;
                                
                                y_2080 = x_2079 + 0.14810397642748008;
                                
                                double x_2081;
                                
                                x_2081 = arg_2060 * y_2080;
                                
                                double y_2082;
                                
                                y_2082 = x_2081 + 0.6897673349851;
                                
                                double x_2083;
                                
                                x_2083 = arg_2060 * y_2082;
                                
                                double y_2084;
                                
                                y_2084 = x_2083 + 1.6763848301838038;
                                
                                double x_2085;
                                
                                x_2085 = arg_2060 * y_2084;
                                
                                double y_2086;
                                
                                y_2086 = x_2085 + 2.053191626637759;
                                
                                double x_2087;
                                
                                x_2087 = arg_2060 * y_2086;
                                
                                double y_2088;
                                
                                y_2088 = x_2087 + 1.0;
                                
                                double arg_2089;
                                
                                arg_2089 = res_2058 - 5.0;
                                
                                double x_2090;
                                
                                x_2090 = arg_2089 * 2.0103343992922881e-7;
                                
                                double y_2091;
                                
                                y_2091 = x_2090 + 2.7115555687434876e-5;
                                
                                double x_2092;
                                
                                x_2092 = arg_2089 * y_2091;
                                
                                double y_2093;
                                
                                y_2093 = x_2092 + 1.2426609473880784e-3;
                                
                                double x_2094;
                                
                                x_2094 = arg_2089 * y_2093;
                                
                                double y_2095;
                                
                                y_2095 = x_2094 + 2.6532189526576124e-2;
                                
                                double x_2096;
                                
                                x_2096 = arg_2089 * y_2095;
                                
                                double y_2097;
                                
                                y_2097 = x_2096 + 0.29656057182850487;
                                
                                double x_2098;
                                
                                x_2098 = arg_2089 * y_2097;
                                
                                double y_2099;
                                
                                y_2099 = x_2098 + 1.7848265399172913;
                                
                                double x_2100;
                                
                                x_2100 = arg_2089 * y_2099;
                                
                                double y_2101;
                                
                                y_2101 = x_2100 + 5.463784911164114;
                                
                                double x_2102;
                                
                                x_2102 = arg_2089 * y_2101;
                                
                                double x_2103;
                                
                                x_2103 = x_2102 + 6.657904643501103;
                                
                                double x_2104;
                                
                                x_2104 = arg_2089 * 2.0442631033899398e-5;
                                
                                double y_2105;
                                
                                y_2105 = x_2104 + 1.421511758316446e-7;
                                
                                double x_2106;
                                
                                x_2106 = arg_2089 * y_2105;
                                
                                double y_2107;
                                
                                y_2107 = x_2106 + 1.8463183175100548e-5;
                                
                                double x_2108;
                                
                                x_2108 = arg_2089 * y_2107;
                                
                                double y_2109;
                                
                                y_2109 = x_2108 + 7.868691311456133e-4;
                                
                                double x_2110;
                                
                                x_2110 = arg_2089 * y_2109;
                                
                                double y_2111;
                                
                                y_2111 = x_2110 + 1.4875361290850615e-2;
                                
                                double x_2112;
                                
                                x_2112 = arg_2089 * y_2111;
                                
                                double y_2113;
                                
                                y_2113 = x_2112 + 0.1369298809227358;
                                
                                double x_2114;
                                
                                x_2114 = arg_2089 * y_2113;
                                
                                double y_2115;
                                
                                y_2115 = x_2114 + 0.599832206555888;
                                
                                double x_2116;
                                
                                x_2116 = arg_2089 * y_2115;
                                
                                double y_2117;
                                
                                y_2117 = x_2116 + 1.0;
                                
                                double res_2120;
                                
                                if (cond_2059) {
                                    double res_2118;
                                    
                                    res_2118 = x_2074 / y_2088;
                                    res_2120 = res_2118;
                                } else {
                                    double res_2119;
                                    
                                    res_2119 = x_2103 / y_2117;
                                    res_2120 = res_2119;
                                }
                                
                                double res_2121;
                                
                                res_2121 = 0.0 - res_2120;
                                
                                double res_2122;
                                
                                if (x_1982) {
                                    res_2122 = res_2121;
                                } else {
                                    res_2122 = res_2120;
                                }
                                res_2123 = res_2122;
                            }
                            *(double*) &mem_656[i_241 * 1 * 8 + 0] = res_2123;
                        }
                        memmove(mem_662 + 0, mem_404 + 0, num_dates_120 *
                                sizeof(double));
                        
                        double y_2127;
                        
                        y_2127 = *(double*) &mem_656[0 * 1 * 8 + 0];
                        
                        double lw_val_2128;
                        
                        lw_val_2128 = x_178 * y_2127;
                        *(double*) &mem_662[idx_175 * 1 * 8 + 0] = lw_val_2128;
                        for (int i_521 = 0; i_521 < num_dates_120; i_521++) {
                            *(double*) &mem_652[8 * (num_dates_120 * i_245) +
                                                8 * i_521] =
                                *(double*) &mem_662[8 * i_521];
                        }
                        for (int i_2131 = 0; i_2131 < bound_185; i_2131++) {
                            int i_2132;
                            
                            i_2132 = i_2131 + 1;
                            
                            char y_2133;
                            
                            y_2133 = i_2132 < size_115;
                            
                            int bounds_check_2134;
                            
                            if (!y_2133) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_2133",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:160:22-160:24");
                                abort();
                            }
                            
                            int x_2135;
                            
                            x_2135 = *(int*) &bb_inds_mem_371[(1 * size_115 +
                                                               i_2132 * 1) * 4 +
                                                              0];
                            
                            int res_2136;
                            
                            res_2136 = x_2135 - 1;
                            
                            int x_2137;
                            
                            x_2137 = *(int*) &bb_inds_mem_371[(2 * size_115 +
                                                               i_2132 * 1) * 4 +
                                                              0];
                            
                            int res_2138;
                            
                            res_2138 = x_2137 - 1;
                            
                            int x_2139;
                            
                            x_2139 = *(int*) &bb_inds_mem_371[(0 * size_115 +
                                                               i_2132 * 1) * 4 +
                                                              0];
                            
                            int res_2140;
                            
                            res_2140 = x_2139 - 1;
                            
                            char x_2141;
                            
                            x_2141 = 0 <= res_2138;
                            
                            char y_2142;
                            
                            y_2142 = res_2138 < num_dates_120;
                            
                            char assert_arg_2143;
                            
                            assert_arg_2143 = x_2141 && y_2142;
                            
                            int bounds_check_2144;
                            
                            if (!assert_arg_2143) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "assert_arg_2143",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:164:22-164:27");
                                abort();
                            }
                            
                            double res_2145;
                            
                            res_2145 = *(double*) &mem_652[(i_245 *
                                                            num_dates_120 +
                                                            res_2138 * 1) * 8 +
                                                           0];
                            
                            double res_2146;
                            
                            res_2146 = *(double*) &mem_656[i_2132 * 1 * 8 + 0];
                            
                            char y_2147;
                            
                            y_2147 = i_2132 < size_117;
                            
                            int bounds_check_2148;
                            
                            if (!y_2147) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_2147",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:166:22-166:24");
                                abort();
                            }
                            
                            double x_2149;
                            
                            x_2149 = *(double*) &bb_data_mem_381[(2 * size_117 +
                                                                  i_2132 * 1) *
                                                                 8 + 0];
                            
                            double x_2150;
                            
                            x_2150 = x_2149 * res_2145;
                            
                            double x_2151;
                            
                            x_2151 = *(double*) &bb_data_mem_381[(0 * size_117 +
                                                                  i_2132 * 1) *
                                                                 8 + 0];
                            
                            double y_2152;
                            
                            y_2152 = x_2151 * res_2146;
                            
                            double res_2153;
                            
                            res_2153 = x_2150 + y_2152;
                            
                            int x_2154;
                            
                            x_2154 = res_2136 + 1;
                            
                            char cond_2155;
                            
                            cond_2155 = x_2154 == 0;
                            
                            char x_2156;
                            
                            x_2156 = 0 <= res_2136;
                            
                            char y_2157;
                            
                            y_2157 = res_2136 < num_dates_120;
                            
                            char assert_arg_2158;
                            
                            assert_arg_2158 = x_2156 && y_2157;
                            
                            double lw_val_2164;
                            
                            if (cond_2155) {
                                lw_val_2164 = res_2153;
                            } else {
                                double x_2159;
                                
                                x_2159 = *(double*) &bb_data_mem_381[(1 *
                                                                      size_117 +
                                                                      i_2132 *
                                                                      1) * 8 +
                                                                     0];
                                
                                int bounds_check_2160;
                                
                                if (!assert_arg_2158) {
                                    fprintf(stderr,
                                            "Assertion %s at %s failed.\n",
                                            "assert_arg_2158",
                                            "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:170:51-170:56");
                                    abort();
                                }
                                
                                double y_2161;
                                
                                y_2161 = *(double*) &mem_652[(i_245 *
                                                              num_dates_120 +
                                                              res_2136 * 1) *
                                                             8 + 0];
                                
                                double y_2162;
                                
                                y_2162 = x_2159 * y_2161;
                                
                                double res_2163;
                                
                                res_2163 = res_2153 + y_2162;
                                lw_val_2164 = res_2163;
                            }
                            
                            char x_2165;
                            
                            x_2165 = 0 <= res_2140;
                            
                            char y_2166;
                            
                            y_2166 = res_2140 < num_dates_120;
                            
                            char assert_arg_2167;
                            
                            assert_arg_2167 = x_2165 && y_2166;
                            
                            int bounds_check_2168;
                            
                            if (!assert_arg_2167) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "assert_arg_2167",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:168:13-168:16");
                                abort();
                            }
                            *(double*) &mem_652[(i_245 * num_dates_120 +
                                                 res_2140 * 1) * 8 + 0] =
                                lw_val_2164;
                        }
                        for (int ii_2172 = 0; ii_2172 < bound_185; ii_2172++) {
                            int y_2173;
                            
                            y_2173 = ii_2172 + 1;
                            
                            int res_2174;
                            
                            res_2174 = num_dates_120 - y_2173;
                            
                            double x_2175;
                            
                            x_2175 = *(double*) &mem_652[(i_245 *
                                                          num_dates_120 +
                                                          res_2174 * 1) * 8 +
                                                         0];
                            
                            int i_2176;
                            
                            i_2176 = res_2174 - 1;
                            
                            double y_2177;
                            
                            y_2177 = *(double*) &mem_652[(i_245 *
                                                          num_dates_120 +
                                                          i_2176 * 1) * 8 + 0];
                            
                            double lw_val_2178;
                            
                            lw_val_2178 = x_2175 - y_2177;
                            *(double*) &mem_652[(i_245 * num_dates_120 +
                                                 res_2174 * 1) * 8 + 0] =
                                lw_val_2178;
                        }
                    }
                    for (int i_257 = 0; i_257 < size_155; i_257++) {
                        int x_384;
                        
                        x_384 = i_257 * size_92;
                        
                        int x_407;
                        
                        x_407 = i_257 * size_89;
                        
                        int x_417;
                        
                        x_417 = i_257 * num_und_121;
                        for (int i_253 = 0; i_253 < size_92; i_253++) {
                            int flat_index_400;
                            
                            flat_index_400 = x_384 + i_253;
                            
                            int i_401;
                            
                            i_401 = flat_index_400 / size_92 - (flat_index_400 %
                                                                size_92 != 0 &&
                                                                flat_index_400 %
                                                                size_92 < 0 !=
                                                                size_92 <
                                                                0 ? 1 : 0);
                            
                            int y_402;
                            
                            y_402 = i_401 * size_92;
                            
                            int i_remnant_403;
                            
                            i_remnant_403 = flat_index_400 - y_402;
                            
                            double x_2217;
                            
                            x_2217 = *(double*) &all_md_drifts_mem_350[(0 *
                                                                        (size_91 *
                                                                         size_92) +
                                                                        i_401 *
                                                                        size_92 +
                                                                        i_remnant_403 *
                                                                        1) * 8 +
                                                                       0];
                            
                            int flat_index_409;
                            
                            flat_index_409 = x_407 + i_253;
                            
                            int i_410;
                            
                            i_410 = flat_index_409 / size_89 - (flat_index_409 %
                                                                size_89 != 0 &&
                                                                flat_index_409 %
                                                                size_89 < 0 !=
                                                                size_89 <
                                                                0 ? 1 : 0);
                            
                            int y_411;
                            
                            y_411 = i_410 * size_89;
                            
                            int i_remnant_412;
                            
                            i_remnant_412 = flat_index_409 - y_411;
                            
                            double y_2219;
                            
                            y_2219 = *(double*) &all_md_vols_mem_268[(0 *
                                                                      (size_86 *
                                                                       size_89) +
                                                                      i_410 *
                                                                      size_89 +
                                                                      i_remnant_412 *
                                                                      1) * 8 +
                                                                     0];
                            
                            int arg_2220;
                            
                            arg_2220 = i_253 + 1;
                            
                            char y_2224;
                            
                            y_2224 = i_253 < size_20;
                            
                            int bounds_check_2225;
                            
                            if (!y_2224) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "y_2224",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:205:58-205:62");
                                abort();
                            }
                            
                            char split_index_cnd_2227;
                            
                            split_index_cnd_2227 = arg_2220 <= size_21;
                            
                            int split_index_assert_2228;
                            
                            if (!split_index_cnd_2227) {
                                fprintf(stderr, "Assertion %s at %s failed.\n",
                                        "split_index_cnd_2227",
                                        "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:201:56-201:61");
                                abort();
                            }
                            
                            double res_2240;
                            double acc_250;
                            
                            acc_250 = 0.0;
                            for (int i_249 = 0; i_249 < arg_2220; i_249++) {
                                int flat_index_419;
                                
                                flat_index_419 = x_417 + i_249;
                                
                                int i_420;
                                
                                i_420 = flat_index_419 / num_und_121 -
                                    (flat_index_419 % num_und_121 != 0 &&
                                     flat_index_419 % num_und_121 < 0 !=
                                     num_und_121 < 0 ? 1 : 0);
                                
                                int y_421;
                                
                                y_421 = i_420 * num_und_121;
                                
                                int i_remnant_422;
                                
                                i_remnant_422 = flat_index_419 - y_421;
                                
                                double x_2236;
                                
                                x_2236 = *(double*) &mem_652[(i_remnant_422 *
                                                              num_dates_120 +
                                                              i_420 * 1) * 8 +
                                                             0];
                                
                                double x_2237;
                                
                                x_2237 = *(double*) &all_md_c_mem_84[(0 *
                                                                      (size_20 *
                                                                       size_21) +
                                                                      i_253 *
                                                                      size_21 +
                                                                      i_249 *
                                                                      1) * 8 +
                                                                     0];
                                
                                double res_2238;
                                
                                res_2238 = x_2236 * x_2237;
                                
                                double res_2239;
                                
                                res_2239 = acc_250 + res_2238;
                                
                                double acc_tmp_529;
                                
                                acc_tmp_529 = res_2239;
                                acc_250 = acc_tmp_529;
                            }
                            res_2240 = acc_250;
                            
                            double res_2241;
                            
                            res_2241 = res_2240 * y_2219;
                            
                            double res_2242;
                            
                            res_2242 = x_2217 + res_2241;
                            *(double*) &mem_722[(i_257 * size_92 + i_253 * 1) *
                                                8 + 0] = res_2242;
                        }
                    }
                    for (int i_271 = 0; i_271 < size_92; i_271++) {
                        double x_2248;
                        
                        x_2248 = *(double*) &all_md_st_mem_359[(0 * size_94 +
                                                                i_271 * 1) * 8 +
                                                               0];
                        for (int i_261 = 0; i_261 < size_155; i_261++) {
                            double not_curried_2249;
                            
                            not_curried_2249 = *(double*) &mem_722[(i_261 *
                                                                    size_92 +
                                                                    i_271 * 1) *
                                                                   8 + 0];
                            
                            double curried_fun_result_2250;
                            
                            curried_fun_result_2250 =
                                futhark_exp64(not_curried_2249);
                            *(double*) &mem_777[i_261 * 1 * 8 + 0] =
                                curried_fun_result_2250;
                        }
                        
                        double acc_266;
                        
                        acc_266 = x_2248;
                        for (int i_265 = 0; i_265 < size_155; i_265++) {
                            double y_2253;
                            
                            y_2253 = *(double*) &mem_777[i_265 * 1 * 8 + 0];
                            
                            double res_2254;
                            
                            res_2254 = acc_266 * y_2253;
                            *(double*) &mem_771[(i_271 * size_155 + i_265 * 1) *
                                                8 + 0] = res_2254;
                            
                            double row_270;
                            
                            row_270 = *(double*) &mem_771[(i_271 * size_155 +
                                                           i_265 * 1) * 8 + 0];
                            
                            double acc_tmp_532;
                            
                            acc_tmp_532 = row_270;
                            acc_266 = acc_tmp_532;
                        }
                    }
                    
                    char res_2273;
                    char acc_276;
                    
                    acc_276 = 0;
                    for (int i_275 = 0; i_275 < size_155; i_275++) {
                        int x_427;
                        
                        x_427 = i_275 * size_190;
                        
                        int i_430;
                        
                        i_430 = x_427 / size_92 - (x_427 % size_92 != 0 &&
                                                   x_427 % size_92 < 0 !=
                                                   size_92 < 0 ? 1 : 0);
                        
                        int y_431;
                        
                        y_431 = i_430 * size_92;
                        
                        int i_remnant_432;
                        
                        i_remnant_432 = x_427 - y_431;
                        
                        double x_2264;
                        
                        x_2264 = *(double*) &mem_771[(i_remnant_432 * size_155 +
                                                      i_430 * 1) * 8 + 0];
                        
                        char x_2265;
                        
                        x_2265 = x_2264 <= 2630.6349999999998;
                        
                        int flat_index_477;
                        
                        flat_index_477 = x_427 + 1;
                        
                        int i_478;
                        
                        i_478 = flat_index_477 / size_92 - (flat_index_477 %
                                                            size_92 != 0 &&
                                                            flat_index_477 %
                                                            size_92 < 0 !=
                                                            size_92 <
                                                            0 ? 1 : 0);
                        
                        int y_479;
                        
                        y_479 = i_478 * size_92;
                        
                        int i_remnant_480;
                        
                        i_remnant_480 = flat_index_477 - y_479;
                        
                        double x_2266;
                        
                        x_2266 = *(double*) &mem_771[(i_remnant_480 * size_155 +
                                                      i_478 * 1) * 8 + 0];
                        
                        char y_2267;
                        
                        y_2267 = x_2266 <= 8288.0;
                        
                        char x_2268;
                        
                        x_2268 = x_2265 || y_2267;
                        
                        int flat_index_486;
                        
                        flat_index_486 = x_427 + 2;
                        
                        int i_487;
                        
                        i_487 = flat_index_486 / size_92 - (flat_index_486 %
                                                            size_92 != 0 &&
                                                            flat_index_486 %
                                                            size_92 < 0 !=
                                                            size_92 <
                                                            0 ? 1 : 0);
                        
                        int y_488;
                        
                        y_488 = i_487 * size_92;
                        
                        int i_remnant_489;
                        
                        i_remnant_489 = flat_index_486 - y_488;
                        
                        double x_2269;
                        
                        x_2269 = *(double*) &mem_771[(i_remnant_489 * size_155 +
                                                      i_487 * 1) * 8 + 0];
                        
                        char y_2270;
                        
                        y_2270 = x_2269 <= 840.0;
                        
                        char res_2271;
                        
                        res_2271 = x_2268 || y_2270;
                        
                        char res_2272;
                        
                        res_2272 = acc_276 || res_2271;
                        
                        char acc_tmp_537;
                        
                        acc_tmp_537 = res_2272;
                        acc_276 = acc_tmp_537;
                    }
                    res_2273 = acc_276;
                    
                    double x_2274;
                    
                    x_2274 = *(double*) &mem_771[(i_remnant_288 * size_155 +
                                                  i_283 * 1) * 8 + 0];
                    
                    char x_2275;
                    
                    x_2275 = x_2274 < 3758.05;
                    
                    double x_2276;
                    
                    x_2276 = *(double*) &mem_771[(i_remnant_296 * size_155 +
                                                  i_291 * 1) * 8 + 0];
                    
                    char y_2277;
                    
                    y_2277 = x_2276 < 1200.0;
                    
                    char x_2278;
                    
                    x_2278 = x_2275 || y_2277;
                    
                    double x_2279;
                    
                    x_2279 = *(double*) &mem_771[(i_remnant_300 * size_155 +
                                                  i_298 * 1) * 8 + 0];
                    
                    char y_2280;
                    
                    y_2280 = x_2279 < 11840.0;
                    
                    char y_2281;
                    
                    y_2281 = x_2278 || y_2280;
                    
                    char res_2282;
                    
                    res_2282 = res_2273 && y_2281;
                    
                    double x_2283;
                    
                    x_2283 = x_2276 / 1200.0;
                    
                    double arg_2284;
                    
                    arg_2284 = x_2283 - 1.0;
                    
                    double x_2285;
                    
                    x_2285 = x_2274 / 3758.05;
                    
                    double arg_2286;
                    
                    arg_2286 = x_2285 - 1.0;
                    
                    char cond_2287;
                    
                    cond_2287 = arg_2284 < arg_2286;
                    
                    double res_2288;
                    
                    if (cond_2287) {
                        res_2288 = arg_2284;
                    } else {
                        res_2288 = arg_2286;
                    }
                    
                    double x_2289;
                    
                    x_2289 = x_2279 / 11840.0;
                    
                    double arg_2290;
                    
                    arg_2290 = x_2289 - 1.0;
                    
                    char cond_2291;
                    
                    cond_2291 = arg_2290 < res_2288;
                    
                    double res_2292;
                    
                    if (cond_2291) {
                        res_2292 = arg_2290;
                    } else {
                        res_2292 = res_2288;
                    }
                    
                    double y_2293;
                    
                    y_2293 = 1.0 + res_2292;
                    
                    double res_2294;
                    
                    res_2294 = 1000.0 * y_2293;
                    
                    double res_2301;
                    
                    if (res_2282) {
                        int bounds_check_2295;
                        
                        if (!y_219) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_219",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                            abort();
                        }
                        
                        double y_2296;
                        
                        y_2296 = *(double*) &all_md_disc_mem_368[(0 * size_113 +
                                                                  1 * 1) * 8 +
                                                                 0];
                        
                        double res_2297;
                        
                        res_2297 = res_2294 * y_2296;
                        res_2301 = res_2297;
                    } else {
                        int bounds_check_2298;
                        
                        if (!y_219) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_219",
                                    "/home/athas/repos/finpar/benchmarks/OptionPricing/implementations/futhark/GenericPricing.fut:331:66-331:70");
                            abort();
                        }
                        
                        double y_2299;
                        
                        y_2299 = *(double*) &all_md_disc_mem_368[(0 * size_113 +
                                                                  1 * 1) * 8 +
                                                                 0];
                        
                        double res_2300;
                        
                        res_2300 = 1000.0 * y_2299;
                        res_2301 = res_2300;
                    }
                    
                    double res_2302;
                    
                    res_2302 = res_1862 + res_2301;
                    
                    double res_2303;
                    
                    res_2303 = acc_279 + res_2302;
                    
                    double acc_tmp_485;
                    
                    acc_tmp_485 = res_2303;
                    acc_279 = acc_tmp_485;
                }
                res_2304 = acc_279;
                res_2305 = res_2304;
            } else {
                res_2305 = 0.0;
            }
            res_2306 = res_2305;
        }
        res_2307 = res_2306;
    }
    
    double y_2308;
    
    y_2308 = futhark_toFloat64(num_mc_it_119);
    
    double arr_elem_2309;
    
    arr_elem_2309 = res_2307 / y_2308;
    
    unsigned char* mem_799;
    
    mem_799 = malloc(8);
    *(double*) &mem_799[0 * 1 * 8 + 0] = arr_elem_2309;
    out_mem_23 = mem_799;
    out_memsize_24 = 8;
    out_arrsize_26 = 1;
    
    struct tuple_int_mem_int retval_0;
    
    retval_0.elem_0 = out_memsize_24;
    retval_0.elem_1 = out_mem_23;
    retval_0.elem_2 = out_arrsize_26;
    return retval_0;
}
int main(int argc, char** argv)
{
    struct timeval t_start,  t_end,  t_diff;
    unsigned long elapsed_usec;
    
    {
        int out_memsize_24;
        unsigned char* out_mem_23;
        int out_arrsize_26;
        int dir_vs_mem_size_27;
        int all_md_c_mem_size_71;
        int all_md_vols_mem_size_255;
        int all_md_drifts_mem_size_345;
        int all_md_st_mem_size_354;
        int all_md_detvals_mem_size_360;
        int all_md_disc_mem_size_363;
        int bb_inds_mem_size_369;
        int bb_data_mem_size_376;
        unsigned char* dir_vs_mem_40;
        unsigned char* all_md_c_mem_84;
        unsigned char* all_md_vols_mem_268;
        unsigned char* all_md_drifts_mem_350;
        unsigned char* all_md_st_mem_359;
        unsigned char* all_md_detvals_mem_361;
        unsigned char* all_md_disc_mem_368;
        unsigned char* bb_inds_mem_371;
        unsigned char* bb_data_mem_381;
        int size_14;
        int size_18;
        int size_19;
        int size_20;
        int size_21;
        int size_22;
        int size_86;
        int size_89;
        int size_90;
        int size_91;
        int size_92;
        int size_93;
        int size_94;
        int size_95;
        int size_111;
        int size_112;
        int size_113;
        int size_114;
        int size_115;
        int size_116;
        int size_117;
        int contract_number_118;
        int num_mc_it_119;
        int num_dates_120;
        int num_und_121;
        int num_models_122;
        int num_bits_123;
        struct tuple_int_mem_int main_ret_1;
        
        if (read_int(&contract_number_118) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&num_mc_it_119) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&num_dates_120) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&num_und_121) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&num_models_122) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&num_bits_123) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(int), read_int, (void**) &dir_vs_mem_40,
                           shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n", "[[int]]");
                exit(1);
            }
            size_14 = shape[0];
            size_18 = shape[1];
            dir_vs_mem_size_27 = sizeof(int) * shape[0] * shape[1];
        }
        {
            int64_t shape[3];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_c_mem_84, shape, 3) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[[float64]]]");
                exit(1);
            }
            size_19 = shape[0];
            size_20 = shape[1];
            size_21 = shape[2];
            all_md_c_mem_size_71 = sizeof(double) * shape[0] * shape[1] *
                shape[2];
        }
        {
            int64_t shape[3];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_vols_mem_268, shape, 3) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[[float64]]]");
                exit(1);
            }
            size_22 = shape[0];
            size_86 = shape[1];
            size_89 = shape[2];
            all_md_vols_mem_size_255 = sizeof(double) * shape[0] * shape[1] *
                shape[2];
        }
        {
            int64_t shape[3];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_drifts_mem_350, shape, 3) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[[float64]]]");
                exit(1);
            }
            size_90 = shape[0];
            size_91 = shape[1];
            size_92 = shape[2];
            all_md_drifts_mem_size_345 = sizeof(double) * shape[0] * shape[1] *
                shape[2];
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_st_mem_359, shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[float64]]");
                exit(1);
            }
            size_93 = shape[0];
            size_94 = shape[1];
            all_md_st_mem_size_354 = sizeof(double) * shape[0] * shape[1];
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_detvals_mem_361, shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[float64]]");
                exit(1);
            }
            size_95 = shape[0];
            size_111 = shape[1];
            all_md_detvals_mem_size_360 = sizeof(double) * shape[0] * shape[1];
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &all_md_disc_mem_368, shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[float64]]");
                exit(1);
            }
            size_112 = shape[0];
            size_113 = shape[1];
            all_md_disc_mem_size_363 = sizeof(double) * shape[0] * shape[1];
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(int), read_int, (void**) &bb_inds_mem_371,
                           shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n", "[[int]]");
                exit(1);
            }
            size_114 = shape[0];
            size_115 = shape[1];
            bb_inds_mem_size_369 = sizeof(int) * shape[0] * shape[1];
        }
        {
            int64_t shape[2];
            
            if (read_array(sizeof(double), read_double,
                           (void**) &bb_data_mem_381, shape, 2) != 0) {
                fprintf(stderr, "Syntax error when reading %s.\n",
                        "[[float64]]");
                exit(1);
            }
            size_116 = shape[0];
            size_117 = shape[1];
            bb_data_mem_size_376 = sizeof(double) * shape[0] * shape[1];
        }
        gettimeofday(&t_start, NULL);
        main_ret_1 = futhark_main(dir_vs_mem_size_27, all_md_c_mem_size_71,
                                  all_md_vols_mem_size_255,
                                  all_md_drifts_mem_size_345,
                                  all_md_st_mem_size_354,
                                  all_md_detvals_mem_size_360,
                                  all_md_disc_mem_size_363,
                                  bb_inds_mem_size_369, bb_data_mem_size_376,
                                  dir_vs_mem_40, all_md_c_mem_84,
                                  all_md_vols_mem_268, all_md_drifts_mem_350,
                                  all_md_st_mem_359, all_md_detvals_mem_361,
                                  all_md_disc_mem_368, bb_inds_mem_371,
                                  bb_data_mem_381, size_14, size_18, size_19,
                                  size_20, size_21, size_22, size_86, size_89,
                                  size_90, size_91, size_92, size_93, size_94,
                                  size_95, size_111, size_112, size_113,
                                  size_114, size_115, size_116, size_117,
                                  contract_number_118, num_mc_it_119,
                                  num_dates_120, num_und_121, num_models_122,
                                  num_bits_123);
        gettimeofday(&t_end, NULL);
        out_memsize_24 = main_ret_1.elem_0;
        out_mem_23 = main_ret_1.elem_1;
        out_arrsize_26 = main_ret_1.elem_2;
        {
            if (out_arrsize_26 == 0) {
                printf("empty(%s)", "float64");
            } else {
                int print_i_2;
                
                putchar('[');
                for (print_i_2 = 0; print_i_2 < out_arrsize_26; print_i_2++) {
                    double* print_elem_3 = (double*) out_mem_23 + print_i_2 * 1;
                    
                    printf("%.6f", *print_elem_3);
                    if (print_i_2 != out_arrsize_26 - 1) {
                        printf(", ");
                    }
                }
                putchar(']');
            }
        }
        printf("\n");
    }
    ;
    if (argc == 3 && strcmp(argv[1], "-t") == 0) {
        FILE * runtime_file;
        runtime_file = fopen(argv[2], "w");
        if (runtime_file == NULL) {
            fprintf(stderr, "Cannot open %s: %s\n", argv[2], strerror(errno));
            exit(1);
        }
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed_usec = t_diff.tv_sec * 1e6 + t_diff.tv_usec;
        fprintf(runtime_file, "%ld\n", elapsed_usec / 1000);
        fclose(runtime_file);
    }
    return 0;
}
