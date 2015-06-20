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
static struct tuple_int_mem_int futhark_main(int outer_loop_count_8,
                                             int numX_12, int numY_13,
                                             int numT_14, double s0_15,
                                             double strike_17, double t_25,
                                             double alpha_26, double nu_31,
                                             double beta_32);
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
static struct tuple_int_mem_int futhark_main(int outer_loop_count_8,
                                             int numX_12, int numY_13,
                                             int numT_14, double s0_15,
                                             double strike_17, double t_25,
                                             double alpha_26, double nu_31,
                                             double beta_32)
{
    int out_memsize_78;
    unsigned char* out_mem_77;
    int out_arrsize_98;
    double res_37;
    
    res_37 = futhark_log64(alpha_26);
    
    double x_61;
    
    x_61 = futhark_toFloat64(numT_14);
    
    double y_97;
    
    y_97 = x_61 - 1.0;
    
    int bytes_42;
    
    bytes_42 = 8 * numT_14;
    
    unsigned char* mem_51;
    
    mem_51 = malloc(bytes_42);
    for (int i_0 = 0; i_0 < numT_14; i_0++) {
        double y_99;
        
        y_99 = futhark_toFloat64(i_0);
        
        double x_100;
        
        x_100 = t_25 * y_99;
        
        double res_101;
        
        res_101 = x_100 / y_97;
        *(double*) &mem_51[i_0 * 1 * 8 + 0] = res_101;
    }
    
    double x_103;
    
    x_103 = 20.0 * alpha_26;
    
    double x_104;
    
    x_104 = x_103 * s0_15;
    
    double y_105;
    
    y_105 = futhark_sqrt64(t_25);
    
    double res_106;
    
    res_106 = x_104 * y_105;
    
    double x_107;
    
    x_107 = 10.0 * nu_31;
    
    double res_108;
    
    res_108 = x_107 * y_105;
    
    double y_109;
    
    y_109 = futhark_toFloat64(numX_12);
    
    double res_110;
    
    res_110 = res_106 / y_109;
    
    double y_111;
    
    y_111 = futhark_toFloat64(numY_13);
    
    double res_112;
    
    res_112 = res_108 / y_111;
    
    double arg_113;
    
    arg_113 = s0_15 / res_110;
    
    int res_115;
    
    res_115 = futhark_trunc64(arg_113);
    
    int res_116;
    
    res_116 = numY_13 / 2 - (numY_13 % 2 != 0 && numY_13 % 2 < 0 != 2 <
                             0 ? 1 : 0);
    
    double x_118;
    
    x_118 = futhark_toFloat64(res_115);
    
    double y_119;
    
    y_119 = x_118 * res_110;
    
    int bytes_157;
    
    bytes_157 = 8 * numX_12;
    
    unsigned char* mem_159;
    
    mem_159 = malloc(bytes_157);
    for (int i_4 = 0; i_4 < numX_12; i_4++) {
        double x_121;
        
        x_121 = futhark_toFloat64(i_4);
        
        double x_122;
        
        x_122 = x_121 * res_110;
        
        double x_123;
        
        x_123 = x_122 - y_119;
        
        double res_124;
        
        res_124 = x_123 + s0_15;
        *(double*) &mem_159[i_4 * 1 * 8 + 0] = res_124;
    }
    
    double x_127;
    
    x_127 = futhark_toFloat64(res_116);
    
    double y_128;
    
    y_128 = x_127 * res_112;
    
    int bytes_173;
    
    bytes_173 = 8 * numY_13;
    
    unsigned char* mem_174;
    
    mem_174 = malloc(bytes_173);
    for (int i_9 = 0; i_9 < numY_13; i_9++) {
        double x_130;
        
        x_130 = futhark_toFloat64(i_9);
        
        double x_131;
        
        x_131 = x_130 * res_112;
        
        double x_132;
        
        x_132 = x_131 - y_128;
        
        double res_133;
        
        res_133 = x_132 + res_37;
        *(double*) &mem_174[i_9 * 1 * 8 + 0] = res_133;
    }
    
    unsigned char* mem_186;
    
    mem_186 = malloc(24);
    *(double*) &mem_186[0 * 1 * 8 + 0] = 0.0;
    *(double*) &mem_186[1 * 1 * 8 + 0] = 0.0;
    *(double*) &mem_186[2 * 1 * 8 + 0] = 0.0;
    
    unsigned char* mem_190;
    
    mem_190 = malloc(24);
    for (int i_147 = 0; i_147 < 3; i_147++) {
        *(double*) &mem_190[8 * i_147] = *(double*) &mem_186[8 * i_147];
    }
    
    int n_215;
    
    n_215 = numX_12 - 2;
    
    int x_192;
    
    x_192 = 8 * n_215;
    
    int bytes_191;
    
    bytes_191 = x_192 * 3;
    
    unsigned char* mem_193;
    
    mem_193 = malloc(bytes_191);
    for (int i_18 = 0; i_18 < n_215; i_18++) {
        int res_220;
        
        res_220 = i_18 + 1;
        
        double x_337;
        
        x_337 = *(double*) &mem_159[res_220 * 1 * 8 + 0];
        
        int i_340;
        
        i_340 = res_220 - 1;
        
        double y_341;
        
        y_341 = *(double*) &mem_159[i_340 * 1 * 8 + 0];
        
        double res_342;
        
        res_342 = x_337 - y_341;
        
        int i_343;
        
        i_343 = res_220 + 1;
        
        double x_345;
        
        x_345 = *(double*) &mem_159[i_343 * 1 * 8 + 0];
        
        double res_349;
        
        res_349 = x_345 - x_337;
        
        double y_387;
        
        y_387 = res_342 + res_349;
        
        double x_423;
        
        x_423 = 2.0 / res_342;
        
        double arr_elem_424;
        
        arr_elem_424 = x_423 / y_387;
        
        double x_425;
        
        x_425 = 1.0 / res_342;
        
        double y_426;
        
        y_426 = 1.0 / res_349;
        
        double y_427;
        
        y_427 = x_425 + y_426;
        
        double x_428;
        
        x_428 = -2.0 * y_427;
        
        double arr_elem_429;
        
        arr_elem_429 = x_428 / y_387;
        
        double x_430;
        
        x_430 = 2.0 / res_349;
        
        double arr_elem_431;
        
        arr_elem_431 = x_430 / y_387;
        *(double*) &mem_193[(i_18 * 3 + 0 * 1) * 8 + 0] = arr_elem_424;
        *(double*) &mem_193[(i_18 * 3 + 1 * 1) * 8 + 0] = arr_elem_429;
        *(double*) &mem_193[(i_18 * 3 + 2 * 1) * 8 + 0] = arr_elem_431;
    }
    
    int i_435;
    
    i_435 = numX_12 - 1;
    
    char x_436;
    
    x_436 = 0 <= i_435;
    
    int bounds_check_437;
    
    if (!x_436) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "x_436",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:27:18-27:19");
        abort();
    }
    
    int conc_tmp_447;
    
    conc_tmp_447 = 1 + n_215;
    
    int conc_tmp_449;
    
    conc_tmp_449 = conc_tmp_447 + 1;
    
    int x_209;
    
    x_209 = 8 * conc_tmp_447;
    
    int bytes_204;
    
    bytes_204 = x_209 * 3;
    
    unsigned char* mem_210;
    
    mem_210 = malloc(bytes_204);
    
    int tmp_offs_153;
    
    tmp_offs_153 = 0;
    memmove(mem_210 + (0 + tmp_offs_153 * 3 * 8), mem_190 + 0, 3 *
            sizeof(double));
    tmp_offs_153 = tmp_offs_153 + 1;
    memmove(mem_210 + (0 + tmp_offs_153 * 3 * 8), mem_193 + 0, n_215 * 3 *
            sizeof(double));
    tmp_offs_153 = tmp_offs_153 + n_215;
    
    int x_212;
    
    x_212 = 8 * conc_tmp_449;
    
    int bytes_211;
    
    bytes_211 = x_212 * 3;
    
    unsigned char* mem_216;
    
    mem_216 = malloc(bytes_211);
    
    int tmp_offs_154;
    
    tmp_offs_154 = 0;
    memmove(mem_216 + (0 + tmp_offs_154 * 3 * 8), mem_210 + 0, conc_tmp_447 *
            3 * sizeof(double));
    tmp_offs_154 = tmp_offs_154 + conc_tmp_447;
    memmove(mem_216 + (0 + tmp_offs_154 * 3 * 8), mem_190 + 0, 3 *
            sizeof(double));
    tmp_offs_154 = tmp_offs_154 + 1;
    
    int n_464;
    
    n_464 = numY_13 - 2;
    
    int x_218;
    
    x_218 = 8 * n_464;
    
    int bytes_217;
    
    bytes_217 = x_218 * 3;
    
    unsigned char* mem_219;
    
    mem_219 = malloc(bytes_217);
    for (int i_22 = 0; i_22 < n_464; i_22++) {
        int res_467;
        
        res_467 = i_22 + 1;
        
        double x_468;
        
        x_468 = *(double*) &mem_174[res_467 * 1 * 8 + 0];
        
        int i_469;
        
        i_469 = res_467 - 1;
        
        double y_470;
        
        y_470 = *(double*) &mem_174[i_469 * 1 * 8 + 0];
        
        double res_471;
        
        res_471 = x_468 - y_470;
        
        int i_472;
        
        i_472 = res_467 + 1;
        
        double x_473;
        
        x_473 = *(double*) &mem_174[i_472 * 1 * 8 + 0];
        
        double res_474;
        
        res_474 = x_473 - x_468;
        
        double y_480;
        
        y_480 = res_471 + res_474;
        
        double x_489;
        
        x_489 = 2.0 / res_471;
        
        double arr_elem_490;
        
        arr_elem_490 = x_489 / y_480;
        
        double x_491;
        
        x_491 = 1.0 / res_471;
        
        double y_494;
        
        y_494 = 1.0 / res_474;
        
        double y_495;
        
        y_495 = x_491 + y_494;
        
        double x_496;
        
        x_496 = -2.0 * y_495;
        
        double arr_elem_497;
        
        arr_elem_497 = x_496 / y_480;
        
        double x_498;
        
        x_498 = 2.0 / res_474;
        
        double arr_elem_499;
        
        arr_elem_499 = x_498 / y_480;
        *(double*) &mem_219[(i_22 * 3 + 0 * 1) * 8 + 0] = arr_elem_490;
        *(double*) &mem_219[(i_22 * 3 + 1 * 1) * 8 + 0] = arr_elem_497;
        *(double*) &mem_219[(i_22 * 3 + 2 * 1) * 8 + 0] = arr_elem_499;
    }
    
    int i_503;
    
    i_503 = numY_13 - 1;
    
    char x_504;
    
    x_504 = 0 <= i_503;
    
    int bounds_check_505;
    
    if (!x_504) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "x_504",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:27:18-27:19");
        abort();
    }
    
    int conc_tmp_515;
    
    conc_tmp_515 = 1 + n_464;
    
    int conc_tmp_517;
    
    conc_tmp_517 = conc_tmp_515 + 1;
    
    int x_226;
    
    x_226 = 8 * conc_tmp_515;
    
    int bytes_225;
    
    bytes_225 = x_226 * 3;
    
    unsigned char* mem_227;
    
    mem_227 = malloc(bytes_225);
    
    int tmp_offs_156;
    
    tmp_offs_156 = 0;
    memmove(mem_227 + (0 + tmp_offs_156 * 3 * 8), mem_190 + 0, 3 *
            sizeof(double));
    tmp_offs_156 = tmp_offs_156 + 1;
    memmove(mem_227 + (0 + tmp_offs_156 * 3 * 8), mem_219 + 0, n_464 * 3 *
            sizeof(double));
    tmp_offs_156 = tmp_offs_156 + n_464;
    
    int x_229;
    
    x_229 = 8 * conc_tmp_517;
    
    int bytes_228;
    
    bytes_228 = x_229 * 3;
    
    unsigned char* mem_230;
    
    mem_230 = malloc(bytes_228);
    
    int tmp_offs_158;
    
    tmp_offs_158 = 0;
    memmove(mem_230 + (0 + tmp_offs_158 * 3 * 8), mem_227 + 0, conc_tmp_515 *
            3 * sizeof(double));
    tmp_offs_158 = tmp_offs_158 + conc_tmp_515;
    memmove(mem_230 + (0 + tmp_offs_158 * 3 * 8), mem_190 + 0, 3 *
            sizeof(double));
    tmp_offs_158 = tmp_offs_158 + 1;
    
    int bound_521;
    
    bound_521 = numT_14 - 1;
    
    int x_522;
    
    x_522 = numT_14 - 2;
    
    double replicate_v_524;
    
    replicate_v_524 = nu_31 * nu_31;
    
    double x_527;
    
    x_527 = 0.5 * nu_31;
    
    double x_528;
    
    x_528 = x_527 * nu_31;
    
    char zip_cmp_530;
    
    zip_cmp_530 = conc_tmp_449 == numX_12;
    
    int zip_assert_531;
    
    if (!zip_cmp_530) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_530",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:165:21-165:24");
        abort();
    }
    
    int y_532;
    
    y_532 = conc_tmp_449 - 1;
    
    int y_536;
    
    y_536 = conc_tmp_517 - 1;
    
    char zip_cmp_537;
    
    zip_cmp_537 = numY_13 == conc_tmp_517;
    
    int zip_assert_538;
    
    if (!zip_cmp_537) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_537",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:219:18-219:21");
        abort();
    }
    
    char zip_cmp_542;
    
    zip_cmp_542 = numX_12 == conc_tmp_449;
    
    int zip_assert_543;
    
    if (!zip_cmp_542) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "zip_cmp_542",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:188:25-188:28");
        abort();
    }
    
    char x_549;
    
    x_549 = 0 <= res_116;
    
    char x_555;
    
    x_555 = 0 <= res_115;
    
    char y_556;
    
    y_556 = res_115 < numX_12;
    
    char assert_arg_557;
    
    assert_arg_557 = x_555 && y_556;
    
    int bounds_check_558;
    
    if (!assert_arg_557) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_557",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:254:5-254:13");
        abort();
    }
    
    double x_566;
    
    x_566 = 0.5 * replicate_v_524;
    
    char y_569;
    
    y_569 = res_116 < numY_13;
    
    char assert_arg_570;
    
    assert_arg_570 = x_549 && y_569;
    
    int bounds_check_571;
    
    if (!assert_arg_570) {
        fprintf(stderr, "Assertion %s at %s failed.\n", "assert_arg_570",
                "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:254:5-254:13");
        abort();
    }
    
    int bytes_231;
    
    bytes_231 = 8 * outer_loop_count_8;
    
    unsigned char* mem_232;
    
    mem_232 = malloc(bytes_231);
    
    int bytes_235;
    
    bytes_235 = bytes_157 * numY_13;
    
    unsigned char* mem_237;
    
    mem_237 = malloc(bytes_235);
    
    int bytes_242;
    
    bytes_242 = bytes_173 * numX_12;
    
    unsigned char* mem_244;
    
    mem_244 = malloc(bytes_242);
    
    int bytes_247;
    
    bytes_247 = bytes_173 * conc_tmp_449;
    
    unsigned char* mem_249;
    
    mem_249 = malloc(bytes_247);
    
    unsigned char* mem_253;
    
    mem_253 = malloc(bytes_157);
    
    unsigned char* mem_266;
    
    mem_266 = malloc(bytes_247);
    
    int bytes_267;
    
    bytes_267 = bytes_157 * conc_tmp_517;
    
    unsigned char* mem_269;
    
    mem_269 = malloc(bytes_267);
    
    unsigned char* mem_280;
    
    mem_280 = malloc(bytes_267);
    
    unsigned char* mem_283;
    
    mem_283 = malloc(bytes_247);
    
    unsigned char* mem_287;
    
    mem_287 = malloc(x_212);
    
    unsigned char* mem_293;
    
    mem_293 = malloc(x_212);
    
    unsigned char* mem_295;
    
    mem_295 = malloc(bytes_157);
    
    unsigned char* mem_297;
    
    mem_297 = malloc(bytes_157);
    
    unsigned char* mem_299;
    
    mem_299 = malloc(bytes_157);
    
    unsigned char* mem_313;
    
    mem_313 = malloc(bytes_157);
    
    int bytes_329;
    
    bytes_329 = x_212 * numY_13;
    
    unsigned char* mem_331;
    
    mem_331 = malloc(bytes_329);
    
    unsigned char* mem_348;
    
    mem_348 = malloc(bytes_329);
    
    unsigned char* mem_351;
    
    mem_351 = malloc(bytes_173);
    
    unsigned char* mem_353;
    
    mem_353 = malloc(bytes_173);
    
    unsigned char* mem_355;
    
    mem_355 = malloc(bytes_173);
    
    unsigned char* mem_370;
    
    mem_370 = malloc(bytes_235);
    
    unsigned char* mem_374;
    
    mem_374 = malloc(bytes_173);
    
    unsigned char* mem_376;
    
    mem_376 = malloc(bytes_173);
    
    unsigned char* double_buffer_mem_41;
    
    double_buffer_mem_41 = malloc(bytes_242);
    
    unsigned char* mem_395;
    
    mem_395 = malloc(bytes_242);
    for (int i_136 = 0; i_136 < outer_loop_count_8; i_136++) {
        double y_575;
        
        y_575 = futhark_toFloat64(i_136);
        
        double res_578;
        
        res_578 = 1.0e-3 * y_575;
        for (int i_28 = 0; i_28 < numX_12; i_28++) {
            double xi_581;
            
            xi_581 = *(double*) &mem_159[i_28 * 1 * 8 + 0];
            
            double arg_587;
            
            arg_587 = xi_581 - res_578;
            
            char cond_588;
            
            cond_588 = 0.0 < arg_587;
            
            double res_591;
            
            if (cond_588) {
                res_591 = arg_587;
            } else {
                res_591 = 0.0;
            }
            for (int i_166 = 0; i_166 < numY_13; i_166++) {
                *(double*) &mem_237[(i_28 * numY_13 + i_166 * 1) * 8 + 0] =
                    res_591;
            }
        }
        for (int i_167 = 0; i_167 < numY_13; i_167++) {
            for (int i_168 = 0; i_168 < numX_12; i_168++) {
                *(double*) &mem_244[8 * (numX_12 * i_167) + 8 * i_168] =
                    *(double*) &mem_237[8 * (numY_13 * i_168) + 8 * i_167];
            }
        }
        
        unsigned char* myResult_mem_398;
        unsigned char* myResult_mem_246;
        
        myResult_mem_246 = mem_244;
        for (int i_601 = 0; i_601 < bound_521; i_601++) {
            int res_602;
            
            res_602 = x_522 - i_601;
            
            double y_603;
            
            y_603 = *(double*) &mem_51[res_602 * 1 * 8 + 0];
            
            double y_604;
            
            y_604 = x_528 * y_603;
            
            int i_605;
            
            i_605 = res_602 + 1;
            
            double x_606;
            
            x_606 = *(double*) &mem_51[i_605 * 1 * 8 + 0];
            
            double y_618;
            
            y_618 = x_606 - y_603;
            
            double res_619;
            
            res_619 = 1.0 / y_618;
            for (int i_44 = 0; i_44 < numY_13; i_44++) {
                double yj_620;
                
                yj_620 = *(double*) &mem_174[i_44 * 1 * 8 + 0];
                for (int i_34 = 0; i_34 < numX_12; i_34++) {
                    double xi_639;
                    
                    xi_639 = *(double*) &mem_159[i_34 * 1 * 8 + 0];
                    
                    double y_640;
                    
                    y_640 = futhark_log64(xi_639);
                    
                    double x_642;
                    
                    x_642 = beta_32 * y_640;
                    
                    double x_643;
                    
                    x_643 = x_642 + yj_620;
                    
                    double y_644;
                    
                    y_644 = x_643 - y_604;
                    
                    double arg_645;
                    
                    arg_645 = 2.0 * y_644;
                    
                    double res_646;
                    
                    res_646 = futhark_exp64(arg_645);
                    *(double*) &mem_253[i_34 * 1 * 8 + 0] = res_646;
                }
                for (int i_40 = 0; i_40 < conc_tmp_449; i_40++) {
                    double tup_653;
                    
                    tup_653 = *(double*) &mem_253[i_40 * 1 * 8 + 0];
                    
                    char cond_657;
                    
                    cond_657 = 0 < i_40;
                    
                    double x_658;
                    
                    x_658 = 0.5 * tup_653;
                    
                    int i_660;
                    
                    i_660 = i_40 - 1;
                    
                    double res_668;
                    
                    if (cond_657) {
                        double y_662;
                        
                        y_662 = *(double*) &mem_216[(i_40 * 3 + 0 * 1) * 8 + 0];
                        
                        double y_663;
                        
                        y_663 = x_658 * y_662;
                        
                        double y_665;
                        
                        y_665 = *(double*) &myResult_mem_246[(i_44 * numX_12 +
                                                              i_660 * 1) * 8 +
                                                             0];
                        
                        double res_666;
                        
                        res_666 = y_663 * y_665;
                        res_668 = res_666;
                    } else {
                        res_668 = 0.0;
                    }
                    
                    char cond_669;
                    
                    cond_669 = i_40 < y_532;
                    
                    int i_671;
                    
                    i_671 = i_40 + 1;
                    
                    char y_672;
                    
                    y_672 = i_671 < numX_12;
                    
                    double res_679;
                    
                    if (cond_669) {
                        double y_674;
                        
                        y_674 = *(double*) &mem_216[(i_40 * 3 + 2 * 1) * 8 + 0];
                        
                        double y_675;
                        
                        y_675 = x_658 * y_674;
                        
                        int bounds_check_676;
                        
                        if (!y_672) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_672",
                                    "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:161:69-161:79");
                            abort();
                        }
                        
                        double y_677;
                        
                        y_677 = *(double*) &myResult_mem_246[(i_44 * numX_12 +
                                                              i_671 * 1) * 8 +
                                                             0];
                        
                        double res_678;
                        
                        res_678 = y_675 * y_677;
                        res_679 = res_678;
                    } else {
                        res_679 = 0.0;
                    }
                    
                    double y_682;
                    
                    y_682 = *(double*) &mem_216[(i_40 * 3 + 1 * 1) * 8 + 0];
                    
                    double y_683;
                    
                    y_683 = x_658 * y_682;
                    
                    double y_686;
                    
                    y_686 = *(double*) &myResult_mem_246[(i_44 * numX_12 +
                                                          i_40 * 1) * 8 + 0];
                    
                    double res_687;
                    
                    res_687 = y_683 * y_686;
                    
                    double x_688;
                    
                    x_688 = res_668 + res_687;
                    
                    double res_689;
                    
                    res_689 = x_688 + res_679;
                    *(double*) &mem_249[(i_44 * conc_tmp_449 + i_40 * 1) * 8 +
                                        0] = res_689;
                }
            }
            memmove(mem_266 + 0, mem_249 + 0, numY_13 * conc_tmp_449 *
                    sizeof(double));
            for (int i_53 = 0; i_53 < numX_12; i_53++) {
                for (int i_49 = 0; i_49 < conc_tmp_517; i_49++) {
                    char cond_707;
                    
                    cond_707 = 0 < i_49;
                    
                    int i_709;
                    
                    i_709 = i_49 - 1;
                    
                    double res_748;
                    
                    if (cond_707) {
                        double y_711;
                        
                        y_711 = *(double*) &mem_230[(i_49 * 3 + 0 * 1) * 8 + 0];
                        
                        double y_717;
                        
                        y_717 = x_566 * y_711;
                        
                        double y_746;
                        
                        y_746 = *(double*) &myResult_mem_246[(i_709 * numX_12 +
                                                              i_53 * 1) * 8 +
                                                             0];
                        
                        double res_747;
                        
                        res_747 = y_717 * y_746;
                        res_748 = res_747;
                    } else {
                        res_748 = 0.0;
                    }
                    
                    char cond_749;
                    
                    cond_749 = i_49 < y_536;
                    
                    int i_750;
                    
                    i_750 = i_49 + 1;
                    
                    char y_751;
                    
                    y_751 = i_750 < numY_13;
                    
                    double res_776;
                    
                    if (cond_749) {
                        double y_754;
                        
                        y_754 = *(double*) &mem_230[(i_49 * 3 + 2 * 1) * 8 + 0];
                        
                        double y_763;
                        
                        y_763 = x_566 * y_754;
                        
                        int bounds_check_764;
                        
                        if (!y_751) {
                            fprintf(stderr, "Assertion %s at %s failed.\n",
                                    "y_751",
                                    "/home/athas/repos/finpar/benchmarks/LocVolCalib/implementations/futhark/CalibVolDiff.fut:161:69-161:79");
                            abort();
                        }
                        
                        double y_767;
                        
                        y_767 = *(double*) &myResult_mem_246[(i_750 * numX_12 +
                                                              i_53 * 1) * 8 +
                                                             0];
                        
                        double res_771;
                        
                        res_771 = y_763 * y_767;
                        res_776 = res_771;
                    } else {
                        res_776 = 0.0;
                    }
                    
                    double y_782;
                    
                    y_782 = *(double*) &mem_230[(i_49 * 3 + 1 * 1) * 8 + 0];
                    
                    double y_792;
                    
                    y_792 = x_566 * y_782;
                    
                    double y_796;
                    
                    y_796 = *(double*) &myResult_mem_246[(i_49 * numX_12 +
                                                          i_53 * 1) * 8 + 0];
                    
                    double res_797;
                    
                    res_797 = y_792 * y_796;
                    
                    double x_805;
                    
                    x_805 = res_748 + res_797;
                    
                    double res_809;
                    
                    res_809 = x_805 + res_776;
                    *(double*) &mem_269[(i_53 * conc_tmp_517 + i_49 * 1) * 8 +
                                        0] = res_809;
                }
            }
            memmove(mem_280 + 0, mem_269 + 0, numX_12 * conc_tmp_517 *
                    sizeof(double));
            for (int i_72 = 0; i_72 < numY_13; i_72++) {
                double yj_846;
                
                yj_846 = *(double*) &mem_174[i_72 * 1 * 8 + 0];
                
                int x_140;
                
                x_140 = i_72 * numX_12;
                for (int i_57 = 0; i_57 < conc_tmp_449; i_57++) {
                    double tup_906;
                    
                    tup_906 = *(double*) &mem_266[(i_72 * conc_tmp_449 + i_57 *
                                                   1) * 8 + 0];
                    
                    double tup_907;
                    
                    tup_907 = *(double*) &myResult_mem_246[(i_72 * numX_12 +
                                                            i_57 * 1) * 8 + 0];
                    
                    int flat_index_142;
                    
                    flat_index_142 = x_140 + i_57;
                    
                    int i_143;
                    
                    i_143 = flat_index_142 / numX_12 - (flat_index_142 %
                                                        numX_12 != 0 &&
                                                        flat_index_142 %
                                                        numX_12 < 0 != numX_12 <
                                                        0 ? 1 : 0);
                    
                    int y_144;
                    
                    y_144 = i_143 * numX_12;
                    
                    int i_remnant_145;
                    
                    i_remnant_145 = flat_index_142 - y_144;
                    
                    double y_908;
                    
                    y_908 = *(double*) &mem_280[(i_remnant_145 * conc_tmp_517 +
                                                 i_143 * 1) * 8 + 0];
                    
                    double x_909;
                    
                    x_909 = res_619 * tup_907;
                    
                    double y_910;
                    
                    y_910 = 0.5 * tup_906;
                    
                    double res_911;
                    
                    res_911 = x_909 + y_910;
                    
                    double res_912;
                    
                    res_912 = res_911 + y_908;
                    *(double*) &mem_287[i_57 * 1 * 8 + 0] = res_912;
                }
                memmove(mem_293 + 0, mem_287 + 0, conc_tmp_449 *
                        sizeof(double));
                for (int i_62 = 0; i_62 < numX_12; i_62++) {
                    double xi_916;
                    
                    xi_916 = *(double*) &mem_159[i_62 * 1 * 8 + 0];
                    
                    double y_917;
                    
                    y_917 = futhark_log64(xi_916);
                    
                    double x_918;
                    
                    x_918 = beta_32 * y_917;
                    
                    double x_919;
                    
                    x_919 = x_918 + yj_846;
                    
                    double y_920;
                    
                    y_920 = x_919 - y_604;
                    
                    double arg_921;
                    
                    arg_921 = 2.0 * y_920;
                    
                    double res_922;
                    
                    res_922 = futhark_exp64(arg_921);
                    
                    double x_923;
                    
                    x_923 = 0.5 * res_922;
                    
                    int x_149;
                    
                    x_149 = i_62 * 3;
                    
                    int i_152;
                    
                    i_152 = x_149 / 3 - (x_149 % 3 != 0 && x_149 % 3 < 0 != 3 <
                                         0 ? 1 : 0);
                    
                    double y_924;
                    
                    y_924 = *(double*) &mem_216[(i_152 * 3 + 0 * 1) * 8 + 0];
                    
                    double y_925;
                    
                    y_925 = x_923 * y_924;
                    
                    double y_926;
                    
                    y_926 = 0.5 * y_925;
                    
                    double res_927;
                    
                    res_927 = 0.0 - y_926;
                    
                    int flat_index_160;
                    
                    flat_index_160 = x_149 + 1;
                    
                    int i_161;
                    
                    i_161 = flat_index_160 / 3 - (flat_index_160 % 3 != 0 &&
                                                  flat_index_160 % 3 < 0 != 3 <
                                                  0 ? 1 : 0);
                    
                    int y_162;
                    
                    y_162 = i_161 * 3;
                    
                    int i_remnant_163;
                    
                    i_remnant_163 = flat_index_160 - y_162;
                    
                    double y_928;
                    
                    y_928 = *(double*) &mem_216[(i_161 * 3 + i_remnant_163 *
                                                 1) * 8 + 0];
                    
                    double y_929;
                    
                    y_929 = x_923 * y_928;
                    
                    double y_930;
                    
                    y_930 = 0.5 * y_929;
                    
                    double res_931;
                    
                    res_931 = res_619 - y_930;
                    
                    int flat_index_169;
                    
                    flat_index_169 = x_149 + 2;
                    
                    int i_170;
                    
                    i_170 = flat_index_169 / 3 - (flat_index_169 % 3 != 0 &&
                                                  flat_index_169 % 3 < 0 != 3 <
                                                  0 ? 1 : 0);
                    
                    int y_171;
                    
                    y_171 = i_170 * 3;
                    
                    int i_remnant_172;
                    
                    i_remnant_172 = flat_index_169 - y_171;
                    
                    double y_932;
                    
                    y_932 = *(double*) &mem_216[(i_170 * 3 + i_remnant_172 *
                                                 1) * 8 + 0];
                    
                    double y_933;
                    
                    y_933 = x_923 * y_932;
                    
                    double y_934;
                    
                    y_934 = 0.5 * y_933;
                    
                    double res_935;
                    
                    res_935 = 0.0 - y_934;
                    *(double*) &mem_295[i_62 * 1 * 8 + 0] = res_927;
                    *(double*) &mem_297[i_62 * 1 * 8 + 0] = res_931;
                    *(double*) &mem_299[i_62 * 1 * 8 + 0] = res_935;
                }
                memmove(mem_313 + 0, mem_297 + 0, numX_12 * sizeof(double));
                for (int i_945 = 0; i_945 < i_435; i_945++) {
                    int res_946;
                    
                    res_946 = i_945 + 1;
                    
                    double x_947;
                    
                    x_947 = *(double*) &mem_295[res_946 * 1 * 8 + 0];
                    
                    int i_948;
                    
                    i_948 = res_946 - 1;
                    
                    double y_949;
                    
                    y_949 = *(double*) &mem_313[i_948 * 1 * 8 + 0];
                    
                    double res_950;
                    
                    res_950 = x_947 / y_949;
                    
                    double x_951;
                    
                    x_951 = *(double*) &mem_313[res_946 * 1 * 8 + 0];
                    
                    double y_952;
                    
                    y_952 = *(double*) &mem_299[i_948 * 1 * 8 + 0];
                    
                    double y_953;
                    
                    y_953 = res_950 * y_952;
                    
                    double lw_val_954;
                    
                    lw_val_954 = x_951 - y_953;
                    *(double*) &mem_313[res_946 * 1 * 8 + 0] = lw_val_954;
                    
                    double x_956;
                    
                    x_956 = *(double*) &mem_293[res_946 * 1 * 8 + 0];
                    
                    double y_957;
                    
                    y_957 = *(double*) &mem_293[i_948 * 1 * 8 + 0];
                    
                    double y_958;
                    
                    y_958 = res_950 * y_957;
                    
                    double lw_val_959;
                    
                    lw_val_959 = x_956 - y_958;
                    *(double*) &mem_293[res_946 * 1 * 8 + 0] = lw_val_959;
                }
                
                double x_963;
                
                x_963 = *(double*) &mem_293[i_435 * 1 * 8 + 0];
                
                double y_964;
                
                y_964 = *(double*) &mem_313[i_435 * 1 * 8 + 0];
                
                double lw_val_965;
                
                lw_val_965 = x_963 / y_964;
                *(double*) &mem_293[i_435 * 1 * 8 + 0] = lw_val_965;
                for (int i_222 = 0; i_222 < conc_tmp_449; i_222++) {
                    *(double*) &mem_283[8 * (i_72 * conc_tmp_449) + 8 * i_222] =
                        *(double*) &mem_293[8 * i_222];
                }
                for (int j_968 = 0; j_968 < i_435; j_968++) {
                    int res_969;
                    
                    res_969 = n_215 - j_968;
                    
                    double x_970;
                    
                    x_970 = *(double*) &mem_283[(i_72 * conc_tmp_449 + res_969 *
                                                 1) * 8 + 0];
                    
                    double x_971;
                    
                    x_971 = *(double*) &mem_299[res_969 * 1 * 8 + 0];
                    
                    int i_972;
                    
                    i_972 = res_969 + 1;
                    
                    double y_973;
                    
                    y_973 = *(double*) &mem_283[(i_72 * conc_tmp_449 + i_972 *
                                                 1) * 8 + 0];
                    
                    double y_974;
                    
                    y_974 = x_971 * y_973;
                    
                    double x_975;
                    
                    x_975 = x_970 - y_974;
                    
                    double y_976;
                    
                    y_976 = *(double*) &mem_313[res_969 * 1 * 8 + 0];
                    
                    double lw_val_977;
                    
                    lw_val_977 = x_975 / y_976;
                    *(double*) &mem_283[(i_72 * conc_tmp_449 + res_969 * 1) *
                                        8 + 0] = lw_val_977;
                }
            }
            for (int i_80 = 0; i_80 < conc_tmp_449; i_80++) {
                int x_176;
                
                x_176 = i_80 * conc_tmp_517;
                for (int i_76 = 0; i_76 < numY_13; i_76++) {
                    double uv_987;
                    
                    uv_987 = *(double*) &mem_283[(i_76 * conc_tmp_449 + i_80 *
                                                  1) * 8 + 0];
                    
                    int flat_index_178;
                    
                    flat_index_178 = x_176 + i_76;
                    
                    int i_179;
                    
                    i_179 = flat_index_178 / conc_tmp_517 - (flat_index_178 %
                                                             conc_tmp_517 !=
                                                             0 &&
                                                             flat_index_178 %
                                                             conc_tmp_517 < 0 !=
                                                             conc_tmp_517 <
                                                             0 ? 1 : 0);
                    
                    int y_180;
                    
                    y_180 = i_179 * conc_tmp_517;
                    
                    int i_remnant_181;
                    
                    i_remnant_181 = flat_index_178 - y_180;
                    
                    double uv_988;
                    
                    uv_988 = *(double*) &mem_280[(i_179 * conc_tmp_517 +
                                                  i_remnant_181 * 1) * 8 + 0];
                    
                    double x_989;
                    
                    x_989 = res_619 * uv_987;
                    
                    double y_990;
                    
                    y_990 = 0.5 * uv_988;
                    
                    double res_991;
                    
                    res_991 = x_989 - y_990;
                    *(double*) &mem_331[(i_80 * numY_13 + i_76 * 1) * 8 + 0] =
                        res_991;
                }
            }
            memmove(mem_348 + 0, mem_331 + 0, conc_tmp_449 * numY_13 *
                    sizeof(double));
            for (int i_84 = 0; i_84 < numY_13; i_84++) {
                int x_185;
                
                x_185 = i_84 * 3;
                
                int i_188;
                
                i_188 = x_185 / 3 - (x_185 % 3 != 0 && x_185 % 3 < 0 != 3 <
                                     0 ? 1 : 0);
                
                double y_997;
                
                y_997 = *(double*) &mem_230[(i_188 * 3 + 0 * 1) * 8 + 0];
                
                double y_998;
                
                y_998 = x_566 * y_997;
                
                double y_999;
                
                y_999 = 0.5 * y_998;
                
                double res_1000;
                
                res_1000 = 0.0 - y_999;
                
                int flat_index_196;
                
                flat_index_196 = x_185 + 1;
                
                int i_197;
                
                i_197 = flat_index_196 / 3 - (flat_index_196 % 3 != 0 &&
                                              flat_index_196 % 3 < 0 != 3 <
                                              0 ? 1 : 0);
                
                int y_198;
                
                y_198 = i_197 * 3;
                
                int i_remnant_199;
                
                i_remnant_199 = flat_index_196 - y_198;
                
                double y_1001;
                
                y_1001 = *(double*) &mem_230[(i_197 * 3 + i_remnant_199 * 1) *
                                             8 + 0];
                
                double y_1002;
                
                y_1002 = x_566 * y_1001;
                
                double y_1003;
                
                y_1003 = 0.5 * y_1002;
                
                double res_1004;
                
                res_1004 = res_619 - y_1003;
                
                int flat_index_205;
                
                flat_index_205 = x_185 + 2;
                
                int i_206;
                
                i_206 = flat_index_205 / 3 - (flat_index_205 % 3 != 0 &&
                                              flat_index_205 % 3 < 0 != 3 <
                                              0 ? 1 : 0);
                
                int y_207;
                
                y_207 = i_206 * 3;
                
                int i_remnant_208;
                
                i_remnant_208 = flat_index_205 - y_207;
                
                double y_1005;
                
                y_1005 = *(double*) &mem_230[(i_206 * 3 + i_remnant_208 * 1) *
                                             8 + 0];
                
                double y_1006;
                
                y_1006 = x_566 * y_1005;
                
                double y_1007;
                
                y_1007 = 0.5 * y_1006;
                
                double res_1008;
                
                res_1008 = 0.0 - y_1007;
                *(double*) &mem_351[i_84 * 1 * 8 + 0] = res_1000;
                *(double*) &mem_353[i_84 * 1 * 8 + 0] = res_1004;
                *(double*) &mem_355[i_84 * 1 * 8 + 0] = res_1008;
            }
            for (int i_94 = 0; i_94 < numX_12; i_94++) {
                for (int i_240 = 0; i_240 < numY_13; i_240++) {
                    *(double*) &mem_374[8 * i_240] = *(double*) &mem_348[8 *
                                                                         (numY_13 *
                                                                          i_94) +
                                                                         8 *
                                                                         i_240];
                }
                memmove(mem_376 + 0, mem_353 + 0, numY_13 * sizeof(double));
                for (int i_1016 = 0; i_1016 < i_503; i_1016++) {
                    int res_1017;
                    
                    res_1017 = i_1016 + 1;
                    
                    double x_1018;
                    
                    x_1018 = *(double*) &mem_351[res_1017 * 1 * 8 + 0];
                    
                    int i_1019;
                    
                    i_1019 = res_1017 - 1;
                    
                    double y_1020;
                    
                    y_1020 = *(double*) &mem_376[i_1019 * 1 * 8 + 0];
                    
                    double res_1021;
                    
                    res_1021 = x_1018 / y_1020;
                    
                    double x_1022;
                    
                    x_1022 = *(double*) &mem_376[res_1017 * 1 * 8 + 0];
                    
                    double y_1023;
                    
                    y_1023 = *(double*) &mem_355[i_1019 * 1 * 8 + 0];
                    
                    double y_1024;
                    
                    y_1024 = res_1021 * y_1023;
                    
                    double lw_val_1025;
                    
                    lw_val_1025 = x_1022 - y_1024;
                    *(double*) &mem_376[res_1017 * 1 * 8 + 0] = lw_val_1025;
                    
                    double x_1027;
                    
                    x_1027 = *(double*) &mem_374[res_1017 * 1 * 8 + 0];
                    
                    double y_1028;
                    
                    y_1028 = *(double*) &mem_374[i_1019 * 1 * 8 + 0];
                    
                    double y_1029;
                    
                    y_1029 = res_1021 * y_1028;
                    
                    double lw_val_1030;
                    
                    lw_val_1030 = x_1027 - y_1029;
                    *(double*) &mem_374[res_1017 * 1 * 8 + 0] = lw_val_1030;
                }
                
                double x_1034;
                
                x_1034 = *(double*) &mem_374[i_503 * 1 * 8 + 0];
                
                double y_1035;
                
                y_1035 = *(double*) &mem_376[i_503 * 1 * 8 + 0];
                
                double lw_val_1036;
                
                lw_val_1036 = x_1034 / y_1035;
                *(double*) &mem_374[i_503 * 1 * 8 + 0] = lw_val_1036;
                for (int i_245 = 0; i_245 < numY_13; i_245++) {
                    *(double*) &mem_370[8 * (numY_13 * i_94) + 8 * i_245] =
                        *(double*) &mem_374[8 * i_245];
                }
                for (int j_1039 = 0; j_1039 < i_503; j_1039++) {
                    int res_1040;
                    
                    res_1040 = n_464 - j_1039;
                    
                    double x_1041;
                    
                    x_1041 = *(double*) &mem_370[(i_94 * numY_13 + res_1040 *
                                                  1) * 8 + 0];
                    
                    double x_1042;
                    
                    x_1042 = *(double*) &mem_355[res_1040 * 1 * 8 + 0];
                    
                    int i_1043;
                    
                    i_1043 = res_1040 + 1;
                    
                    double y_1044;
                    
                    y_1044 = *(double*) &mem_370[(i_94 * numY_13 + i_1043 * 1) *
                                                 8 + 0];
                    
                    double y_1045;
                    
                    y_1045 = x_1042 * y_1044;
                    
                    double x_1046;
                    
                    x_1046 = x_1041 - y_1045;
                    
                    double y_1047;
                    
                    y_1047 = *(double*) &mem_376[res_1040 * 1 * 8 + 0];
                    
                    double lw_val_1048;
                    
                    lw_val_1048 = x_1046 / y_1047;
                    *(double*) &mem_370[(i_94 * numY_13 + res_1040 * 1) * 8 +
                                        0] = lw_val_1048;
                }
            }
            for (int i_250 = 0; i_250 < numY_13; i_250++) {
                for (int i_251 = 0; i_251 < numX_12; i_251++) {
                    *(double*) &mem_395[8 * (numX_12 * i_250) + 8 * i_251] =
                        *(double*) &mem_370[8 * (numY_13 * i_251) + 8 * i_250];
                }
            }
            memmove(double_buffer_mem_41 + 0, mem_395 + 0, numY_13 * numX_12 *
                    sizeof(double));
            
            unsigned char* myResult_mem_tmp_175;
            
            myResult_mem_tmp_175 = double_buffer_mem_41;
            myResult_mem_246 = myResult_mem_tmp_175;
        }
        myResult_mem_398 = myResult_mem_246;
        
        double res_1054;
        
        res_1054 = *(double*) &myResult_mem_398[(res_116 * numX_12 + res_115 *
                                                 1) * 8 + 0];
        *(double*) &mem_232[i_136 * 1 * 8 + 0] = res_1054;
    }
    out_mem_77 = mem_232;
    out_memsize_78 = bytes_231;
    out_arrsize_98 = outer_loop_count_8;
    
    struct tuple_int_mem_int retval_0;
    
    retval_0.elem_0 = out_memsize_78;
    retval_0.elem_1 = out_mem_77;
    retval_0.elem_2 = out_arrsize_98;
    return retval_0;
}
int main(int argc, char** argv)
{
    struct timeval t_start,  t_end,  t_diff;
    unsigned long elapsed_usec;
    
    {
        int out_memsize_78;
        unsigned char* out_mem_77;
        int out_arrsize_98;
        int outer_loop_count_8;
        int numX_12;
        int numY_13;
        int numT_14;
        double s0_15;
        double strike_17;
        double t_25;
        double alpha_26;
        double nu_31;
        double beta_32;
        struct tuple_int_mem_int main_ret_1;
        
        if (read_int(&outer_loop_count_8) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&numX_12) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&numY_13) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_int(&numT_14) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "int");
            exit(1);
        }
        if (read_double(&s0_15) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        if (read_double(&strike_17) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        if (read_double(&t_25) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        if (read_double(&alpha_26) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        if (read_double(&nu_31) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        if (read_double(&beta_32) != 0) {
            fprintf(stderr, "Syntax error when reading %s.\n", "float64");
            exit(1);
        }
        gettimeofday(&t_start, NULL);
        main_ret_1 = futhark_main(outer_loop_count_8, numX_12, numY_13, numT_14,
                                  s0_15, strike_17, t_25, alpha_26, nu_31,
                                  beta_32);
        gettimeofday(&t_end, NULL);
        out_memsize_78 = main_ret_1.elem_0;
        out_mem_77 = main_ret_1.elem_1;
        out_arrsize_98 = main_ret_1.elem_2;
        {
            if (out_arrsize_98 == 0) {
                printf("empty(%s)", "float64");
            } else {
                int print_i_2;
                
                putchar('[');
                for (print_i_2 = 0; print_i_2 < out_arrsize_98; print_i_2++) {
                    double* print_elem_3 = (double*) out_mem_77 + print_i_2 * 1;
                    
                    printf("%.6f", *print_elem_3);
                    if (print_i_2 != out_arrsize_98 - 1) {
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
