// reduction with +
real_t to_solve(           bool   iddle,
                         uint   n_schedi,
                __local  uchar* flags,
                __local  real_t * sh_mem,
                
                         real_t   scale,
                         real_t   bbi,
                         real_t   yhat
) {
    real_t accum;
    sh_mem[TH_ID] = (iddle) ?  0.0 : scale * exp( - bbi * yhat );
    barrier(CLK_LOCAL_MEM_FENCE);

    segm_reduce_plus( sh_mem, flags ); // ToDo 2

    int last_ind  = TH_ID - getIotaInd( (__local short*)flags );
    if(!iddle) last_ind += n_schedi - 1;

    accum = sh_mem[last_ind];
    barrier(CLK_LOCAL_MEM_FENCE);

    return accum - 1.0;
}


/**
 * root finding by Brent's method.
 */ 
//real2_t tmp2 = rootFinding_Brent(iddle, shape, sh_mem4, scales, uplo.y, uplo.x, 1.0e-4, 1000);
//rootFinding_Brent(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);
inline
real2_t rootFinding_Brent(         bool   iddle,
                                 uint   n_schedi,
                        __local  uchar* flags,
                        __local  real_t * sh_mem,
                                 real_t   bbi,  //sh_mem4[TH_ID].z;
                                 real_t   scale, 
                                 real_t   a, //lb,
                                 real_t   b, //ub, 
                                 real_t   tol, 
                                 uint   iter_max
                                 // (root, fb) : (real_t, real_t)
) { 
    real2_t res;

    real_t  fa = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, a);
    real_t  fb = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, b);

    if( (!iddle) && fa*fb >= 0.0 ) {
        res.x = 0.0;
        res.y = ( a >= 0.0 ) ? INFTY : -INFTY;
        iddle = true;
    } 

    if( fabs(fa) < fabs(fb) ) { real_t tmp = fa; fa = fb; fb = tmp; tmp = a; a = b; b = tmp; }
    
    real_t c = a, fc = fa;
    bool mflag = true;
    real_t d     = 0.0;

    uchar ok = 0; int i = 0;
    do {
        { // reduce with && global termination 
            ok = ( fb == 0.0 || fabs(b-a) < tol || i >= iter_max );

            if( ok ) {
                if( !iddle ) {
                    res.x = b;
                    res.y = fb;
                }
                iddle = true;
            }

            flags[TH_ID+get_local_size(0)] = iddle;
            barrier(CLK_LOCAL_MEM_FENCE);            
            ok = reduce_reg_and( flags + get_local_size(0) );
            //barrier(CLK_LOCAL_MEM_FENCE);
        }
        i ++;

        // the loop's body 
        if( !ok ) {
            real_t s;
            if( !iddle ) { 
                if( fa == fc || fb == fc ) {
                    s = b - fb * (b - a) / (fb - fa);
                } else {
                    real_t s1 = (a*fb*fc)/( (fa-fb)*(fa-fc) );
                    real_t s2 = (b*fa*fc)/( (fb-fa)*(fb-fc) );
                    real_t s3 = (c*fa*fb)/( (fc-fa)*(fc-fb) );
                    s = s1 + s2 + s3;
                }

                bool big_cond = 
                     ( (3.0 * a + b) /4.0 > s || s > b)        ||
                     (  mflag  && fabs(b-c)/2.0 <= fabs(s-b) ) ||
                     ( !mflag  && fabs(c-d)/2.0 <= fabs(s-b) ) ||
                     (  mflag  && fabs(b-c)     <= fabs(tol) ) ||
                     ( !mflag  && fabs(c-d)     <= fabs(tol) )  ;
    
                if ( big_cond ) {
                    mflag = true;
                    s     = (a + b) / 2.0;
                } else {
                    mflag = false;
                }
            }
            real_t fs = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, s);

            if( !iddle ) { 
                // d is assigned for the first time here:
                // it's not used above because mflag is set
                d = c;
                c = b; fc = fb;
            
                if( fa*fs < 0.0 ) { b = s; fb = fs; }
                else              { a = s; fa = fs; }
    
                if( fabs(fa) < fabs(fb) ) { 
                    real_t tmp;
                    tmp = a;   a =  b;  b = tmp;
                    tmp = fa; fa = fb; fb = tmp;
                }
            }
        }
    } while (!ok);

    return res;
}


/*********************************************/
/*********************************************/

inline
real2_t rootBisection(             bool   iddle,
                                 uint   n_schedi,
                        __local  uchar* flags,
                        __local  real_t * sh_mem,
                                 real_t   bbi,  //sh_mem4[TH_ID].z;
                                 real_t   scale, 
                                 real_t   a, //lb,
                                 real_t   b, //ub, 
                                 real_t   tol, 
                                 uint   iter_max
                                 // (root, fb) : (real_t, real_t)
) { 
    real2_t res;

    real_t  fa = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, a);
    real_t  fb = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, b);

    if( (!iddle) && fa*fb >= 0.0 ) {
        res.x = 0.0;
        res.y = ( a >= 0.0 ) ? INFTY : -INFTY;
        iddle = true;
    } 
    
    real_t x = b, fx = fb;

    uchar ok = 0; int i = 0;
    do {
        { // reduce with && global termination 
            ok = ( fx == 0.0 || fabs(b-a) < tol || i >= iter_max );

            if( ok ) {
                if( !iddle ) {
                    res.x = x;
                    res.y = fx;
                }
                iddle = true;
            }

            flags[TH_ID+get_local_size(0)] = iddle;
            barrier(CLK_LOCAL_MEM_FENCE); 
            ok = reduce_reg_and( flags + get_local_size(0) );
            //barrier(CLK_LOCAL_MEM_FENCE);
        }
        i ++;

        // the loop's body 
        if( !ok ) {
            real_t s;

            x=(a+b)/2.0;

            real_t fx = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, x);

            if( !iddle ) { 
                if ( fa*fx > 0.0 ) { a = x; fa = fx; }
                else               { b = x; }
            }
        }
    } while (!ok);

    return res;
}


/********************************************/
/********************************************/


/**
 * Computes Yhat
 */
real_t exactYhat(          bool   iddle,
                __local  uchar* flags, 
                         uint   n_schedi, 
                         real_t   y0, 
                         real_t   y1, 
                __local  real4_t* sh_mem4,
                         real_t   x 
) {
    real_t  res = 0.0, scales;
    real2_t uplo;

    { // the first reduction (+,max)
        bool is_last; short lind;
        __local real2_t* sh_mem2 = (__local  real2_t*) (sh_mem4 + LWG_EG);
        
        if( iddle ) {
            uplo = (real2_t) ( 0.0, -INFTY );
        } else {
            real4_t tmp4 = sh_mem4[ TH_ID ];
            real_t  baix = tmp4.y * x;

            uplo.x = tmp4.x * exp( -baix );
            uplo.y = ( tmp4.w - baix ) / tmp4.z;
        }
        scales = uplo.x;

        barrier(CLK_LOCAL_MEM_FENCE);
        sh_mem2[ TH_ID ] = uplo;  
        barrier(CLK_LOCAL_MEM_FENCE);

        segm_reduce_plusmax( sh_mem2, flags ); // ToDo: Fix this please!

        int last_ind  = TH_ID - getIotaInd( (__local short*) flags );
        if(!iddle) last_ind += n_schedi - 1;
        uplo = sh_mem2[last_ind];      // the reduction result!
        barrier(CLK_LOCAL_MEM_FENCE);  
    }

// CHECKING uplo

    if ( n_schedi == 1 ) { 
        res = uplo.y; 
        iddle = true;
    }

    if(!iddle) { // adjust `up'   
        int  lfind = TH_ID - getIotaInd( (__local short*)flags ) + n_schedi - 1;
        real_t log_s = log(uplo.x);
        real_t tmp   = log_s / sh_mem4[ lfind ].z;  //bbi[n_schedi-1];
        lfind     += 1 - n_schedi;

        if ( tmp <= 0.0 ) {
            uplo.x = tmp;
        } else {
            tmp = log_s / sh_mem4[lfind].z;  // bbi[0];
            if ( 0.0 <= tmp ) uplo.x = tmp;
            else              uplo.x = - INFTY;
        }
    }

//    const real_t yl = lo - EPS;
//    const real_t yu = up + EPS;
    if        ( (!iddle) && (y1 <= uplo.y - EPS) ) {
        // yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < EPS

        res = y1 + 1.0;
        iddle = true;  
    } else if ( (!iddle) && (uplo.x + EPS <= y0) ) {
        // yhat is lower than y0 => phi(h_i(x, yhat)) < EPS)

        res = y0 - 1.0;
        iddle = true;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    {
        uplo.x = minR( uplo.x + EPS, y1 );
        uplo.y = maxR( uplo.y - EPS, y0 );
        // tmp2 = (root, error); //UINT iter;

#if 1
        real2_t tmp2 = rootFinding_Brent(
                            iddle,  n_schedi, flags,  (__local real_t*) (sh_mem4+LWG_EG), 
                            sh_mem4[TH_ID].z, scales, uplo.y, uplo.x, 1.0e-4, 1000
                        );
#else

        real2_t tmp2 = rootBisection(
                            iddle,  n_schedi, flags,  (__local real_t*) (sh_mem4+LWG_EG), 
                            sh_mem4[TH_ID].z, scales, uplo.y, uplo.x, 1.0e-4, 1000
                        );

        
#endif

        if( !iddle ) {
            res = ( tmp2.y == -INFTY ) ?  
                        y0 - 1.0 : 
                        ( tmp2.y ==  INFTY ) ? y1 + 1.0 : tmp2.x;
        }
    }

    return res;
}


