// reduction with +
REAL to_solve(           bool   iddle,
                         uint   n_schedi,
                __local  uchar* flags,
                __local  REAL * sh_mem,
                
                         REAL   scale,
                         REAL   bbi,
                         REAL   yhat
) {
    REAL accum;
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
//REAL2 tmp2 = rootFinding_Brent(iddle, shape, sh_mem4, scales, uplo.y, uplo.x, 1.0e-4, 1000);
//rootFinding_Brent(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);
inline
REAL2 rootFinding_Brent(         bool   iddle,
                                 uint   n_schedi,
                        __local  uchar* flags,
                        __local  REAL * sh_mem,
                                 REAL   bbi,  //sh_mem4[TH_ID].z;
                                 REAL   scale, 
                                 REAL   a, //lb,
                                 REAL   b, //ub, 
                                 REAL   tol, 
                                 uint   iter_max
                                 // (root, fb) : (REAL, REAL)
) { 
    REAL2 res;

    REAL  fa = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, a);
    REAL  fb = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, b);

    __local uchar* bools = flags + LWG_EG;

    if( bools[TH_ID] == 0 && fa*fb >= 0.0 ) {
        res.x = 0.0;
        res.y = ( a >= 0.0 ) ? INFTY : -INFTY;
        bools[TH_ID] = 1;
    } 

    if( fabs(fa) < fabs(fb) ) { REAL tmp = fa; fa = fb; fb = tmp; }
    
    REAL c = a, fc = fa;
    bool mflag = true;
    REAL d     = 0.0;

    uchar ok = 0; int i = 0;
    do {
        { // reduce with && global termination 
            ok = ( fb == 0.0 || fabs(b-a) < tol || i >= iter_max );

            if( ok ) {
                if( bools[ TH_ID ] == 0 ) {
                    res.x = b;
                    res.y = fb;
                }
                bools[ TH_ID ] = 1;
            }
            bool prev_b = bools[TH_ID];
            bools[TH_ID] = iddle || prev_b;
            barrier(CLK_LOCAL_MEM_FENCE);
            
            ok = reduce_reg_and( bools );

            //barrier(CLK_LOCAL_MEM_FENCE);
            bools[TH_ID]= prev_b;
        }
        i ++;

        // the loop's body 
        if( !ok ) {
            REAL s;
            if( bools[TH_ID] == 0 && (!iddle) ) { 
                if( fa == fc || fb == fc ) {
                    s = b - fb * (b - a) / (fb - fa);
                } else {
                    REAL s1 = (a*fb*fc)/( (fa-fb)*(fa-fc) );
                    REAL s2 = (b*fa*fc)/( (fb-fa)*(fb-fc) );
                    REAL s3 = (c*fa*fb)/( (fc-fa)*(fc-fb) );
                    s = s1 + s2 + s3;
                }

                REAL smb = fabs(s - b), bmc = fabs(b - c), cmd = fabs(c - d), atol = fabs(tol);

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
            REAL fs = to_solve(iddle, n_schedi, flags, sh_mem, scale, bbi, s);

            if( bools[TH_ID] == 0 && (!iddle) ) { 
                // d is assigned for the first time here:
                // it's not used above because mflag is set
                d = c;
                c = b; fc = fb;
            
                if( fa*fs < 0.0 ) { b = s; fb = fs; }
                else              { a = s; fa = fs; }
    
                if( fabs(fa) < fabs(fb) ) { 
                    REAL tmp;
                    tmp = a;   a =  b;  b = tmp;
                    tmp = fa; fa = fb; fb = tmp;
                }
            }
        }
    } while (!ok);

    return res;
}


/**
 * Computes Yhat
 */
REAL exactYhat(          bool   iddle,
                __local  uchar* flags, 
                         uint   n_schedi, 
                         REAL   y0, 
                         REAL   y1, 
                __local  REAL4* sh_mem4,
                         REAL   x 
) {
    REAL  res = 0.0, scales;
    REAL2 uplo;
    __local  uchar* bools = flags + LWG_EG;

    { // the first reduction (+,max)
        bool is_last; short lind;
        __local REAL2* sh_mem2 = (__local  REAL2*) (sh_mem4 + LWG_EG);
        
        if( iddle ) {
            uplo = (REAL2) ( 0.0, -INFTY );
        } else {
            REAL4 tmp4 = sh_mem4[ TH_ID ];
            REAL  baix = tmp4.y * x;

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
        bools[ TH_ID ] = 1;
    }

    { // adjust `up'   
        int  lfind = TH_ID - getIotaInd( (__local short*)flags ) + n_schedi - 1;
        REAL log_s = log(uplo.x);
        REAL tmp   = log_s / sh_mem4[ lfind ].z;  //bbi[n_schedi-1];
        lfind     += 1 - n_schedi;

        if ( tmp <= 0.0 ) {
            uplo.x = tmp;
        } else {
            tmp = log_s / sh_mem4[lfind].z;  // bbi[0];
            if ( 0.0 <= tmp ) uplo.x = tmp;
            else              uplo.x = - INFTY;
        }
    }
return res;
//    const REAL yl = lo - EPS;
//    const REAL yu = up + EPS;
    if        ( (bools[ TH_ID ] == 0) && (y1 <= uplo.y - EPS) ) {
        // yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < EPS

        res = y1 + 1.0;
        bools[ TH_ID ] = 1;  
    } else if ( (bools[ TH_ID ] == 0) && (uplo.x + EPS <= y0) ) {
        // yhat is lower than y0 => phi(h_i(x, yhat)) < EPS)

        res = y0 - 1.0;
        bools[ TH_ID ] = 1;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    {        
        uchar prev_flag = bools[TH_ID];

        //REAL root_lb = std::max( uplo.y - EPS, y0 );
        //REAL root_ub = std::min( uplo.x + EPS, y1 );
        uplo.x = min( uplo.x + EPS, y1 );
        uplo.y = max( uplo.y - EPS, y0 );
        // tmp2 = (root, error); //UINT iter;
        REAL2 tmp2 = rootFinding_Brent(
                            iddle,  n_schedi, flags,  (__local REAL*) (sh_mem4+LWG_EG), 
                            sh_mem4[TH_ID].z, scales, uplo.y, uplo.x, 1.0e-4, 1000
                        );

        if( prev_flag == 0 ) {
            res = ( tmp2.y == -INFTY ) ?  
                        y0 - 1.0 : 
                        ( tmp2.y ==  INFTY ) ? y1 + 1.0 : tmp2.x;
        }
    }

    bools[ TH_ID ] = 0; //flags[ TH_ID + LWG_EG ] = 0;
    return res;
}


