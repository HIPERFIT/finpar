fun {int,int,[real],[real],[real],[[real]],[[real]],[[real]],[[real]]} initGrid
    (real s0, real alpha, real nu, real t, int numX, int numY, int numT) =
    let logAlpha = log(alpha) in
    let myTimeline = map(fn real (int i) => t * toReal(i) / (toReal(numT) - 1.0), iota(numT)) in
    let {stdX, stdY} = {20.0 * alpha * s0 * sqrt(t),
                        10.0 * nu         * sqrt(t)} in
    let {dx, dy} = {stdX / toReal(numX), stdY / toReal(numY)} in
    let {myXindex, myYindex} = {trunc(s0 / dx), numY / 2} in
    let myX = map(fn real (int i) => toReal(i) * dx - toReal(myXindex) * dx + s0, iota(numX)) in
    let myY = map(fn real (int i) => toReal(i) * dy - toReal(myYindex) * dy + logAlpha, iota(numY)) in
    let xXy = replicate(numX,replicate(numY, 0.0)) in
    let {myMuX, myVarX, myMuY, myVarY} = {xXy, xXy, xXy, xXy} in
    {myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY}

fun {[[real]],[[real]]} initOperator([real] x) =
    let n = size(0, x) in
    let dxu = x[1] - x[0] in
    let dxl = 0.0 in
    let Dxlow = [[0.0, -1.0 / dxu, 1.0 / dxu]] in
    let Dxxlow = [[0.0, 0.0, 0.0]] in
    let Dxmid = map(fn [real] (int i) => let dxl = x[i] - x[i-1] in
                                         let dxu = x[i+1] - x[i] in
                                         [-dxu/dxl/(dxl+dxu),
                                         (dxu/dxl - dxl/dxu)/(dxl+dxu),
                                         dxl/dxu/(dxl+dxu)],
                   map (op + (1), iota(n-2))) in
    let Dxxmid = map(fn [real] (int i) => let dxl = x[i] - x[i-1] in
                                          let dxu = x[i+1] - x[i] in
                                          [2.0/dxl/(dxl+dxu),
                                          -2.0*(1.0/dxl + 1.0/dxu)/(dxl+dxu),
                                          2.0/dxu/(dxl+dxu)],
                    map (op + (1), iota(n-2))) in
    let dxl = x[n-1] - x[n-2] in
    let dxu = 0.0 in
    let Dxhigh = [[-1.0 / dxl, 1.0 / dxl, 0.0 ]] in
    let Dxxhigh = [[0.0, 0.0, 0.0 ]] in
    let Dx = concat(concat(Dxlow, Dxmid), Dxhigh) in
    let Dxx = concat(concat(Dxxlow, Dxxmid), Dxxhigh) in
    {Dx, Dxx}

fun real max(real x, real y) = if y < x then x else y
fun int maxInt(int x, int y) = if y < x then x else y

fun *[[real]] setPayoff(real strike, [real] myX, [real] myY) =
    let n = size(0, myY) in
    copy(map(fn [real] (real xi) => replicate(n, max(xi-strike,0.0)), myX))

// Returns new myMuX, myVarX, myMuY, myVarY.
fun {[[real]] , [[real]] , [[real]] , [[real]]} updateParams
    ([real] myX, [real] myY, [real] myTimeline, int g, real alpha, real beta, real nu) =
    unzip (map(fn {[real],[real],[real],[real]} (real xi) =>
           unzip (map (fn {real,real,real,real} (real yj) =>
                  {0.0,
                   exp(2.0*(beta*log(xi) + yj - 0.5*nu*nu*myTimeline[g])),
                   0.0,
                   nu * nu}, myY)), myX))

fun {[real],[real]} tridag
    ([real] a, [real] b, [real] c, [real] r, int n) =
    let bet = 1.0/b[0] in
    let {u, uu} = {copy(replicate(n,0.0)),
                   copy(replicate(n,0.0))} in
    let u[0] = r[0] * bet in
    loop ({u, uu, bet}) =
      for j < n-1 do
        let j = j + 1 in
        let uu[j] = c[j-1]*bet in
        let bet = 1.0/(b[j] - a[j] * uu[j]) in
        let u[j] = (r[j] - a[j]*u[j-1]) * bet in
        {u, uu, bet} in
    loop (u) = for j < n - 1 do
                 let j = n - 2 - j in
                 let u[j] = u[j] - uu[j+1]*u[j+1] in
                 u
    in {u, uu}

fun *[[real]] explicitX
    (int numX, int numY, real dtInv,
     [[real]] myResult, [[real]] myMuX, [[real]] myDx, [[real]] myDxx, [[real]] myVarX) =
    copy(map(fn [real] (int j) =>
               map(fn real (int i) =>
                     let kl = if i == 0 then 1 else 0 in
                     let ku = 2 - if i==numX-1 then 1 else 0 in
                     reduce(op +, dtInv*myResult[i,j],
                            map(fn real (int k) =>
                                  let k = k + kl in
                                  0.5 * (myMuX[i,j]*myDx[i,k]+0.5*myVarX[i,j]*myDxx[i,k])
                                * myResult[i+k-1,j],
                                iota(ku-kl+1))),
                   iota(numX)),
               iota(numY)))

fun *[[real]] explicitY
    (int numX, int numY, real dtInv,
     [[real]] myResult, [[real]] myMuY, [[real]] myDy, [[real]] myDyy, [[real]] myVarY) =
    copy(map(fn [real] (int i) =>
               map(fn real (int j) =>
                     let ll = 1 * (if j == 0 then 1 else 0) in
                     let lu = 2 - 1*(if j==numY-1 then 1 else 0) in
                     reduce(op +, 0.0,
                            map(fn real (int l) =>
                                  let l = l + ll in
                                  (myMuY[i,j]*myDy[j,l]+0.5*myVarY[i,j]*myDyy[j,l])
                                * myResult[i,j+l-1],
                                iota(lu-ll+1))),
                   iota(numY)),
               iota(numX)))

fun *[[real]] rollback
    ([real] myX, [real] myY, [real] myTimeline, *[[real]] myResult,
     [[real]] myMuX, [[real]] myDx, [[real]] myDxx, [[real]] myVarX,
     [[real]] myMuY, [[real]] myDy, [[real]] myDyy, [[real]] myVarY, int g) =
    let {numX, numY} = {size(0, myX), size(0, myY)} in
    let numZ = maxInt(numX, numY) in
    let dtInv = 1.0/(myTimeline[g+1]-myTimeline[g]) in
    let u = explicitX(numX, numY, dtInv, myResult, myMuX, myDx, myDxx, myVarX) in
    let v = explicitY(numX, numY, 0.0, myResult, myMuY, myDy, myDyy, myVarY) in
    let u = map(fn [real] ([real] us, [real] vs) => map(op +, zip(us, vs)),
                zip(u, transpose(v))) in
    let u = map(fn [real] ({[real], int} t) =>
                let {uj, j} = t in
                let {a,b,c} = unzip(map(fn {real,real,real} (int i) =>
                                        {-0.5*(myMuX[i,j]*myDx[i,0] + 0.5*myVarX[i,j]*myDxx[i,0]),
                                         dtInv - 0.5*(myMuX[i,j]*myDx[i,1] + 0.5*myVarX[i,j]*myDxx[i,1]),
                                         -0.5*(myMuX[i,j]*myDx[i,2]+0.5*myVarX[i,j]*myDxx[i,2])},
                                         iota(numX))) in
                let {uj, yy} = tridag(a,b,c,uj,numX) in uj,
            zip(u, iota(numY))) in
    loop (myResult) =
      for i < numX do
       let {a,b,c} = unzip(map(fn {real,real,real} (int j) =>
                                 {-0.5*(myMuY[i,j]*myDy[j,0]+0.5*myVarY[i,j]*myDyy[j,0]),
                                  dtInv - 0.5*(myMuY[i,j]*myDy[j,1]+0.5*myVarY[i,j]*myDyy[j,1]),
                                  -0.5*(myMuY[i,j]*myDy[j,2]+0.5*myVarY[i,j]*myDyy[j,2])},
                                  iota(numY))) in
       let y = copy(replicate(numY, 0.0)) in
       loop (y) = for j < numY do
                    let y[j] = dtInv * u[j,i] - 0.5*v[i,j] in
                    y in
       let {ri, yy} = tridag(a,b,c,y,numY) in
       let myResult[i] = ri in myResult
    in myResult

fun real value(int numX, int numY, int numT, real s0, real strike, real t, real alpha, real nu, real beta) =
    let {myXindex, myYindex, myX, myY, myTimeline, myMuX, myVarX, myMuY, myVarY} =
        initGrid(s0, alpha, nu, t, numX, numY, numT) in
    let {myDx, myDxx} = initOperator(myX) in
    let {myDy, myDyy} = initOperator(myY) in
    let myResult = setPayoff(strike, myX, myY) in
    //loop ((myResult, myMuX, myVarX, myMuY, myVarY)) =
    loop (myResult) =
        for i < numT - 1 do
            let i = numT-2-i in
            let {myMuX, myVarX, myMuY, myVarY} =
                updateParams(myX, myY, myTimeline, i, alpha, beta, nu) in
            let myResult = rollback(myX, myY, myTimeline, myResult,
                                    myMuX, myDx, myDxx, myVarX,
                                    myMuY, myDy, myDyy, myVarY, i) in
            //(myResult, myMuX, myVarX, myMuY, myVarY) in
            myResult in
    myResult[myXindex,myYindex]

fun [real] main (int outer_loop_count, int numX, int numY, int numT) =
    let s0 = 0.03 in
    let strike = 0.03 in
    let t = 5.0 in
    let alpha = 0.2 in
    let nu = 0.6 in
    let beta = 0.5 in
    let strikes = map(fn real (int i) => 0.001*toReal(i), iota(outer_loop_count)) in
    let res = map(fn real (real x) => value(numX, numY, numT, s0, x, t, alpha, nu, beta), strikes) in
    res
