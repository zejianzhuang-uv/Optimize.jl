



function brentqv(f::Function, init_x::Vector{Float64}, n::Int64)
    x0 = init_x[1]
    solv = Float64[]
    for x in init_x[2:end]
        xx = f(x0) * f(x)
        if xx < 0
            sol = brentq(f, x0, x)
            if isnan(sol) == false
                push!(solv, sol)
                if length(solv) == n
                    break
                end
                x0 = x
            else
                x0 = x
            end
        else
            continue
        end
    end
    return solv
end

function brentqv(f::Function, init_x::AbstractRange, n::Int64)
    x0 = init_x[1]
    solv = Float64[]
    for x in init_x[2:end]
        xx = f(x0) * f(x)
        if xx < 0
            sol = brentq(f, x0, x)
            if isnan(sol) == false
                push!(solv, sol)
                if length(solv) == n
                    break
                end
                x0 = x
            else
                x0 = x
            end
        else
            continue
        end
    end
    return solv
end




















"""
At the top of the loop the situation is the following:

    1. the root is bracketed between xa and xb
    2. xa is the most recent estimate
    3. xp is the previous estimate
    4. |fp| < |fb|

  The order of xa and xp doesn't matter, but assume xp < xb. Then xa lies to
  the right of xp and the assumption is that xa is increasing towards the root.
  In this situation we will attempt quadratic extrapolation as long as the
  condition

  *  |fa| < |fp| < |fb|

  is satisfied. That is, the function value is decreasing as we go along.
  Note the 4 above implies that the right inequlity already holds.

  The first check is that xa is still to the left of the root. If not, xb is
  replaced by xp and the interval reverses, with xb < xa. In this situation
  we will try linear interpolation. That this has happened is signaled by the
  equality xb == xp;

  The second check is that |fa| < |fb|. If this is not the case, we swap
  xa and xb and resort to bisection.
"""
function brentq(f::Function, xa::Real, xb::Real; iter::Int64=1000, xtol=2e-12, rtol=1e-10, ftol=1e-6)::Float64
    xpre, xcur = xa, xb
    xblk = 0e0
    fpre, fcur, fblk, spre, scur, sbis = 0e0, 0e0, 0e0, 0e0, 0e0, 0e0
    delta, stry, dpre, dblk = 0e0, 0e0, 0e0, 0e0
    fpre = f(xpre)
    fcur = f(xcur)
    
    if abs(fpre) <= ftol
        return xpre
    end
    if abs(fcur) <= ftol
        return xcur
    end
    if signbit(fpre) == signbit(fcur)
        error("f(xa) and f(xb) have the same sign.")
        return 0e0
    end
    i = 1
    while i < iter
        if (fpre != 0) && (fcur != 0) && signbit(fpre) != signbit(fcur)
            xblk = xpre
            fblk = fpre
            spre = scur = xcur - xpre
        end
        if (abs(fblk) < abs(fcur))
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre
        end
        delta = (xtol + rtol*abs(xcur) ) / 2
        sbis = (xblk - xcur) / 2
        if fcur == 0 || abs(sbis) < delta
            # return xcur
            break
        end
        if abs(spre) > delta && abs(fcur) < abs(fpre)
            if (xpre == xblk)
                # interpolate
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                #/* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk*dblk - fpre*dpre) /(dblk*dpre*(fblk - fpre))
            end

            if (2*abs(stry) < minimum([abs(spre), 3*abs(sbis) - delta]))
                spre = scur
                scur = stry
            else
                spre = sbis
                scur = sbis
            end
        else
            spre = sbis
            scur = sbis
        end
        xpre = xcur
        fpre = fcur
        if (abs(scur) > delta)
            xcur += scur
        else
            scur += (sbis > 0 ? delta : -delta)
        end
        fcur = f(xcur)
        i += 1
    end
    if i >= iter
        println("iteration reaches its maximum.")
    end
    xcur = (abs(xcur) < 1e-6 ? 0e0 : xcur)
    if abs(f(xcur) ) > 0.1
        return NaN64
    else
        return xcur
    end
end
