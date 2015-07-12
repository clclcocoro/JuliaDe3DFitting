function planeFit(X)
    n = size(X)[1]
    X′ = zeros(n, 3)
    Xmean = mean(X, 1)
    for i in 1:n
        X′[i, :] = X[i, :] - Xmean
    end
    A = zeros(3, 3)
    B = zeros(3)
    for i in 1:n
        A[1, 1] += X′[i, 1] ^ 2
        A[1, 2] += X′[i, 1] * X′[i, 2]
        A[1, 3] += X′[i, 1]
        A[2, 2] += X′[i, 2] ^ 2
        A[2, 3] += X′[i, 2]
        A[3, 3] += 1.0
        B[1] += X′[i, 1] * X′[i, 3]
        B[2] += X′[i, 2] * X′[i, 3]
        B[3] += X′[i, 3]
    end
    A[2, 1] = A[1, 2]
    A[3, 1] = A[1, 3]
    A[3, 2] = A[2, 3]
    
    Y = \(A, B)
    Y[3] = Y[3] - Y[1] * Xmean[1] - Y[2] * Xmean[2] + Xmean[3]
    return Y, Xmean
end


function sphereFit(X, iternum)
    n = size(X)[1]
    Xmean = mean(X, 1)

    function update(a, b, c, Xmean)
        L̄ = 0.0
        L̄a = 0.0
        L̄b = 0.0
        L̄c = 0.0
        function L(x, y, z, a, b, c)
            return sqrt((x - a)^2 + (y - b)^2 + (z - c)^2)
        end
        for i in 1:n
            L̄ += L(X[i, 1], X[i, 2], X[i, 3], a, b, c)
            L̄a += (a - X[i, 1]) / L(X[i, 1], X[i, 2], X[i, 3], a, b, c)
            L̄b += (b - X[i, 2]) / L(X[i, 1], X[i, 2], X[i, 3], a, b, c)
            L̄c += (c - X[i, 3]) / L(X[i, 1], X[i, 2], X[i, 3], a, b, c)
        end
        L̄ /= n
        L̄a /= n
        L̄b /= n
        L̄c /= n

        return Xmean[1] + L̄ * L̄a, Xmean[2] + L̄ * L̄b, Xmean[3] + L̄ * L̄c, L̄
    end

    a = Xmean[1]
    b = Xmean[2]
    c = Xmean[3]
    r = 0.0
    for i in 1:iternum
        a, b, c, r = update(a, b, c, Xmean)
    end
    return a, b, c, r
end


type param
    θ₀::Float64
    θ₁::Float64
    C::Vector{Float64}
    sampleSize::Int64
end


function cylinderFit(X)

    function G(X, W, rSqr, params)
        S = zeros(3, 3)
        P = eye(3, 3)
        sinθ₀ = sin(params.θ₀)
        cosθ₀ = cos(params.θ₀)
        sinθ₁ = sin(params.θ₁)
        cosθ₁ = cos(params.θ₁)
        W[1] = cosθ₀ * sinθ₁
        W[2] = sinθ₀ * sinθ₁
        W[3] = cosθ₁
        S[1, 1] = 0.0
        S[1, 2] = -W[3]
        S[1, 3] = W[2]
        S[2, 1] = W[3]
        S[2, 2] = 0.0
        S[2, 3] = W[1]
        S[3, 1] = -W[2]
        S[3, 2] = W[1]
        S[3, 3] = 0.0
        P -= W * W'
        A = zeros(3, 3)
        B = zeros(3)
        Y = zeros(params.sampleSize, 3)
        averageSqrLength = 0.0
        sqrLength = zeros(params.sampleSize)
        for i in 1:params.sampleSize
            Y[i, :] = P * vec(X[i, :])
            sqrLength[i] = dot(vec(Y[i, :]), vec(Y[i, :]))
            A += vec(Y[i, :]) * vec(Y[i, :])'
            B += sqrLength[i] * vec(Y[i, :])
            averageSqrLength += sqrLength[i]
        end
        A /= params.sampleSize
        B /= params.sampleSize
        averageSqrLength /= params.sampleSize
        Ahat = - S * A * S
        params.C = (Ahat * B) / trace(Ahat * A)
        error = 0.0
        rSqr = 0.0
        for i in 1:params.sampleSize
            term = sqrLength[i] - averageSqrLength - 2 * dot(vec(Y[i, :]), params.C)
            error += term * term
            diff = params.C - vec(Y[i, :])
            rSqr += dot(diff, diff)
        end
        error /= params.sampleSize
        rSqr /= params.sampleSize
        return  error, W, rSqr, params
    end

    params = param(0, 0, zeros(3), size(X)[1])
    average = mean(X, 1)
    for i in 1:params.sampleSize
        X[i, :] -= average
    end
    jMax = 30 
    iMax = 120
    W = zeros(3)
    C = zeros(3)
    θ₀ = 0.0
    θ₁ = 0.0
    rSqr = 0.0
    minError = Inf
    for i in 1:iMax
        params.θ₀ = 2 * pi * (i / iMax)
        for j in 1:jMax
            params.θ₁ = pi / 2.0 * (j / jMax)
            currentW = zeros(3)
            currentRSqr = 0.0
            error, currentW, currentRSqr, params = G(X, currentW, currentRSqr, params)
            if error < minError
                minError = error
                W =  currentW
                C = params.C
                θ₀ = params.θ₀
                θ₁ = params.θ₁
                rSqr = currentRSqr
            end
        end
    end
    return  W, C + vec(average), sqrt(rSqr), θ₀, θ₁
end
