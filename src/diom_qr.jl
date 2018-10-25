export diom_qr

function diom_qr{T <: Number}(A :: AbstractLinearOperator, b :: AbstractVector{T};
                              M :: AbstractLinearOperator=opEye(size(A,1)),
                              atol :: Float64=1.0e-8, rtol :: Float64=1.0e-6,
                              itmax :: Int=0, memory :: Int=20, verbose :: Bool=false)

  m, n = size(A)
  m == n || error("System must be square")
  size(b, 1) == m || error("Inconsistent problem size")
  verbose && @printf("DQGMRES: system of size %d\n", n)

  # Initial solution x₀ and residual r₀.
  x = zeros(T, n)
  r = copy(b)
  # Compute β.
  rNorm = @knrm2(n, r) # rNorm = ‖r₀‖₂
  rNorm ≈ 0.0 && return x, SimpleStats(true, false, [rNorm], [], "x = 0 is a zero-residual solution")

  iter = 0
  itmax == 0 && (itmax = 2*n)

  rNorms = [rNorm;]
  ε = atol + rtol * rNorm
  verbose && @printf("%5d  %7.1e\n", iter, rNorm)

  # Set up workspace.
  mem = min(memory, itmax) # Memory.
  V = zeros(n, mem) # Preconditioned Krylov vectors, orthogonal basis for {r₀, AM⁻¹r₀, (AM⁻¹)²r₀, ..., (AM⁻¹)ᵐ⁻¹r₀}.
  W = zeros(n, mem) # Directions for x : Wₘ = Vₘ(Rₘ)⁻¹.
  s = zeros(mem)    # Last mem Givens sines used for the factorization QₘRₘ = Hₘ.
  c = zeros(mem)    # Last mem Givens cosines used for the factorization QₘRₘ = Hₘ.
  H = zeros(mem+2)  # Last column of the band hessenberg matrix Hₘ.
  # Each column has at most mem + 1 nonzero elements. hᵢ.ₘ is stored as H[m-i+2].
  # m-i+2 represents the indice of the diagonal where hᵢ.ₘ is located.
  # In addition of that, the last column of Rₘ is also stored in H.

  # Initial γ₁ and V₁.
  γₘ     = rNorm
  V[:,1] = r / rNorm

  # The following stopping criterion compensates for the lag in the
  # residual, but usually increases the number of iterations.
  # solved = sqrt(max(1, iter-mem+1)) * |γₘ₊₁| ≤ ε
  solved = rNorm ≤ ε # less accurate, but acceptable.
  tired = iter ≥ itmax
  status = "unknown"

  while !(solved || tired)
    
    # Update iteration index.
    iter = iter + 1 
    γₘ₋₁ = γₘ

    # Set position in circulars stacks.
    last_pos = mod(iter-2, mem) + 1 # Position corresponding to wₘ₋₁ and vₘ₋₁ in circular stacks W and V.
    pos = mod(iter-1, mem) + 1 # Position corresponding to wₘ and vₘ in circular stacks W and V.
    next_pos = mod(iter, mem) + 1 # Position corresponding to wₘ₊₁ and vₘ₊₁ in circular stacks W and V.

    # Incomplete Arnoldi procedure.
    z = M * V[:,pos] # Forms pₘ
    w = A * z # Forms vₘ₊₁
    for i = max(1, iter-mem+1) : iter
      ipos = mod(i-1, mem) + 1 # Position corresponding to vᵢ in the circular stack V.
      diag = iter - i + 2
      H[diag] = @kdot(n, w, V[:,ipos]) # hᵢ.ₘ = < A * vₘ , vᵢ >
      @kaxpy!(n, -H[diag], V[:,ipos], w) # w ← w - hᵢ.ₘ * vᵢ
    end
    # Compute hₘ₊₁.ₘ and vₘ₊₁.
    H[1] = @knrm2(n, w) # hₘ₊₁.ₘ = ‖vₘ₊₁‖₂
    if H[1] ≉ 0.0 # hₘ₊₁.ₘ ≈ 0 ⇒ "lucky breakdown"
      V[:,next_pos] = w / H[1] # vₘ₊₁ = w / hₘ₊₁.ₘ
    end
    # rₘ₋ₘₑₘ.ₘ ≠ 0 when m ≥ mem + 1
    if iter ≥ mem + 2
      H[mem+2] = 0.0 # hₘ₋ₘₑₘ.ₘ = 0
    end

    # Update the QR factorization of H.
    # Apply mem previous Givens reflections Ωᵢ.
    if iter ≥ 2
    	for i = max(2,iter-mem+1) : iter
      	irot_pos = mod(i-1, mem) + 1 # Position corresponding to cᵢ and sᵢ in circular stacks c and s.
      	diag = iter - i + 2
      	next_diag = diag + 1
      	H_aux        = c[irot_pos] * H[next_diag] + s[irot_pos] * H[diag]
      	H[diag]      = s[irot_pos] * H[next_diag] - c[irot_pos] * H[diag]
      	H[next_diag] = H_aux
    	end

    	# Update γ
    	γₘ   = s[pos] * γₘ₋₁
    	γₘ₋₁ = c[pos] * γₘ₋₁
    end
    
    # Compute next Givens reflection Ωₘ.
    # [cₘ  sₘ] [ ȓₘ.ₘ ] = [rₘ.ₘ]
    # [sₘ -cₘ] [hₘ₊₁.ₘ]   [ 0  ]
    (c[next_pos], s[next_pos], r) = sym_givens(H[2], H[1])

    # Compute the direction ŵₘ, the last column of Wₘ = Vₘ(Rₘ)⁻¹.
    for i = max(1,iter-mem) : iter - 1
      ipos = mod(i-1, mem) + 1 # Position corresponding to wᵢ in the circular stack W.
      diag = iter - i + 2
      # z ← z - rᵢ.ₘ * wᵢ
      @kaxpy!(n, -H[diag], W[:,ipos], z)
    end
    # ŵₘ = z / ȓₘ.ₘ
    W[:,pos] = z / H[2]

    # Update residual norm.
    # ‖ b - Axₘ ‖₂ = hₘ₊₁.ₘ * |γₘ / ȓₘ.ₘ|
    rNorm = H[1] * abs(γₘ / H[2])
    push!(rNorms, rNorm)

    if iter ≥ 2
      # xₜₘₚ = xₜₘₚ + γₘ₋₁ * wₘ₋₁
      @kaxpy!(n, γₘ₋₁, W[:,last_pos], x)
    end

    # Update stopping criterion.
    solved = rNorm ≤ ε
    tired = iter ≥ itmax
    
    if solved || tired
    	# xₘ = xₜₘₚ + (γₘ)_barre * ŵₘ
    	@kaxpy!(n, γₘ, W[:,pos], x)
    end
    # wₘ = ŵₘ * ȓₘ.ₘ / rₘ.ₘ
    W[:,pos] = W[:,pos] * H[2] / r

    verbose && @printf("%5d  %7.1e\n", iter, rNorm)
  end
  verbose && @printf("\n")

  status = tired ? "maximum number of iterations exceeded" : "solution good enough given atol and rtol"
  stats = SimpleStats(solved, false, rNorms, T[], status)
  return (x, stats)
end